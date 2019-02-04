/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "PositionAndForceTrackingModifier.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellCycleModelOdeHandler.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellId.hpp"

template<unsigned DIM>
PositionAndForceTrackingModifier<DIM>::PositionAndForceTrackingModifier()
: AbstractCellBasedSimulationModifier<DIM>()
  {
  }

template<unsigned DIM>
PositionAndForceTrackingModifier<DIM>::~PositionAndForceTrackingModifier()
{
}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	UpdateCellData(rCellPopulation);

	// Create output file
	OutputFileHandler output_file_handler( outputDirectory + "/", false);
	mpDataFile = output_file_handler.OpenOutputFile("positionsandforces.dat");

	//Initialise method
	CalculateModifierData(rCellPopulation);

	*mpDataFile << "\n";
}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	rCellPopulation.Update();
}

template<unsigned DIM>
std::set<unsigned> PositionAndForceTrackingModifier<DIM>::GetNeighbouringNodeIndices(AbstractCellPopulation<DIM, DIM>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

	assert(!(p_tissue->IsGhostNode(nodeIndex)));

	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	// Find the indices of the elements owned by this node
	std::set<unsigned> containing_elem_indices = p_tissue->GetNode(nodeIndex)->rGetContainingElementIndices();

	// Iterate over these elements
	for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin();
			elem_iter != containing_elem_indices.end();
			++elem_iter)
	{
		// Get all the nodes contained in this element
		// Don't want to include the current node
		unsigned neighbour_global_index;

		for (unsigned local_index=0; local_index<3; local_index++)
		{
			neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);

			if( (neighbour_global_index != nodeIndex) && (!p_tissue->IsGhostNode(neighbour_global_index)) )
			{
				neighbouring_node_indices.insert(neighbour_global_index);
			}
		}
	}

	return neighbouring_node_indices;
}

template<unsigned DIM>
std::vector<unsigned> PositionAndForceTrackingModifier<DIM>::GetEpitheliumInArcLengthOrder(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{

	// Initialise the vector
	std::vector<unsigned> epithelium_in_order;

	// We sort the indices in order of x-coordinate
	std::vector<std::pair<double, unsigned> > epithelial_x_coordinates_and_indices;

	// First, we obtain the epithelial indices and sort them by their x-coordinate. We assume that the left-most cell is the 'starting' cell.
	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{

			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			// Only consider epithelial cells
			if  ( p_type->IsType<DifferentiatedCellProliferativeType>()==false )
			{
				Node<DIM>* p_node = p_cell_population->GetNodeCorrespondingToCell(*cell_iter);

				double x = p_node->rGetLocation()[0]; // Get the x coordinate
				unsigned node_index = p_node->GetIndex(); // Get index

				std::pair<double, unsigned> x_coordinate_and_index = std::make_pair(x, node_index);

				// Update the vector
				epithelial_x_coordinates_and_indices.push_back(x_coordinate_and_index);

			}
		}

		//Sort indices by the x-coordinate
		std::sort(epithelial_x_coordinates_and_indices.begin(), epithelial_x_coordinates_and_indices.end());

		//We now go through the epithelial_indices until we've accounted for every cell
		unsigned num_epithelial_cells = epithelial_x_coordinates_and_indices.size();
		unsigned current_index = epithelial_x_coordinates_and_indices[0].second;
		unsigned cell_count = 1;

		epithelium_in_order.push_back(current_index);

		// The way this method works: we start from the first index, and work through its neighbours until
		// we find another epithelial cell. Note that due to the VT model, we should have at most two neighbouring
		// epithelial cells.
		while (cell_count < num_epithelial_cells)
		{
			std::set<unsigned> neighbour_indices = GetNeighbouringNodeIndices(rCellPopulation, current_index);

			for(std::set<unsigned>::iterator neighbour_iter=neighbour_indices.begin();
					neighbour_iter != neighbour_indices.end();
					++neighbour_iter)
			{

				boost::shared_ptr<AbstractCellProperty> p_type = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType();

				// If the neighbour is an epithelial cell and not in the epithelium index vector.
				if ( (!p_type->IsType<DifferentiatedCellProliferativeType>())&&
						(std::find(epithelium_in_order.begin(), epithelium_in_order.end(), *neighbour_iter) == epithelium_in_order.end()))
				{
					// Update epithelium
					epithelium_in_order.push_back(*neighbour_iter);

					// Update iteration
					current_index = *neighbour_iter;
					cell_count += 1;
				}
			}
		}

	}

	return epithelium_in_order;
}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::CalculateModifierData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Get the epithelium, ordered by arc length.
	std::vector<unsigned> epithelial_indices = GetEpitheliumInArcLengthOrder(rCellPopulation);

	unsigned num_epithelial_cells = epithelial_indices.size();

	// We will record the arc length, x coordinate, y coordinate, horizontal applied force, vertical applied forces,
	// and tangent angles.

	// We're collecting the data in this extremely inefficient	 format for the sake of writing it to the data.
	std::vector<double> x_coordinates;
	std::vector<double> y_coordinates;
	std::vector<double> arclengths;
	std::vector<double> turning_angles;

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		for (unsigned i = 0; i < num_epithelial_cells; i++)
		{
			// Get the location and force via the index
			unsigned epithelial_node_index = epithelial_indices[i];

			c_vector<double, 2> cell_location = p_cell_population->GetNode(epithelial_node_index)->rGetLocation();

			// Add the coordinates
			x_coordinates.push_back(cell_location[0]);
			y_coordinates.push_back(cell_location[1]);

		}

		// Use the position to calculate the
		double current_arclength = 0.0;
		arclengths.push_back(current_arclength);

		// Let's calculate the arc lengths and turning angles for this
		for (unsigned i = 0; i < (num_epithelial_cells - 1); i++)
		{
			// Get the x coordinates
			double x_first = x_coordinates[i];
			double x_second = x_coordinates[i + 1];

			// Get the y_coordinates
			double y_first = y_coordinates[i];
			double y_second =  y_coordinates[i + 1];

			// Calculate the arclength increment and update the current arclength
			double arclength_increment = sqrt(pow(x_second - x_first, 2.0) + pow(y_second - y_first, 2.0));

			current_arclength += arclength_increment;
			arclengths.push_back(current_arclength);

			//Let's calculate the tangent angle as well (so we can work out the tangential and normal stresses)
			c_vector<double, 2> tangent_vector;
			tangent_vector[0] = x_second - x_first;
			tangent_vector[1] = y_second - y_first;

			tangent_vector /= norm_2(tangent_vector); // Normalise vector

			// Calculate the angle and adjust for the quadrant
			double angle = atan(tangent_vector[1]/tangent_vector[0]); //Get initial angle argument
			turning_angles.push_back(angle);

		}

		*mpDataFile << SimulationTime::Instance()->GetTime() << "\n";

		// Write everything to the file now

		// arclength
		for (unsigned i = 0; i < num_epithelial_cells; i++)
		{
			*mpDataFile << arclengths[i] << "\t";
		}
		*mpDataFile << "\n";

		// x coordinates
		for (unsigned i = 0; i < num_epithelial_cells; i++)
		{
			*mpDataFile << x_coordinates[i] << "\t";
		}
		*mpDataFile << "\n";

		// y coordinates
		for (unsigned i = 0; i < num_epithelial_cells; i++)
		{
			*mpDataFile << y_coordinates[i] << "\t";
		}
		*mpDataFile << "\n";

		// vertical forces
		for (unsigned i = 0; i < (num_epithelial_cells - 1); i++)
		{
			*mpDataFile << turning_angles[i] << "\t";
		}
		*mpDataFile << "\n";

	}

}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{

	UpdateCellData(rCellPopulation);

	CalculateModifierData(rCellPopulation);

	// Close output file.
	mpDataFile->close();

}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

	// No parameters to output, so just call method on direct parent class
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PositionAndForceTrackingModifier<1>;
template class PositionAndForceTrackingModifier<2>;
template class PositionAndForceTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PositionAndForceTrackingModifier)
