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
#include "Debug.hpp"

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
	//	UpdateCellData(rCellPopulation);

	//	CalculateModifierData(rCellPopulation);

}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	CalculateModifierData(rCellPopulation);

}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
	 * We must update CellData in SetupSolve(), otherwise it will not have been
	 * fully initialised by the time we enter the main time loop.
	 */
	//	UpdateCellData(rCellPopulation);

	// Create output file
	OutputFileHandler output_file_handler(outputDirectory + "/", false);
	mpDataFile = output_file_handler.OpenOutputFile("positionsandforces.dat");

	//	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	//		{
	//			MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);
	//
	//			for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
	//					cell_iter != rCellPopulation.End();
	//					++cell_iter)
	//			{
	//
	//				boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();
	//				unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get index
	//
	//				// Only consider epithelial cells
	//				if (!p_tissue->IsGhostNode(node_index))
	//				{
	//					Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
	//					p_node->ClearAppliedForce();
	//				}
	//			}
	//
	//		}

	//Initialise method
	//	CalculateModifierData(rCellPopulation);

	//	*mpDataFile << "\n";
}

template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	rCellPopulation.Update();
}

template<unsigned DIM>
std::set<unsigned> PositionAndForceTrackingModifier<DIM>::GetNeighbouringEpithelialIndices(AbstractCellPopulation<DIM, DIM>& rCellPopulation, unsigned nodeIndex)
{

	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{

		MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		p_tissue->CreateVoronoiTessellation();

		assert(!(p_tissue->IsGhostNode(nodeIndex)));

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
					// Only take epithelial indices
					boost::shared_ptr<AbstractCellProperty> p_type = p_tissue->GetCellUsingLocationIndex(nodeIndex)->GetCellProliferativeType();

					if (!p_type->template IsType<DifferentiatedCellProliferativeType>() )
					{
						neighbouring_node_indices.insert(neighbour_global_index);
					}
				}
			}
		}
	}
	else if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{

		MeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		p_tissue->CreateVoronoiTessellation();

		assert(!(p_tissue->IsGhostNode(nodeIndex)));

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
					// Only take epithelial indices
					boost::shared_ptr<AbstractCellProperty> p_type = p_tissue->GetCellUsingLocationIndex(nodeIndex)->GetCellProliferativeType();

					if (!p_type->template IsType<DifferentiatedCellProliferativeType>() )
					{
						neighbouring_node_indices.insert(neighbour_global_index);
					}
				}
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
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		p_tissue->CreateVoronoiTessellation();

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{

			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get index

			// Only consider epithelial cells
			if (!p_tissue->IsGhostNode(node_index))
			{
				if  ( p_type->template IsType<DifferentiatedCellProliferativeType>()==false )
				{

					double x = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[0]; // Get the x coordinate

					std::pair<double, unsigned> x_coordinate_and_index = std::make_pair(x, node_index);

					// Update the vector
					epithelial_x_coordinates_and_indices.push_back(x_coordinate_and_index);

				}
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
			std::set<unsigned> neighbour_indices = GetNeighbouringEpithelialIndices(rCellPopulation, current_index);

			for(std::set<unsigned>::iterator neighbour_iter=neighbour_indices.begin();
					neighbour_iter != neighbour_indices.end();
					++neighbour_iter)
			{

				if (!p_tissue->IsGhostNode(*neighbour_iter))
				{
					boost::shared_ptr<AbstractCellProperty> p_type = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType();

					// If the neighbour is an epithelial cell and not already in the vector, add it.
					if ( (!p_type->template IsType<DifferentiatedCellProliferativeType>())&&
							(std::find(epithelium_in_order.begin(), epithelium_in_order.end(), *neighbour_iter) == epithelium_in_order.end()))
					{
						// Update epithelium
						epithelium_in_order.push_back(*neighbour_iter);

						// Update iteration
						current_index = *neighbour_iter;
						cell_count += 1;
						break;
					}
				}
			}
		}

	}
	else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		p_tissue->CreateVoronoiTessellation();

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{

			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get index

			// Only consider epithelial cells
			if (!p_tissue->IsGhostNode(node_index))
			{
				if  ( p_type->template IsType<DifferentiatedCellProliferativeType>()==false )
				{

					double x = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[0]; // Get the x coordinate

					std::pair<double, unsigned> x_coordinate_and_index = std::make_pair(x, node_index);

					// Update the vector
					epithelial_x_coordinates_and_indices.push_back(x_coordinate_and_index);

				}
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
			std::set<unsigned> neighbour_indices = GetNeighbouringEpithelialIndices(rCellPopulation, current_index);

			for(std::set<unsigned>::iterator neighbour_iter=neighbour_indices.begin();
					neighbour_iter != neighbour_indices.end();
					++neighbour_iter)
			{

				if (!p_tissue->IsGhostNode(*neighbour_iter))
				{
					boost::shared_ptr<AbstractCellProperty> p_type = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType();

					// If the neighbour is an epithelial cell and not already in the vector, add it.
					if ( (!p_type->template IsType<DifferentiatedCellProliferativeType>())&&
							(std::find(epithelium_in_order.begin(), epithelium_in_order.end(), *neighbour_iter) == epithelium_in_order.end()))
					{
						// Update epithelium
						epithelium_in_order.push_back(*neighbour_iter);

						// Update iteration
						current_index = *neighbour_iter;
						cell_count += 1;
						break;
					}
				}
			}
		}

	}
	else if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
//		NodeBasedCellPopulation<DIM>* p_tissue = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

		for(typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{

			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get index

			// Only consider epithelial cells
			if  ( p_type->template IsType<DifferentiatedCellProliferativeType>()==false )
			{

				double x = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[0]; // Get the x coordinate

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
			std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(current_index);

			for(std::set<unsigned>::iterator neighbour_iter=neighbour_indices.begin();
					neighbour_iter != neighbour_indices.end();
					++neighbour_iter)
			{


				boost::shared_ptr<AbstractCellProperty> p_type = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType();

				// If the neighbour is an epithelial cell and not already in the vector, add it.
				if ( (!p_type->template IsType<DifferentiatedCellProliferativeType>())&&
						(std::find(epithelium_in_order.begin(), epithelium_in_order.end(), *neighbour_iter) == epithelium_in_order.end()))
				{
					// Update epithelium
					epithelium_in_order.push_back(*neighbour_iter);

					// Update iteration
					current_index = *neighbour_iter;
					cell_count += 1;
					break;
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

	// We're collecting the data in this extremely inefficient format for the sake of writing it to the data.
	std::vector<double> x_coordinates;
	std::vector<double> y_coordinates;
	std::vector<double> horizontal_forces;
	std::vector<double> vertical_forces;

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		p_tissue->CreateVoronoiTessellation();

		for (unsigned i = 0; i < num_epithelial_cells; i++)
		{
			// Get the location and force via the index
			unsigned epithelial_node_index = epithelial_indices[i];

			c_vector<double, 2> cell_location = rCellPopulation.GetNode(epithelial_node_index)->rGetLocation();


			c_vector<double, 2> applied_force =  rCellPopulation.GetNode(epithelial_node_index)->rGetAppliedForce();


			// Add the coordinates
			x_coordinates.push_back(cell_location[0]);
			y_coordinates.push_back(cell_location[1]);
			horizontal_forces.push_back(applied_force[0]);
			vertical_forces.push_back(applied_force[1]);

		}

	}
	else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		p_tissue->CreateVoronoiTessellation();

		for (unsigned i = 0; i < num_epithelial_cells; i++)
		{
			// Get the location and force via the index
			unsigned epithelial_node_index = epithelial_indices[i];
			PRINT_VARIABLE(epithelial_node_index);

			c_vector<double, 2> cell_location = rCellPopulation.GetNode(epithelial_node_index)->rGetLocation();

			PRINT_2_VARIABLES(cell_location[0], cell_location[1]);

			std::vector<double> node_attributes = rCellPopulation.GetNode(epithelial_node_index)->rGetNodeAttributes();

			PRINT_2_VARIABLES(cell_location[0], cell_location[1]);


			for (unsigned j = 0; j < node_attributes.size(); j++)
			{
				PRINT_VARIABLE(node_attributes[j]);
			}

			c_vector<double, 2> applied_force =  rCellPopulation.GetNode(epithelial_node_index)->rGetAppliedForce();


			// Add the coordinates
			x_coordinates.push_back(cell_location[0]);
			y_coordinates.push_back(cell_location[1]);
			horizontal_forces.push_back(applied_force[0]);
			vertical_forces.push_back(applied_force[1]);

		}
	}


	*mpDataFile << SimulationTime::Instance()->GetTime() << "\n";

	// Write everything to the file now

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

	// horizontal forces
	for (unsigned i = 0; i < num_epithelial_cells; i++)
	{
		*mpDataFile << horizontal_forces[i] << "\t";
	}
	*mpDataFile << "\n";

	// vertical forces
	for (unsigned i = 0; i < num_epithelial_cells; i++)
	{
		*mpDataFile << vertical_forces[i] << "\t";
	}
	*mpDataFile << "\n";

}


template<unsigned DIM>
void PositionAndForceTrackingModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{

	//	UpdateCellData(rCellPopulation);

	//	CalculateModifierData(rCellPopulation);

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
