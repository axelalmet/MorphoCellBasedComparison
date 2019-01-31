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

#include "LongSpringDivisionModifier.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellCycleModelOdeHandler.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellId.hpp"

template<unsigned DIM>
LongSpringDivisionModifier<DIM>::LongSpringDivisionModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
LongSpringDivisionModifier<DIM>::~LongSpringDivisionModifier()
{
}

template<unsigned DIM>
std::vector<c_vector<unsigned, 2> > LongSpringDivisionModifier<DIM>::GetEpithelialEpithelialPairs(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

    std::vector<c_vector<unsigned, 2> > epithelial_epithelial_pairs; // Initialise vector
    // Make sure the cell population is updated
    rCellPopulation.Update();

    if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
	{
		MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

		p_cell_population->CreateVoronoiTessellation();

        //Iterate over all the springs
		for (typename MeshBasedCellPopulationWithGhostNodes<DIM>::SpringIterator spring_iter = p_cell_population->SpringsBegin();
				spring_iter != p_cell_population->SpringsEnd();
				++spring_iter)
		{
			unsigned nodeA_global_index = spring_iter.GetNodeA()->GetIndex();
			unsigned nodeB_global_index = spring_iter.GetNodeB()->GetIndex();

			//If they're both real nodes (don't really need this, but paranoia demands it)
			if ( (!p_cell_population->IsGhostNode(nodeA_global_index))&&(!p_cell_population->IsGhostNode(nodeB_global_index)) )
			{   

                // Get the types associated with these cells
				boost::shared_ptr<AbstractCellProperty> p_type_A = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index)->GetCellProliferativeType();
				boost::shared_ptr<AbstractCellProperty> p_type_B = p_cell_population->GetCellUsingLocationIndex(nodeB_global_index)->GetCellProliferativeType();

                //If they're both epithelial cell
				if ( (!p_type_A->IsType<DifferentiatedCellProliferativeType>())&&(p_type_A == p_type_B) )
				{
                    c_vector<unsigned, 2> pair;
                    pair(0) = nodeA_global_index;
                    pair(1) = nodeB_global_index;

                    epithelial_epithelial_pairs.push_back(pair);
				}
			}
		}
    }
    else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
		MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

		p_cell_population->CreateVoronoiTessellation();

        //Iterate over all the springs
		for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iter = p_cell_population->SpringsBegin();
				spring_iter != p_cell_population->SpringsEnd();
				++spring_iter)
		{
			unsigned nodeA_global_index = spring_iter.GetNodeA()->GetIndex();
			unsigned nodeB_global_index = spring_iter.GetNodeB()->GetIndex();

            boost::shared_ptr<AbstractCellProperty> p_type_A = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index)->GetCellProliferativeType();
            boost::shared_ptr<AbstractCellProperty> p_type_B = p_cell_population->GetCellUsingLocationIndex(nodeB_global_index)->GetCellProliferativeType();

            //If they're both epithelial cells
            if ( (!p_type_A->IsType<DifferentiatedCellProliferativeType>())&&(p_type_A == p_type_B) )
            {
                c_vector<unsigned, 2> pair;
                pair(0) = nodeA_global_index;
                pair(1) = nodeB_global_index;

                epithelial_epithelial_pairs.push_back(pair);
            }
        }
    }

    return epithelial_epithelial_pairs;
}

template<unsigned DIM>
bool IsTimeToDivide()
{
    return true;
}

template<unsigned DIM>
c_vector<unsigned, 2> LongSpringDivisionModifier<DIM>::PickPairToDivide(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

    //Initialise vector
    std::vector<c_vector<unsigned, 2> > epithelial_epithelial_pairs = GetEpithelialEpithelialPairs(rCellPopulation);

    //Initise vector of pairs of distances and node index pairs
    std::vector<std::pair<double, c_vector<unsigned, 2> > > distances_and_pairs;

    if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
    {
        MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

        p_cell_population->CreateVoronoiTessellation();

        for (unsigned i = 0; i < epithelial_epithelial_pairs.size(); i++)
        {
            c_vector<unsigned, 2> pair = epithelial_epithelial_pairs[i];

            unsigned nodeA_global_index = pair[0];
            unsigned nodeB_global_index = pair[1];

            // Get the distance between the two
            double distance = p_cell_population->rGetMesh().GetDistanceBetweenNodes(nodeA_global_index, nodeB_global_index);

            std::pair<double, c_vector<unsigned, 2> > distance_and_index_pair = std::make_pair(1.0/distance, pair);

           distances_and_pairs.push_back(distance_and_index_pair);

        }
    }
    else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

        p_cell_population->CreateVoronoiTessellation();

        for (unsigned i = 0; i < epithelial_epithelial_pairs.size(); i++)
        {
            c_vector<unsigned, 2> pair = epithelial_epithelial_pairs[i];

            unsigned nodeA_global_index = pair[0];
            unsigned nodeB_global_index = pair[1];

            // Get the distance between the two
            double distance = p_cell_population->rGetMesh().GetDistanceBetweenNodes(nodeA_global_index, nodeB_global_index);

            std::pair<double, c_vector<unsigned, 2> > distance_and_index_pair = std::make_pair(1.0/distance, pair);

            distances_and_pairs.push_back(distance_and_index_pair);

        }
    }

//    // The first pair should be the one with the maximal distance
//    std::sort(distances_and_pairs.begin(), distances_and_pairs.end());

    std::pair<double, c_vector<unsigned, 2> > maximal_distance_pair = distances_and_pairs[0];

    return maximal_distance_pair.second;

}

template<unsigned DIM>
void LongSpringDivisionModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void LongSpringDivisionModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void LongSpringDivisionModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

//    if (IsTimeToDivide() )
//    {
//        // Pick the pair to diivde
//        c_vector<unsigned, 2> dividing_pair = PickPairToDivide(rCellPopulation);
//
//        unsigned nodeA_global_index = dividing_pair[0];
//        unsigned nodeB_global_index = dividing_pair[1];
//
//        // Split the edge between the two nodes
//        if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
//        {
//            MeshBasedCellPopulationWithGhostNodes<DIM>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);
//
//            p_cell_population->CreateVoronoiTessellation();
//
//            // Get the nodes using the indices
//            Node<DIM>* p_nodeA = p_cell_population->GetNode(nodeA_global_index);
//            Node<DIM>* p_nodeB = p_cell_population->GetNode(nodeB_global_index);
//
//            // Split the edge between the two
//            c_vector<unsigned, 3> new_node_indices = p_cell_population->rGetMesh().SplitEdge(p_nodeA, p_nodeB);
//
//            // As a new node has been added, we need to add a cell there.
//            unsigned new_node_index = new_node_indices[0];
//
//            CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index);
//
//            // Create copy of cell property collection to modify for daughter cell
//            CellPropertyCollection daughter_property_collection = p_neighbour_cell->rGetCellPropertyCollection();
//
//            // Remove the CellId from the daughter cell a new one will be assigned in the constructor
//            daughter_property_collection.RemoveProperty<CellId>();
//
//            CellPtr p_new_cell(new Cell(p_neighbour_cell->GetMutationState(),
//                                        p_neighbour_cell->GetCellCycleModel()->CreateCellCycleModel(),
//                                        p_neighbour_cell->GetSrnModel()->CreateSrnModel(),
//                                        false,
//                                        daughter_property_collection));
//
//            // Add new cell to cell population
//            p_cell_population->mCells.push_back(p_new_cell);
//            p_cell_population->AddCellUsingLocationIndex(new_node_index, p_new_cell);
//
//        }
//        else if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
//        {
//            MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);
//
//            p_cell_population->CreateVoronoiTessellation();
//
//            // Get the nodes using the indices
//            Node<DIM>* p_nodeA = p_cell_population->GetNode(nodeA_global_index);
//            Node<DIM>* p_nodeB = p_cell_population->GetNode(nodeB_global_index);
//
//            // Split the edge between the two
//            c_vector<unsigned, 3> new_node_indices = p_cell_population->rGetMesh().SplitEdge(p_nodeA, p_nodeB);
//
//            // As a new node has been added, we need to add a cell there.
//            unsigned new_node_index = new_node_indices[0];
//
//            CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index);
//
//            // Create copy of cell property collection to modify for daughter cell
//            CellPropertyCollection daughter_property_collection = p_neighbour_cell->rGetCellPropertyCollection();
//
//            // Remove the CellId from the daughter cell a new one will be assigned in the constructor
//            daughter_property_collection.RemoveProperty<CellId>();
//
//            CellPtr p_new_cell(new Cell(p_neighbour_cell->GetMutationState(),
//                                        p_neighbour_cell->GetCellCycleModel()->CreateCellCycleModel(),
//                                        p_neighbour_cell->GetSrnModel()->CreateSrnModel(),
//                                        false,
//                                        daughter_property_collection));
//
//            // Add new cell to cell population
//            p_cell_population->mCells.push_back(p_new_cell);
//            p_cell_population->AddCellUsingLocationIndex(new_node_index, p_new_cell);
//
//        }
//
//    }
}

template<unsigned DIM>
void LongSpringDivisionModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class LongSpringDivisionModifier<1>;
template class LongSpringDivisionModifier<2>;
template class LongSpringDivisionModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LongSpringDivisionModifier)
