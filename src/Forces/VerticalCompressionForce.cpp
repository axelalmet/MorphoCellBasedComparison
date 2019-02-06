#include "VerticalCompressionForce.hpp"
#include "AbstractCellProperty.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"

/*
 * Created on: 21/12/2014
 * Last modified: 02/10/2015
 */

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
VerticalCompressionForce::VerticalCompressionForce()
:  AbstractForce<2>(),
   mForceMagnitude(DOUBLE_UNSET),
   mApplyOnlyToEpithelialCells(false)
   {
   }

VerticalCompressionForce::~VerticalCompressionForce()
{

}

void VerticalCompressionForce::SetForceMagnitude(double forceMagnitude)
{
    mForceMagnitude = forceMagnitude;
}

double VerticalCompressionForce::GetForceMagnitude()
{
    return mForceMagnitude;
}

void VerticalCompressionForce::ApplyForceOnlyToEpithelialCells(bool applyOnlyToEpithelialCells)
{
    mApplyOnlyToEpithelialCells = applyOnlyToEpithelialCells;
}

bool IsForceAppliedOnlyToEpithelialCells()
{
    return mApplyOnlyToEpithelialCells;
}

//Method overriding the virtual method for AbstractForce. The crux of what really needs to be done.
void VerticalCompressionForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
    double force_magnitude = mForceMagnitude;

    c_vector<double, 2> vertical_force;
    vertical_force[0] = 0.0;
    vertical_force[1] = -force_magnitude;

    if (dynamic_cast<NodeBasedCellPopulation<2>*>(pCellPopulation))
    {
        NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);
    }
    else if (dynamic_cast<MeshBasedCellPopulation<2>*>(pCellPopulation))
    {
        MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);
    }
    else if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(pCellPopulation))
    {	
        MeshBasedCellPopulationWithGhostNodes<2>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation);
    }

    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
	{

        // Get the node index
        Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);	// Pointer to node
		unsigned node_index = p_node->GetIndex();

        // Get the cell type
        boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

        bool apply_only_to_epithelial_cells = IsForceAppliedOnlyToEpithelialCells();

        if ( apply_only_to_epithelial_cells )
        {
            if (p_type->IsType<DifferentiatedCellProliferativeType>()==false)
            {
                rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(vertical_force);
            }
        }
        else
        {

            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(vertical_force);

        }

    }

}

void VerticalCompressionForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<ForceMagnitude>"<<  mForceMagnitude << "</ForceMagnitude> \n";
	*rParamsFile <<  "\t\t\t<ApplyForceOnlyToEpithelialCells>" << mApplyForceOnlyToEpithelialCells << "</ApplyForceOnlyToEpithelialCells> \n";

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VerticalCompressionForce)
