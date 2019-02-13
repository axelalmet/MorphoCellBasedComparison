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
   mForceMagnitude(DOUBLE_UNSET)
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

//Method overriding the virtual method for AbstractForce. The crux of what really needs to be done.
void VerticalCompressionForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
	double force_magnitude = mForceMagnitude;

	c_vector<double, 2> vertical_force;
	vertical_force[0] = 0.0;
	vertical_force[1] = -1.0*force_magnitude;

	if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation))
	{
		for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{

			MeshBasedCellPopulationWithGhostNodes<2>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&rCellPopulation);

			// Get the node index
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			// Get the cell type
			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			//Apply only to epithelial cells
			if (!p_cell_population->IsGhostNode(node_index))
			{
				if (p_type->template IsType<DifferentiatedCellProliferativeType>()==false)
				{
					rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(vertical_force);
				}
			}


		}
	}

}

void VerticalCompressionForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<ForceMagnitude>"<<  mForceMagnitude << "</ForceMagnitude> \n";

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VerticalCompressionForce)
