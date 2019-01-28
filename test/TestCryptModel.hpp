#ifndef TESTCRYPTMODEL_HPP_
#define TESTCRYPTMODEL_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"

#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CylindricalHoneycombMeshGenerator.hpp" //Generates mesh
#include "NodesOnlyMesh.hpp"
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "LongSpringDivisionModifier.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "LinearSpringForceWithVariableRestLength.hpp"
#include "NonPeriodicBasementMembraneForce.hpp" //BM force (based off Dunn et al. (2012))
#include "OverlappingSpheresBasedBasementMembraneForce.hpp" //Overlapping spheres equivalent
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "BoxModelWithNoDivision";

class TestCryptModel : public AbstractCellBasedTestSuite
{
public:

	void TestCrossSectionGeometry()
	{

	}
};

#endif /* TESTBOXMODELWITHNODIVISION_HPP_ */
