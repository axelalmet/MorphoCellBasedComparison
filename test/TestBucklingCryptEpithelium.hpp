#ifndef TESTSTRESSESINDEFORMATIONFORVTMODELS_HPP_
#define TESTSTRESSESINDEFORMATIONFORVTMODELS_HPP_

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
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "LinearSpringForceWithVariableRestLength.hpp"
#include "NonPeriodicBasementMembraneForce.hpp" //BM force (based off Dunn et al. (2012))
#include "OverlappingSpheresBasedBasementMembraneForce.hpp" //Overlapping spheres equivalent
#include "VerticalCompressionForce.hpp" // Vertical compression force
#include "AnoikisCellKiller.hpp" // Anoikis-based cell killer
#include "PositionAndForceTrackingModifier.hpp" // Modifier to track the epithelium and resultant forces
#include "GeneralisedCellAppliedForceWriter.hpp" // Modifier to write forces
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "StressesInBucklingEpithelium";
static const double M_DT = 0.005;
static const double M_END_TIME = 100.0;
static const double M_SAMPLING_TIMESTEP = M_END_TIME/M_DT;

class TestBucklingCryptEpithelium : public AbstractCellBasedTestSuite
{
public:

	void TestEvolingCryptEpithelium()
	{
		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 45.0;
		double epithelial_stromal_stiffness = 45.0;
		double stromal_stromal_stiffness = 45.0;

		//Set the number of cells across and down for the array
		unsigned cells_across = 30;
		unsigned cells_up = 20;
		unsigned ghosts = 2; //Set the number of ghost node layers

		//Set the basement membrane force parameters
		double bm_stiffness = 12.0;
		double target_curvature = 0.3;
		//

		double epithelial_epithelial_resting_spring_length = 1.0;

		// Generate the periodic mesh
		CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

		//Get the initial non-ghost indices
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		//Get the maximum width and height of the real nodes to define the monolayer
		double max_height = 0.0;
		double max_width = 0.0;

		for (unsigned i = 0; i < location_indices.size(); i++)
		{
			unsigned node_index = location_indices[i];
			double x = p_mesh->GetNode(node_index)->rGetLocation()[0];
			double y = p_mesh->GetNode(node_index)->rGetLocation()[1];

			if (y > max_height)
			{
				max_height = y;
			}

			if (x > max_width)
			{
				max_width = x;
			}
		}


		double left_boundary = 0.4*max_width;
		double right_boundary = 0.6*max_width;



		//Create shared pointers for cell and mutation states
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
		boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();

		//Create tissue of cells. Initially we set them all to be differentiated
		std::vector<CellPtr> cells; //Create vector of cells

		for (unsigned i = 0; i<location_indices.size(); i++)
		{

			unsigned node_index = location_indices[i];

			double y = p_mesh->GetNode(node_index)->rGetLocation()[1];

			if (y == max_height)
			{
				//Set stochastic duration based cell cycle
				UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel(); //Uniformly distributed cell cycle times
				//			NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Uniformly distributed cell cycle times
				p_cycle_model->SetDimension(2);
				double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,13) and set in the past
				p_cycle_model->SetBirthTime(-birth_time);

				CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
				p_cell->InitialiseCellCycleModel(); // For paranoia really.

				p_cell->SetCellProliferativeType(p_stem_type);

				cells.push_back(p_cell);
			}
			else if (y != max_height)
			{
				//Set stochastic duration based cell cycle
				//				UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel(); //Uniformly distributed cell cycle times
				NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Uniformly distributed cell cycle times
				p_cycle_model->SetDimension(2);
				double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,13) and set in the past
				p_cycle_model->SetBirthTime(-birth_time);

				CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
				p_cell->InitialiseCellCycleModel(); // For paranoia really.

				p_cell->SetCellProliferativeType(p_diff_type);

				if (y == 0.0)
				{

					p_cell->AddCellProperty(p_label);

				}

				cells.push_back(p_cell);
			}

		}

		//Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

		//Output data to vtk format so we can visualise it in Paraview
		cell_population.SetWriteVtkAsPoints(true);
//		cell_population.AddPopulationWriter<VoronoiDataWriter>();
		//		cell_population.AddCellWriter<GeneralisedCellAppliedForceWriter>();

		OffLatticeSimulation<2> simulator(cell_population);

		// Set steady state output directory
		std::stringstream out;
		out << "/B_" << bm_stiffness << "_TC_" << target_curvature;
		std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
		simulator.SetOutputDirectory(output_directory);

		simulator.SetDt(M_DT);
		simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
		simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

		// Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
		MAKE_PTR(LinearSpringForceWithVariableRestLength<2>, p_spring_force);
		p_spring_force->SetCutOffLength(epithelial_epithelial_resting_spring_length);
		p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness); //Default is 15
		p_spring_force->SetEpithelialStromalSpringStiffness(epithelial_stromal_stiffness); //Default is 15
		p_spring_force->SetStromalStromalSpringStiffness(stromal_stromal_stiffness); //Default is 15
		p_spring_force->SetEpithelialEpithelialRestingSpringLength(epithelial_epithelial_resting_spring_length);
		simulator.AddForce(p_spring_force);

		//Add basement membrane force
		MAKE_PTR(NonPeriodicBasementMembraneForce, p_bm_force);
		p_bm_force->SetBasementMembraneParameter(bm_stiffness); //Equivalent to beta in SJD's papers
		p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
		p_bm_force->SetLeftCryptBoundary(left_boundary);
		p_bm_force->SetRightCryptBoundary(right_boundary);
		p_bm_force->ApplyForceToCrypt(true);
		//		p_bm_force->ApplyVerticallyDependentTargetCurvature(true);
		simulator.AddForce(p_bm_force);

		//Add anoikis-based cell killer
		MAKE_PTR_ARGS(AnoikisCellKiller, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);
//
//		// Add modifier to track positions
//		MAKE_PTR(PositionAndForceTrackingModifier<2>, p_position_tracking_modifier);
//		simulator.AddSimulationModifier(p_position_tracking_modifier);

		//Fix the bottom row of cells
		c_vector<double, 2> point, normal;

		point(0) = 0.0;
		point(1) = 0.25;
		normal(0) = 0.0;
		normal(1) = -1.0;
		MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc1);

		simulator.Solve();

	}
};

#endif /* TESTBOXMODELWITHNODIVISION_HPP_ */
