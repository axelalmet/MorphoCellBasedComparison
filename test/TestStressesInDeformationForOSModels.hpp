#ifndef TESTSTRESSESINDEFORMATIONFOROSMODELS_HPP_
#define TESTSTRESSESINDEFORMATIONFOROSMODELS_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"

#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CylindricalHoneycombMeshGenerator.hpp" //Generates mesh
#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "NodesOnlyMesh.hpp"
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "NodeBasedCellPopulation.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "LinearSpringForceWithVariableRestLength.hpp"
#include "OverlappingSpheresBasedBasementMembraneForce.hpp" //Overlapping spheres equivalent
#include "VerticalCompressionForce.hpp" // Vertical compression force
#include "AnoikisCellKiller.hpp" // Anoikis-based cell killer
#include "PositionAndForceTrackingModifier.hpp" // Modifier to track the epithelium and resultant forces
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "MeasuringOSStresses";
static const double M_DT = 0.005;
static const double M_END_TIME = 0.1;
//static const double M_SECOND_END_TIME = 1.0;
static const double M_SAMPLING_TIMESTEP = M_END_TIME/M_DT;

class TestStressesInDeformationForOSModels : public AbstractCellBasedTestSuite
{
public:
	//	void TestFlatEpithelium()
	//	{
	//
	//		//Set all the spring stiffness variables
	//		double epithelial_epithelial_stiffness = 45.0;
	//		double epithelial_stromal_stiffness = 45.0;
	//		double stromal_stromal_stiffness = 45.0;
	//
	//		//Set the number of cells across and down for the array
	//		unsigned cells_across = 20;
	//		unsigned cells_up = 6;
	//
	//		//Set the basement membrane force parameters
	//		double bm_stiffness = 5.0;
	//		double target_curvature = 0.0;
	//
	//		double epithelial_epithelial_resting_spring_length = 1.0;
	//
	//		double radius_of_interaction = 1.5;
	//		double division_separation = 0.1;
	//
	//		for (double vert = 1.0; vert < 2.0; vert += 1.0)
	//		{
	//			double vertical_force_magnitude = 0.0*vert;
	//
	//			HoneycombMeshGenerator generator(cells_across, cells_up, 0); //Create mesh
	//			MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //Generate mesh
	//
	//			//Create some nodes
	//			std::vector<Node<2>*> nodes;
	//
	//			Cylindrical2dNodesOnlyMesh mesh(20.0);
	//			mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 20.0); //Construct mesh
	//
	//			//Get the real indices
	//			std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
	//
	//			//Create the vector of cells
	//			//Create shared pointers for cell and mutation states
	//			boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
	//			boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
	//			boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
	//
	//			//Create tissue of cells. Initially we set them all to be differentiated
	//			std::vector<CellPtr> cells; //Create vector of cells
	//
	//			for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
	//			{
	//				//Set stochastic duration based cell cycle
	//				NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Don't give them any cell cycle model yet.
	//				p_cycle_model->SetDimension(2);
	//
	//				CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
	//				p_cell->InitialiseCellCycleModel(); // For paranoia really.
	//
	//				p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated
	//
	//				cells.push_back(p_cell);
	//			}
	//
	//			//Create cell population
	//			NodeBasedCellPopulation<2> cell_population(mesh, cells);
	//
	//			//Set the division separation
	//			cell_population.SetMeinekeDivisionSeparation(division_separation);
	//
	//			//Get the maximum width and height of the real nodes to define the monolayer
	//			double max_height = 0.0;
	//			double max_width = 0.0;
	//
	//			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
	//					cell_iter != cell_population.End(); ++cell_iter)
	//			{
	//				double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
	//				double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
	//
	//				if (y > max_height)
	//				{
	//					max_height = y;
	//				}
	//
	//				if (x > max_width)
	//				{
	//					max_width = x;
	//				}
	//			}
	//
	//			double left_boundary = 0.25*max_width;
	//			double right_boundary = 0.75*max_width;
	//
	//			// We now 'un-differentiate' the top row to initialise the epithelium, while
	//			// labelling the bottom row for the fixed boundary condition.
	//			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
	//					cell_iter != cell_population.End(); ++cell_iter)
	//			{
	//
	//				double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];
	//
	//				if (y == 0.0)
	//				{
	//
	//					boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();
	//					cell_iter->AddCellProperty(p_label);
	//
	//				}
	//				else if (y == max_height)
	//				{
	//					cell_iter->SetCellProliferativeType(p_stem_type);
	//
	//				}
	//
	//			}
	//
	//			//Output data to vtk format so we can visualise it in Paraview
	//			//			cell_population.AddCellWriter<GeneralisedCellAppliedForceWriter>();
	//
	//			OffLatticeSimulation<2> simulator(cell_population);
	//
	//			//Set output directory
	//			std::stringstream out;
	//			out << "/FLAT/B_" << bm_stiffness << "_TC_" << target_curvature << "/V_" << vertical_force_magnitude
	//					<< "/L0_" << epithelial_epithelial_resting_spring_length;
	//			std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
	//			simulator.SetOutputDirectory(output_directory);
	//
	//			simulator.SetDt(M_DT);
	//			simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
	//			simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state
	//
	//			//Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
	//			MAKE_PTR(LinearSpringForceWithVariableRestLength<2>, p_spring_force);
	////			p_spring_force->SetCutOffLength(epithelial_epithelial_resting_spring_length);
	//			p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness); //Default is 15
	//			p_spring_force->SetEpithelialStromalSpringStiffness(epithelial_stromal_stiffness); //Default is 15
	//			p_spring_force->SetStromalStromalSpringStiffness(stromal_stromal_stiffness); //Default is 15
	//			p_spring_force->SetEpithelialEpithelialRestingSpringLength(epithelial_epithelial_resting_spring_length);
	//			p_spring_force->SetCutOffLength(radius_of_interaction);
	//			p_spring_force->SetMeinekeDivisionRestingSpringLength(division_separation);
	//			simulator.AddForce(p_spring_force);
	//
	//			//			//Add basement membrane force
	//			MAKE_PTR(OverlappingSpheresBasedBasementMembraneForce, p_bm_force);
	//			p_bm_force->SetBasementMembraneParameter(bm_stiffness); //Equivalent to beta in SJD's papers
	//			p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
	//			p_bm_force->SetCryptGeometry(true); // Apply to crypt
	//			p_bm_force->ApplyPeriodicForce(true); // Apply if we have periodic mesh
	//			p_bm_force->SetLeftCryptBoundary(left_boundary);
	//			p_bm_force->SetRightCryptBoundary(right_boundary);
	//			p_bm_force->SetCutOffRadius(radius_of_interaction); //Set cut off radius for defining neighbours
	//			simulator.AddForce(p_bm_force);
	//
	//			//			Add vertical compression force
	//			MAKE_PTR(VerticalCompressionForce, p_vertical_force);
	//			p_vertical_force->SetForceMagnitude(vertical_force_magnitude);
	//			simulator.AddForce(p_vertical_force);
	//
	//			//Add anoikis-based cell killer
	//			MAKE_PTR_ARGS(AnoikisCellKiller, p_anoikis_killer, (&cell_population));
	//			p_anoikis_killer->SetCutOffRadius(radius_of_interaction);
	//			simulator.AddCellKiller(p_anoikis_killer);
	//
	//			//			// Add modifier to track positions
	//			//			MAKE_PTR(PositionAndForceTrackingModifier<2>, p_position_tracking_modifier);
	//			//			simulator.AddSimulationModifier(p_position_tracking_modifier);
	//
	//			//			Fix the bottom row of cells
	//			c_vector<double, 2> point, normal;
	//
	//			point(0) = 0.0;
	//			point(1) = 0.25;
	//			normal(0) = 0.0;
	//			normal(1) = -1.0;
	//			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
	//			simulator.AddCellPopulationBoundaryCondition(p_bc1);
	//
	//			simulator.Solve();
	//
	//			// Get a pointer to the vertical compression force to unload the force
	//			boost::shared_ptr<VerticalCompressionForce> p_vertical_compression_force =
	//					boost::static_pointer_cast<VerticalCompressionForce>(simulator.rGetForceCollection()[2]);
	//
	//			p_vertical_compression_force->SetForceMagnitude(0.0);
	//
	//			simulator.SetEndTime(M_END_TIME + M_SECOND_END_TIME); // Set new end time
	//
	//			simulator.Solve(); // Re-run simulation
	//
	//			//Tidying up
	//			SimulationTime::Instance()->Destroy();
	//			SimulationTime::Instance()->SetStartTime(0.0);
	//
	//		}
	//
	//	}

	void TestBuckledEpithelium()
	{

		//		Set all the spring stiffness variables
		//				double epithelial_epithelial_stiffness = 45.0;
		//				double epithelial_stromal_stiffness = 45.0;
		//				double stromal_stromal_stiffness = 45.0;

		//Set the number of cells across and down for the array
		unsigned cells_across = 10;
		unsigned cells_up = 20;

		//Set the basement membrane force parameters
		double bm_stiffness = 5.0;
		double target_curvature = 0.0;

		double epithelial_epithelial_resting_spring_length = 1.0;

		double radius_of_interaction = 1.5;
		double division_separation = 0.1;

		for (double vert = 1.0; vert < 2.0; vert += 1.0)
		{
			double vertical_force_magnitude = 0.0*vert;

			HoneycombMeshGenerator generator(cells_across, cells_up, 0); //Create mesh
			MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //Generate mesh

			//Create some nodes
			std::vector<Node<2>*> nodes;

			Cylindrical2dNodesOnlyMesh mesh(20.0);
			mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 20.0); //Construct mesh

			//Get the real indices
			std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

			//Create the vector of cells
			//Create shared pointers for cell and mutation states
			boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
			boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();

			//Create tissue of cells. Initially we set them all to be differentiated
			std::vector<CellPtr> cells; //Create vector of cells

			for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
			{
				//Set stochastic duration based cell cycle
				NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Don't give them any cell cycle model yet.
				p_cycle_model->SetDimension(2);

				CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
				p_cell->InitialiseCellCycleModel(); // For paranoia really.

				p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated

				cells.push_back(p_cell);
			}

			//Create cell population
			NodeBasedCellPopulation<2> cell_population(mesh, cells);

			//Set the division separation
			cell_population.SetMeinekeDivisionSeparation(division_separation);

			//Get the maximum width and height of the real nodes to define the monolayer
			double max_height = 0.0;
			double max_width = 0.0;

			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
					cell_iter != cell_population.End(); ++cell_iter)
			{
				double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
				double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

				if (y > max_height)
				{
					max_height = y;
				}

				if (x > max_width)
				{
					max_width = x;
				}
			}

			// Define the crypt orifice
			double left_boundary = 4.0;
			double right_boundary = 6.0;
			double bottom_boundary = 4.0;

			// Create the crypt-orifice by removing the middle bit of the cell
			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
					cell_iter != cell_population.End(); ++cell_iter)
			{

				double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
				double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

				if (y == 0.0)
				{

					cell_iter->AddCellProperty(p_label);

				}
				else if ( (x >= left_boundary)&&(x <= right_boundary)&&(y > bottom_boundary) )
				{
					cell_iter->Kill();
				}
			}

			cell_population.RemoveDeadCells();

			// Un-Differentiate the epithelium
			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
					cell_iter != cell_population.End(); ++cell_iter)
			{

				double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
				double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

				if (y > bottom_boundary - 1.0)
				{
					if (y < bottom_boundary)
					{
						if ( (x >= left_boundary - 0.5)&&(x <= right_boundary + 0.5) )
						{
							cell_iter->SetCellProliferativeType(p_stem_type);
						}
					}
					else if ( (x >= left_boundary - 1.0)&&(x <= right_boundary + 1.0) )
					{
						cell_iter->SetCellProliferativeType(p_stem_type);
					}
					else if (y == max_height)
					{
						cell_iter->SetCellProliferativeType(p_stem_type);
					}
				}
			}

			//Output data to vtk format so we can visualise it in Paraview
			//			cell_population.AddCellWriter<GeneralisedCellAppliedForceWriter>();

			OffLatticeSimulation<2> simulator(cell_population);

			//Set output directory
			std::stringstream out;
			out << "/CRYPT/B_" << bm_stiffness << "_TC_" << target_curvature << "/V_" << vertical_force_magnitude
					<< "/L0_" << epithelial_epithelial_resting_spring_length;
			std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
			simulator.SetOutputDirectory(output_directory);

			simulator.SetDt(M_DT);
			simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
			simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

			//Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
			//			MAKE_PTR(LinearSpringForceWithVariableRestLength<2>, p_spring_force);
			//			//			p_spring_force->SetCutOffLength(epithelial_epithelial_resting_spring_length);
			//			p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness); //Default is 15
			//			p_spring_force->SetEpithelialStromalSpringStiffness(epithelial_stromal_stiffness); //Default is 15
			//			p_spring_force->SetStromalStromalSpringStiffness(stromal_stromal_stiffness); //Default is 15
			//			p_spring_force->SetEpithelialEpithelialRestingSpringLength(epithelial_epithelial_resting_spring_length);
			//			p_spring_force->SetCutOffLength(radius_of_interaction);
			//			p_spring_force->SetMeinekeDivisionRestingSpringLength(division_separation);
			//			simulator.AddForce(p_spring_force);

			//			//Add basement membrane force
			//			MAKE_PTR(OverlappingSpheresBasedBasementMembraneForce, p_bm_force);
			//			p_bm_force->SetBasementMembraneParameter(bm_stiffness); //Equivalent to beta in SJD's papers
			//			p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
			//			p_bm_force->SetCryptGeometry(true); // Apply to crypt
			//			p_bm_force->ApplyPeriodicForce(true); // Apply if we have periodic mesh
			//			p_bm_force->SetLeftCryptBoundary(left_boundary);
			//			p_bm_force->SetRightCryptBoundary(right_boundary);
			//			p_bm_force->SetCutOffRadius(radius_of_interaction); //Set cut off radius for defining neighbours
			//			simulator.AddForce(p_bm_force);
			//
			//			//			Add vertical compression force
			//			MAKE_PTR(VerticalCompressionForce, p_vertical_force);
			//			p_vertical_force->SetForceMagnitude(vertical_force_magnitude);
			//			simulator.AddForce(p_vertical_force);

			//Add anoikis-based cell killer
			MAKE_PTR_ARGS(AnoikisCellKiller, p_anoikis_killer, (&cell_population));
			p_anoikis_killer->SetCutOffRadius(radius_of_interaction);
			simulator.AddCellKiller(p_anoikis_killer);

			// Add modifier to track positions
			MAKE_PTR(PositionAndForceTrackingModifier<2>, p_position_tracking_modifier);
			simulator.AddSimulationModifier(p_position_tracking_modifier);

			//			Fix the bottom row of cells
			c_vector<double, 2> point, normal;

			point(0) = 0.0;
			point(1) = 0.25;
			normal(0) = 0.0;
			normal(1) = -1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc1);

			simulator.Solve();
			//
			//			// Get a pointer to the vertical compression force to unload the force
			//			boost::shared_ptr<VerticalCompressionForce> p_vertical_compression_force =
			//					boost::static_pointer_cast<VerticalCompressionForce>(simulator.rGetForceCollection()[2]);
			//
			//			p_vertical_compression_force->SetForceMagnitude(0.0);
			//
			//			simulator.SetEndTime(M_END_TIME + M_SECOND_END_TIME); // Set new end time
			//
			//			simulator.Solve(); // Re-run simulation

			//Tidying up
			SimulationTime::Instance()->Destroy();
			SimulationTime::Instance()->SetStartTime(0.0);

		}

	}
};

#endif /* TESTSTRESSESINDEFORMATIONFOROSMODELS_HPP_ */
