#ifndef TESTBOXMODELWITHNODIVISION_HPP_
#define TESTBOXMODELWITHNODIVISION_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"

#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CylindricalHoneycombMeshGenerator.hpp" //Generates mesh
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

class TestBoxModelWithNoDivision : public AbstractCellBasedTestSuite
{
public:
	void TestDivisionEvent() throw(Exception)
	{
		double dt = 0.005; //Set dt
		double end_time = 9.0; //Set end time
		double sampling_timestep = 0.1/dt;

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 15.0;
		double epithelial_stromal_stiffness = 15.0;
		double stromal_stromal_stiffness = 15.0;

		//Set the number of cells across and down for the array
		unsigned cells_across = 10;
		unsigned cells_up = 5;
		unsigned ghosts = 2; //Set the number of ghost node layers

		//Set the basement membrane force parameters
		double bm_stiffness = 10.0;
		double target_curvature = 0.0;
		double epithelial_epithelial_resting_spring_length = 1.0;

		//Generate the mesh
		HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		//Get the real indices
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		//Translate mesh 1.5 units left and 1.5 units down, so that we have a layer of fixed cells around the boundary for the BC.
		c_vector<double, 2> translate_left = zero_vector<double>(2);
		translate_left(0) = -1.0;

		c_vector<double, 2> translate_down = zero_vector<double>(2);
		translate_down(1) = 0.0;

		//Translate mesh appropriately
		p_mesh->Translate(translate_left);
		p_mesh->Translate(translate_down);

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


		//Create the vector of cells
		//Create shared pointers for cell and mutation states
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

		//Create tissue of cells. Initially we set them all to be differentiated
		std::vector<CellPtr> cells; //Create vector of cells

		for (unsigned i = 0; i<location_indices.size(); i++)
		{

			unsigned node_index = location_indices[i];
			double x = p_mesh->GetNode(node_index)->rGetLocation()[0];
			double y = p_mesh->GetNode(node_index)->rGetLocation()[1];

			if (y != max_height)
			{

				//Set stochastic duration based cell cycle
				NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Don't give them any cell cycle model yet.
				double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,13) and set in the past
				p_cycle_model->SetBirthTime(-birth_time);

				CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
				p_cell->InitialiseCellCycleModel(); // For paranoia really.

				p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated

				cells.push_back(p_cell);

			}
			else
			{

				if (x == 3.0)
				{
					UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel(); //Don't give them any cell cycle model yet.
					p_cycle_model->SetDimension(2);
					p_cycle_model->SetBirthTime(-1.0);
					p_cycle_model->SetMaxCellCycleDuration(1.0);

					CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
					p_cell->InitialiseCellCycleModel(); // For paranoia really.

					p_cell->SetCellProliferativeType(p_stem_type);

					cells.push_back(p_cell);
				}
				else
				{
					//Set stochastic duration based cell cycle
					// UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel(); //Uniformly distributed cell cycle times
					NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Don't give them any cell cycle model yet.
					double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,13) and set in the past
					p_cycle_model->SetBirthTime(-birth_time);

					CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
					p_cell->InitialiseCellCycleModel(); // For paranoia really.

					p_cell->SetCellProliferativeType(p_stem_type);

					cells.push_back(p_cell);
				}

			}

		}

		//Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

		//Output data to vtk format so we can visualise it in Paraview
		cell_population.SetWriteVtkAsPoints(true);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		//Any cell with the maximum height is converted to be a proliferative cell

		//Create ring of stem cells
		for (unsigned i = 0; i < location_indices.size(); i++)
		{
			unsigned cell_index = location_indices[i];
			CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
			double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];

			if(y==max_height)
			{
				cell_iter->SetCellProliferativeType(p_stem_type);
			}
		}

		OffLatticeSimulation<2> simulator(cell_population);

		//Set output directory
		simulator.SetOutputDirectory("BoxModelWithDivisionEvent");

		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(sampling_timestep); //Sample the simulation at every hour
		simulator.SetEndTime(end_time); //Hopefully this is long enough for a steady state

		//Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
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
		simulator.AddForce(p_bm_force);

		//Fix the bottom row of cells
		c_vector<double, 2> point, normal;

		point(0) = 0.0;
		point(1) = 0.75;
		normal(0) = 0.0;
		normal(1) = -1.0;
		MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc1);

		point(0) = -0.25;
		point(1) = 0.0;
		normal(0) = -1.0;
		normal(1) = 0.0;
		MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc2);

		point(0) = max_width - 0.75;
		point(1) = 0.0;
		normal(0) = 1.0;
		normal(1) = 0.0;
		MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc3);

		simulator.Solve();

	}
};

#endif /* TESTBOXMODELWITHNODIVISION_HPP_ */
