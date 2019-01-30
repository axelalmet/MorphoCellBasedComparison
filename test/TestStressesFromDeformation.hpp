#ifndef TESTBOXMODELWITHNODIVISION_HPP_
#define TESTBOXMODELWITHNODIVISION_HPP_

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
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "LinearSpringForceWithVariableRestLength.hpp"
#include "NonPeriodicBasementMembraneForce.hpp" //BM force (based off Dunn et al. (2012))
#include "OverlappingSpheresBasedBasementMembraneForce.hpp" //Overlapping spheres equivalent
#include "AnoikisCellKiller.hpp" // Anoikis-based cell killer
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "MeasuringStresses";

class TestBoxModelWithNoDivision : public AbstractCellBasedTestSuite
{
public:
	void TestFlatToBuckledEpithelium()
	{
		double dt = 0.01; //Set dt
		double end_time = 1.0; //Set end time
		double sampling_timestep = end_time/dt;

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 15.0;
		double epithelial_stromal_stiffness = 15.0;
		double stromal_stromal_stiffness = 15.0;

		//Set the number of cells across and down for the array
		unsigned cells_across = 20;
		unsigned cells_up = 10;
		unsigned ghosts = 2; //Set the number of ghost node layers

		//Set the basement membrane force parameters
		double bm_stiffness = 12.0;
		double target_curvature = 0.4;

		double left_boundary = 7.5;
		double right_boundary = 12.5;

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
			double y = p_mesh->GetNode(node_index)->rGetLocation()[1];

			if (y != max_height)
			{

				//Set stochastic duration based cell cycle
				NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Don't give them any cell cycle model yet.
				p_cycle_model->SetDimension(2);

				CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
				p_cell->InitialiseCellCycleModel(); // For paranoia really.

				p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated

				cells.push_back(p_cell);

			}
			else
			{

				//Set stochastic duration based cell cycle
				UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel(); //Uniformly distributed cell cycle times
				double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,13) and set in the past
				p_cycle_model->SetBirthTime(-birth_time);
				p_cycle_model->SetMinCellCycleDuration(10.0);
				p_cycle_model->SetMaxCellCycleDuration(14.0);

				CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
				p_cell->InitialiseCellCycleModel(); // For paranoia really.

				p_cell->SetCellProliferativeType(p_stem_type);

				cells.push_back(p_cell);

			}

		}

		//Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

		//Output data to vtk format so we can visualise it in Paraview
		cell_population.SetWriteVtkAsPoints(true);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		OffLatticeSimulation<2> simulator(cell_population);

		//Set output directory
		std::stringstream out;
		out << "/FLAT/VT/B_" << bm_stiffness << "_TC_" << target_curvature << "/L0_" << epithelial_epithelial_resting_spring_length;
		std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
		simulator.SetOutputDirectory(output_directory);

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
		p_bm_force->ApplyForceToCrypt(true);
		p_bm_force->SetLeftCryptBoundary(left_boundary);
		p_bm_force->SetRightCryptBoundary(right_boundary);
		simulator.AddForce(p_bm_force);

		//Add anoikis-based cell killer
		MAKE_PTR_ARGS(AnoikisCellKiller, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);

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

		//		// Calculate the relevant quantities. First we need to get the epithelium.
		//		std::vector<unsigned> epithelial_indices;
		//
		//		for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
		//				cell_iter != simulator.rGetCellPopulation().End(); ++cell_iter)
		//		{
		//			if (cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
		//			{
		//				unsigned index = simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);
		//
		//				epithelial_indices.push_back(index);
		//			}
		//		}
		//
		//
		////		MeshBasedCellPopulationWithGhostNodes<2>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&simulator.rGetCellPopulation());
		////
		////		p_cell_population->CreateVoronoiTessellation();
		//
		//		//Create results file for model
		//		OutputFileHandler results_handler(output_directory + "/", false); //Create output file handler
		//
		//		// For completeness, we should write the positions of the epithelium to file as well.
		//		std::string epithelium_positions_filename = "epithelium_positions.dat";
		//		out_stream epithelium_positions_file = results_handler.OpenOutputFile(epithelium_positions_filename);
		//
		//		// We also include the applied forces on the epithelium
		//		std::string epithelium_forces_filename = "epithelium_forces.dat";
		//		out_stream epithelium_forces_file = results_handler.OpenOutputFile(epithelium_forces_filename);
		//
		//		//We need to sort the epithelial cells by x-coordinates to determine 'wrinkliness'.
		//		std::vector<std::pair<double, unsigned> > x_coords_and_epithelial_indices;
		//
		//		for (unsigned i = 0; i < epithelial_indices.size(); i++)
		//		{
		//			unsigned epithelial_index = epithelial_indices[i]; // Get the node index
		//
		//			CellPtr epithelial_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(epithelial_index); //Get the cell
		//
		//			//Get the cell location
		//			double x = simulator.rGetCellPopulation().GetLocationOfCellCentre(epithelial_cell)[0];
		//
		//			std::pair<double, unsigned> position_and_index = std::make_pair(x, epithelial_index);
		//
		//			x_coords_and_epithelial_indices.push_back(position_and_index);
		//		}
		//
		//		std::sort(x_coords_and_epithelial_indices.begin(), x_coords_and_epithelial_indices.end());
		//
		//		MeshBasedCellPopulationWithGhostNodes<2>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&simulator.rGetCellPopulation());
		//
		//		p_cell_population->CreateVoronoiTessellation();
		//
		//		// Now iterate through the sorted vector
		//		for (unsigned i = 0; i < x_coords_and_epithelial_indices.size(); i++)
		//		{
		//			unsigned node_index = x_coords_and_epithelial_indices[i].second;
		//
		//			CellPtr epithelial_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(node_index);
		//
		//			// Get the cell position
		//			c_vector<double, 2> location = simulator.rGetCellPopulation().GetLocationOfCellCentre(epithelial_cell);
		//
		//			*epithelium_positions_file << location[0] << "\t" << location[1] << "\n";
		//
		//			// Get the applied force
		//			Node<2>* epithelial_node = p_cell_population->rGetMesh().GetNode(node_index);
		//
		//			c_vector<double, 2> applied_force = epithelial_node->rGetAppliedForce();
		//
		////			*epithelium_forces_file << applied_force[0] << "\t" << applied_force[1] << "\n";
		//		}

	}
	void TestBuckledCryptEpithelium()
	{
		double dt = 0.01; //Set dt
		double end_time = 1.0; //Set end time
		double sampling_timestep = end_time/dt;

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 45.0;
		double epithelial_stromal_stiffness = 15.0;
		double stromal_stromal_stiffness = 15.0;

		//Set the number of cells across and down for the array
		unsigned cells_across = 10;
		unsigned cells_up = 20;
		unsigned ghosts = 2; //Set the number of ghost node layers

		//Set the basement membrane force parameters
		double bm_stiffness = 12.0;
		double target_curvature = 0.4;

		double left_boundary = 2.5;
		double right_boundary = 5.0;
		double bottom_boundary = 5.0;

		double epithelial_epithelial_resting_spring_length = 1.0;

		//Generate the mesh
		HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		//Get the initial non-ghost indices
		std::vector<unsigned> initial_location_indices = generator.GetCellLocationIndices();

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

		// Re-defined the actual non-ghost indices that defines a crypt shape
		std::vector<unsigned> location_indices;

		for (unsigned i = 0; i < initial_location_indices.size(); i++)
		{
			unsigned node_index = initial_location_indices[i];
			double x = p_mesh->GetNode(node_index)->rGetLocation()[0];
			double y = p_mesh->GetNode(node_index)->rGetLocation()[1];

			if ( (x < left_boundary)||(x > right_boundary) )
			{
				location_indices.push_back(node_index);
			}
			else
			{
				if (y < bottom_boundary)
				{
					location_indices.push_back(node_index);
				}
			}
		}

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

			//Set stochastic duration based cell cycle
			NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Uniformly distributed cell cycle times

			CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
			p_cell->InitialiseCellCycleModel(); // For paranoia really.

			p_cell->SetCellProliferativeType(p_diff_type);

			cells.push_back(p_cell);

		}

		//Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

		//Create the epithelium of stem cells
		for (unsigned i = 0; i < location_indices.size(); i++)
		{
			unsigned cell_index = location_indices[i];
			CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
			double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
			double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];

			//"un-differentiate" the epithelium
			if (y == max_height)
			{
				Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);

				//Iterate over the elements (triangles) containing the nodes
				for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
						iter != p_node->ContainingElementsEnd();
						++iter)
				{
					bool element_contains_ghost_nodes = false;

					// Get a pointer to the element
					Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

					// Check whether its triangulation contains a ghost node
					for (unsigned local_index=0; local_index<3; local_index++)
					{
						unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

						if (cell_population.IsGhostNode(nodeGlobalIndex) == true)
						{
							element_contains_ghost_nodes = true;
							break; 				// This should break out of the inner for loop
						}
					}

					if(element_contains_ghost_nodes)
					{
						cell_iter->SetCellProliferativeType(p_stem_type);
					}
				}
			}
			else if ( (x >= left_boundary - 1.0)&&(x <= right_boundary + 1.0)&&(y > bottom_boundary - 1.0) )
			{
				Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);

				//Iterate over the elements (triangles) containing the nodes
				for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
						iter != p_node->ContainingElementsEnd();
						++iter)
				{
					bool element_contains_ghost_nodes = false;

					// Get a pointer to the element
					Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

					// Check whether its triangulation contains a ghost node
					for (unsigned local_index=0; local_index<3; local_index++)
					{
						unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

						if (cell_population.IsGhostNode(nodeGlobalIndex) == true)
						{
							element_contains_ghost_nodes = true;
							break; 				// This should break out of the inner for loop
						}
					}

					if(element_contains_ghost_nodes)
					{
						cell_iter->SetCellProliferativeType(p_stem_type);
					}
				}
			}
		}

		//Output data to vtk format so we can visualise it in Paraview
		cell_population.SetWriteVtkAsPoints(true);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		OffLatticeSimulation<2> simulator(cell_population);

		//Set output directory
		std::stringstream out;
		out << "/CRYPT/VT/B_" << bm_stiffness << "_TC_" << target_curvature << "/L0_" << epithelial_epithelial_resting_spring_length;
		std::string output_directory = M_OUTPUT_DIRECTORY + out.str();
		simulator.SetOutputDirectory(output_directory);

		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(sampling_timestep); //Sample the simulation at every hour
		simulator.SetEndTime(end_time); //Hopefully this is long enough for a steady state

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
		p_bm_force->ApplyForceToCrypt(true);
//		p_bm_force->SetLeftCryptBoundary(left_boundary);
//		p_bm_force->SetRightCryptBoundary(right_boundary);
		simulator.AddForce(p_bm_force);

		//Add anoikis-based cell killer
		MAKE_PTR_ARGS(AnoikisCellKiller, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);

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
