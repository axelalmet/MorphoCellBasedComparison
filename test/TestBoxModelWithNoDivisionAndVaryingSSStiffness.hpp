#ifndef TESTBOXMODELWITHNOPROLIFERATIONANDVARYINGSSSTIFFNESS_HPP_
#define TESTBOXMODELWITHNOPROLIFERATIONANDVARYINGSSSTIFFNESS_HPP_

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
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "NonPeriodicBasementMembraneForce.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "LinearSpringForceWithVariableRestLength.hpp"
#include "NonPeriodicBasementMembraneForce.hpp" //BM force that takes into account pe
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestBoxModelNoDivisionAndVaryingSSStiffness : public AbstractCellBasedTestSuite
{
public:
	void TestVaryingSSStiffness() throw(Exception)
	{
		double dt = 0.005; //Set dt
		double end_time = 20.0; //Set end time
		double sampling_timestep = end_time/dt;

		//Set all the spring stiffness variables
        double epithelial_epithelial_stiffness = 15.0;
        double epithelial_stromal_stiffness = 15.0;

		//Set the number of cells across and down for the array
		unsigned cells_across = 21;
		unsigned cells_up = 5;
		unsigned ghosts = 2; //Set the number of ghost node layers

		for (double SS = 2.0; SS < 4.0; SS++)
		{
            double stromal_stromal_stiffness = 15.0*SS;

			for (double BM = 0.0; BM < 16.0; BM ++)
			{
				//Set the basement membrane force parameters
				double bm_stiffness = BM;
				double target_curvature = 0.0;

				for (double EE = 0.0; EE < 41.0; EE ++)
				{
					double epithelial_epithelial_resting_spring_length = 1.0 + 0.1*EE;

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
							double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,13) and set in the past
							p_cycle_model->SetBirthTime(-birth_time);

							CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
							p_cell->InitialiseCellCycleModel(); // For paranoia really.

							p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated

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
					std::stringstream out;
					out << "SS_" << stromal_stromal_stiffness << "/B_" << bm_stiffness 
									<< "/L0_" << epithelial_epithelial_resting_spring_length;
					std::string output_directory = "BoxModelVaryingSSStiffness/" + out.str();
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

					//Tidying up
					SimulationTime::Instance()->Destroy();
					SimulationTime::Instance()->SetStartTime(0.0);

					// Calculate the relevant quantities. First we need to get the epithelium.
					std::vector<unsigned> epithelial_indices;
					for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
							cell_iter != simulator.rGetCellPopulation().End(); ++cell_iter)
					{
						if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
						{
							unsigned index = simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);

							epithelial_indices.push_back(index);
						}
					}

					unsigned num_epithelial_cells = epithelial_indices.size(); // Get the number of epithelial cells

					//Now we need to check whether the epithelium is still a continuous line.
					bool is_epithelium_confluent = true; // Initialise boolean

					// To check confluence, we iterate over all the springs in the tessellation and count the number of epithelial-epithelial connections.
					unsigned num_epithelial_epithelial_connections = 0;
					
					MeshBasedCellPopulationWithGhostNodes<2>* p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&simulator.rGetCellPopulation());

					p_cell_population->CreateVoronoiTessellation();

					for (MeshBasedCellPopulationWithGhostNodes<2>::SpringIterator spring_iter = p_cell_population->SpringsBegin();
					spring_iter != p_cell_population->SpringsEnd(); ++spring_iter)
					{
						unsigned nodeA_global_index = spring_iter.GetNodeA()->GetIndex();
						unsigned nodeB_global_index = spring_iter.GetNodeB()->GetIndex();

						//If they're both real nodes
						if ( (!p_cell_population->IsGhostNode(nodeA_global_index))&&(!p_cell_population->IsGhostNode(nodeB_global_index)) )
						{
							//Get the cells corresponding to them
							CellPtr p_cell_A = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index);
							CellPtr p_cell_B = p_cell_population->GetCellUsingLocationIndex(nodeB_global_index);
							//If they're both epithelial cells

							boost::shared_ptr<AbstractCellProperty> p_type_A = p_cell_A->GetCellProliferativeType();
							boost::shared_ptr<AbstractCellProperty> p_type_B = p_cell_B->GetCellProliferativeType();
							if ( (p_type_A == p_type_B)&&(!p_type_A->IsType<DifferentiatedCellProliferativeType>()) )
							{
								num_epithelial_epithelial_connections += 1;
							}
						}
					}

					if (num_epithelial_epithelial_connections != (num_epithelial_cells - 1) )	
					{
						is_epithelium_confluent = false;
					}			

					if (!is_epithelium_confluent)
					{
						break;
					}
					else
					{
						//Create results file for model
						OutputFileHandler results_handler(output_directory + "/", false); //Create output file handler

						//Create file
						std::string epithelium_data_filename = "epithelium_data.dat";
						out_stream epithelium_data_file = results_handler.OpenOutputFile(epithelium_data_filename);

						// The initial and growth length are based on the resting spring lengths, which are easy
						double initial_length = 1.0*num_epithelial_epithelial_connections;
						double growth_length = epithelial_epithelial_resting_spring_length*num_epithelial_epithelial_connections;

						// Calculate the current length
						double current_length = 0.0;

						for (MeshBasedCellPopulationWithGhostNodes<2>::SpringIterator spring_iter = p_cell_population->SpringsBegin();
						spring_iter != p_cell_population->SpringsEnd(); ++spring_iter)
						{
							unsigned nodeA_global_index = spring_iter.GetNodeA()->GetIndex();
							unsigned nodeB_global_index = spring_iter.GetNodeB()->GetIndex();

							//If they're both real nodes
							if ( (!p_cell_population->IsGhostNode(nodeA_global_index))&&(!p_cell_population->IsGhostNode(nodeB_global_index)) )
							{
								//Get the cells corresponding to them
								CellPtr p_cell_A = p_cell_population->GetCellUsingLocationIndex(nodeA_global_index);
								CellPtr p_cell_B = p_cell_population->GetCellUsingLocationIndex(nodeB_global_index);
								//If they're both epithelial cells

								boost::shared_ptr<AbstractCellProperty> p_type_A = p_cell_A->GetCellProliferativeType();
								boost::shared_ptr<AbstractCellProperty> p_type_B = p_cell_B->GetCellProliferativeType();
								if ( (p_type_A == p_type_B)&&(!p_type_A->IsType<DifferentiatedCellProliferativeType>()) )
								{
									double distance_between_nodes = p_cell_population->rGetMesh().GetDistanceBetweenNodes(nodeA_global_index, nodeB_global_index);
									current_length += distance_between_nodes;
								}
							}
						}

						//We need to sort the epithelial cells by x-coordinates to determine 'wrinkliness'.
						std::vector<std::pair<double, unsigned> > x_coords_and_epithelial_indices;

						for (unsigned i = 0; i < epithelial_indices.size(); i++)
						{
							unsigned epithelial_index = epithelial_indices[i]; // Get the node index

							CellPtr epithelial_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(epithelial_index); //Get the cell

							//Get the cell location
							double x = simulator.rGetCellPopulation().GetLocationOfCellCentre(epithelial_cell)[0];

							std::pair<double, unsigned> position_and_index = std::make_pair(x, epithelial_index);

							x_coords_and_epithelial_indices.push_back(position_and_index);
						}

						std::sort(x_coords_and_epithelial_indices.begin(), x_coords_and_epithelial_indices.end());

						// Now we can calculate the wrinkliness
						double wrinkliness = 0.0;

						for (unsigned i = 0; i < (x_coords_and_epithelial_indices.size() - 1); i++)
						{
							unsigned first_node_index = x_coords_and_epithelial_indices[i].second;
							unsigned second_node_index = x_coords_and_epithelial_indices[i + 1].second;

							CellPtr first_epithelial_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(first_node_index);
							CellPtr second_epithelial_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(second_node_index);

							c_vector<double, 2> first_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(first_epithelial_cell);
							c_vector<double, 2> second_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(second_epithelial_cell);

							double wrinkliness_increment = fabs((second_location[1] - first_location[1])/(second_location[0] - first_location[0]));

							wrinkliness += wrinkliness_increment;
						}

						wrinkliness /= (num_epithelial_epithelial_connections);

						// Finally, we can write the data to the file
						*epithelium_data_file << initial_length << "\t" << growth_length << "\t" << current_length << "\t" << wrinkliness << "\n";
						
						// For completeness, we should write the positions of the epithelium to file as well.
						std::string epithelium_positions_filename = "epithelium_positions.dat";
						out_stream epithelium_positions_file = results_handler.OpenOutputFile(epithelium_positions_filename);

						for (unsigned i = 0; i < x_coords_and_epithelial_indices.size(); i++)
						{
							unsigned node_index = x_coords_and_epithelial_indices[i].second;

							CellPtr epithelial_cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(node_index);

							c_vector<double, 2> location = simulator.rGetCellPopulation().GetLocationOfCellCentre(epithelial_cell);

							*epithelium_positions_file << location[0] << "\t" << location[1] << "\n";
						}

					}
				}
			}
		}
	}
};

#endif /* TESTBOXMODELWITHNODIVISIONANDVARYINGESSTIFFNESS_HPP_ */
