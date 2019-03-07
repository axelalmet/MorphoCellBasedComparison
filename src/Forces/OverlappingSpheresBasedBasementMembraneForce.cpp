#include "OverlappingSpheresBasedBasementMembraneForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellProperty.hpp"
#include "Debug.hpp"
#include "Exception.hpp"

/*
 * MODIFIED FOR PERSONAL USE BY AXEL ALMET
 * Last modified: Feb 22 2016
 */


/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
OverlappingSpheresBasedBasementMembraneForce::OverlappingSpheresBasedBasementMembraneForce()
:  AbstractForce<2>(),
   mBasementMembraneParameter(DOUBLE_UNSET),
   mTargetCurvature(DOUBLE_UNSET),
   mLeftBoundary(DOUBLE_UNSET),
   mRightBoundary(DOUBLE_UNSET),
   mUsePositionDependentMembraneForce(false),
   mApplyForceToCrypt(true),
   mMembraneForceMultiplier(DOUBLE_UNSET),
   mCutOffRadius(1.5),
   mApplyPeriodicForce(false)
   {
	// Sets up output file
	//	OutputFileHandler output_file_handler("CurvatureData/", false);
	//	mMeinekeOutputFile = output_file_handler.OpenOutputFile("results.curvature");
   }

OverlappingSpheresBasedBasementMembraneForce::~OverlappingSpheresBasedBasementMembraneForce()
{
	//    mMeinekeOutputFile->close();
}

void OverlappingSpheresBasedBasementMembraneForce::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double OverlappingSpheresBasedBasementMembraneForce::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}


void OverlappingSpheresBasedBasementMembraneForce::SetTargetCurvature(double targetCurvature)
{
	mTargetCurvature = targetCurvature;
}


double OverlappingSpheresBasedBasementMembraneForce::GetTargetCurvature()
{
	return mTargetCurvature;
}

void OverlappingSpheresBasedBasementMembraneForce::SetLeftCryptBoundary(double leftBoundary)
{
	mLeftBoundary = leftBoundary;
}

double OverlappingSpheresBasedBasementMembraneForce::GetLeftCryptBoundary()
{
	return mLeftBoundary;
}

void OverlappingSpheresBasedBasementMembraneForce::SetRightCryptBoundary(double rightBoundary)
{
	mRightBoundary = rightBoundary;
}

double OverlappingSpheresBasedBasementMembraneForce::GetRightCryptBoundary()
{
	return mRightBoundary;
}

void OverlappingSpheresBasedBasementMembraneForce::SetCryptGeometry(bool applyForceToCrypt)
{
	mApplyForceToCrypt = applyForceToCrypt;
}

bool OverlappingSpheresBasedBasementMembraneForce::GetCryptGeometryCheck()
{
	return mApplyForceToCrypt;
}

void OverlappingSpheresBasedBasementMembraneForce::SetPositionDependentMultiplier(bool usePositionDependentMembraneForce, double MembraneForceMultiplier)
{
	mUsePositionDependentMembraneForce = usePositionDependentMembraneForce;
	mMembraneForceMultiplier = MembraneForceMultiplier;
}


double OverlappingSpheresBasedBasementMembraneForce::GetPositionDependentMultiplier()
{
	return mMembraneForceMultiplier;
}

double OverlappingSpheresBasedBasementMembraneForce::GetCutOffRadius()
{
	return mCutOffRadius;
}

void OverlappingSpheresBasedBasementMembraneForce::SetCutOffRadius(double cutOffRadius)
{
	mCutOffRadius = cutOffRadius;
}

void OverlappingSpheresBasedBasementMembraneForce::ApplyPeriodicForce(bool applyPeriodicForce)
{
	mApplyPeriodicForce = applyPeriodicForce;
}

bool OverlappingSpheresBasedBasementMembraneForce::IsPeriodicForceApplied()
{
	return mApplyPeriodicForce;
}


void OverlappingSpheresBasedBasementMembraneForce::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
	std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
	rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

c_vector<double,2> OverlappingSpheresBasedBasementMembraneForce::GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	// Create a vector to store the y-coordinates of the lowest point of the crypt base and the highest point of the
	// crypt orifice
	c_vector<double,2> height_extremes;

	double max_height = 0.0;
	double min_height = DBL_MAX;

	double current_height_coordinate;

	// We iterate over all cells in the tissue, and deal only with those that are epithelial cells
	for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
			cell_iter != rCellPopulation.End();
			++cell_iter)
	{

		// Need these to not be labelled cells
		if ( (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) )
		{
			Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);

			current_height_coordinate = p_node->rGetLocation()[1];

			if (current_height_coordinate > max_height)
			{
				max_height = current_height_coordinate;
			}
			else if (current_height_coordinate < min_height)
			{
				min_height = current_height_coordinate;
			}
		}
	}

	height_extremes[0] = max_height;
	height_extremes[1] = min_height;

	return height_extremes;
}

c_vector<double,2> OverlappingSpheresBasedBasementMembraneForce::GetCryptWidthExtremes(AbstractCellPopulation<2>& rCellPopulation)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	// Create a vector to store the y-coordinates of the lowest point of the crypt base and the highest point of the
	// crypt orifice
	c_vector<double,2> width_extremes;

	double max_width = 0.0;
	double min_width = DBL_MAX;

	double current_width_coordinate;

	// We iterate over all cells in the tissue, and deal only with those that are epithelial cells
	for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
			cell_iter != rCellPopulation.End();
			++cell_iter)
	{

		// Need these to not be labelled cells
		if ( (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) )
		{
			Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);

			current_width_coordinate = p_node->rGetLocation()[0];

			if (current_width_coordinate > max_width)
			{
				max_width = current_width_coordinate;
			}
			else if (current_width_coordinate < min_width)
			{
				min_width = current_width_coordinate;
			}
		}
	}

	width_extremes[0] = max_width;
	width_extremes[1] = min_width;

	return width_extremes;
}

/*
 * Method to calculate the tangent line at a point, based on approximating the epithelium positions
 * by a cosine curve.
 */
c_vector<double, 2> OverlappingSpheresBasedBasementMembraneForce::GetCosineBasedTangentVector(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> point)
{
	// Get the minimal and maximal x-coordinates of epithelium
	c_vector<double, 2> width_extremes = GetCryptWidthExtremes(rCellPopulation);
	double max_width = width_extremes[0];
	double min_width = width_extremes[1];

	// Get the minimal and maximal y-coordinates of epithelium
	c_vector<double, 2> height_extremes = GetCryptHeightExtremes(rCellPopulation);
	double max_height = height_extremes[0];
	double min_height = height_extremes[1];

	// Define width of crypt
	double crypt_width = max_width - min_width;

	// Define crypt height
	double crypt_height = max_height - min_height;

	// Tangent defined by derivative of approximated cosine function
	c_vector<double, 2> tangent_vector;
	double x = point[0];

	// y-component is derivative of cos function approximation of epithelium, with appropriate scalings.
	tangent_vector[0] = 1.0;
	tangent_vector[1] = -(M_PI/crypt_width)*crypt_height*sin(2.0*M_PI/crypt_width*(x - min_width));

	tangent_vector /= norm_2(tangent_vector); // Convert to unit vector

	// Return a double-length vector---this works best for tracking neighbours
	tangent_vector *= 2.0;

	return tangent_vector;
}

/*
 * Return vector of epithelial indices that are close to the considered epithelial node,
 * but based on an approximated cosine approximation
 */
std::vector<unsigned> OverlappingSpheresBasedBasementMembraneForce::GetClosestNeighboursBasedOnCosineApproximation(AbstractCellPopulation<2>& rCellPopulation, unsigned epithelialIndex)
{
	// Initialise vector
	std::vector<unsigned> closest_neighbours;

	// Get the epithelial indices
	std::vector<unsigned> epithelial_indices = GetEpithelialIndices(rCellPopulation);

	// Get the location of the epithelial node
	c_vector<double, 2> epithelial_location = rCellPopulation.GetNode(epithelialIndex)->rGetLocation();

	// Sweep through the indices
	for (unsigned i = 0; i < epithelial_indices.size(); i++)
	{
		unsigned neighbour_index = epithelial_indices[i];

		if (neighbour_index != epithelialIndex)
		{

			c_vector<double, 2> neighbour_location = rCellPopulation.GetNode(neighbour_index)->rGetLocation();

			if (norm_2(neighbour_location - epithelial_location) < 2.0 ) // Distance of two is based on length of tangent vector, i.e. long enough to catch enough neighbours
			{
				closest_neighbours.push_back(neighbour_index);
			}

		}
	}

	return closest_neighbours;
}

/*
 * Return the nearest neighbour based on vector projections from the tangent vector
 * at a point along the cosine approximation of the epithelium.
 */
unsigned OverlappingSpheresBasedBasementMembraneForce::GetNearestNeighboursAlongCosineApproximation(AbstractCellPopulation<2>& rCellPopulation, unsigned epithelialIndex)
{
	double min_scalar_projection = DBL_MAX;
	double min_projection_index = 0;

	// Get the closest neighbours, based on the cosine approximation
	std::vector<unsigned> closest_neighbours = GetClosestNeighboursBasedOnCosineApproximation(rCellPopulation, epithelialIndex);

	// Get the location of the considered epithelial node
	c_vector<double, 2> epithelial_location = rCellPopulation.GetNode(epithelialIndex)->rGetLocation();

	// Calculate the tangent vector at the epithelial node
	c_vector<double, 2> tangent_vector = GetCosineBasedTangentVector(rCellPopulation, epithelial_location);

	for (unsigned i = 0; i < closest_neighbours.size(); i++)
	{
		// Get the location of the neighbour
		unsigned neighbour_index = closest_neighbours[i];
		c_vector<double, 2> neighbour_location = rCellPopulation.GetNode(neighbour_index)->rGetLocation();

		// Get the scalar projection of the relative position vector onto the tangent vector
		double scalar_projection = inner_prod(neighbour_location - epithelial_location, tangent_vector)/norm_2(tangent_vector);

		if (scalar_projection < min_scalar_projection)
		{
			min_scalar_projection = scalar_projection;
			min_projection_index = neighbour_index;
		}
	}

	return min_projection_index;
}


/*
 * Function to find the curvature along three points, using the method previously described by SJD
 * Method has been adjusted to account for periodic meshes etc, i.e. heavy use of GetVectorFromAtoB
 */
double OverlappingSpheresBasedBasementMembraneForce::FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
		c_vector<double, 2> leftPoint,
		c_vector<double, 2> centrePoint,
		c_vector<double, 2> rightPoint)
{

	//Get the relevant vectors (all possible differences)
	c_vector<double, 2> left_to_centre = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftPoint, centrePoint);
	c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centrePoint, rightPoint);
	c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftPoint, rightPoint);

	// Firstly find the parametric intervals
	double left_s = pow(pow(left_to_centre[0],2) + pow(left_to_centre[1],2), 0.5);
	double right_s = pow(pow(centre_to_right[0],2) + pow(centre_to_right[1],2), 0.5);
	//	(left_s, right_s);

	double sum_intervals = left_s + right_s;

	//Calculate finite difference of first derivatives
	double x_prime = (left_to_right[0])/sum_intervals;
	double y_prime = (left_to_right[1])/sum_intervals;

	//Calculate finite difference of second derivatives
	double x_double_prime = 2*(left_s*centre_to_right[0] - right_s*left_to_centre[0])/(left_s*right_s*sum_intervals);
	double y_double_prime = 2*(left_s*centre_to_right[1] - right_s*left_to_centre[1])/(left_s*right_s*sum_intervals);

	//Calculate curvature using formula
	double curvature = (x_prime*y_double_prime - y_prime*x_double_prime)/pow((pow(x_prime,2) + pow(y_prime,2)), 1.5);

	return curvature;
}

/*
 * Method to find all the epithelial cells that make up the monolayer
 */
std::vector<unsigned> OverlappingSpheresBasedBasementMembraneForce::GetEpithelialIndices(AbstractCellPopulation<2>& rCellPopulation)
{
	//Create the vector of epithelial cell indices
	std::vector<unsigned> epithelial_indices;

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

		for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
				cell_iter != p_tissue->End();
				++cell_iter)
		{
			//Get the cell type
			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			if(!p_type->IsType<DifferentiatedCellProliferativeType>()) //If we have an epithelial cell
			{
				unsigned node_index = p_tissue->GetLocationIndexUsingCell(*cell_iter);
				epithelial_indices.push_back(node_index);
			}
		}
	}

	return epithelial_indices;
}

/*
 * Method to get the neighbouring nodes that are epithelial cells
 */
std::vector<unsigned> OverlappingSpheresBasedBasementMembraneForce::GetNeighbouringEpithelialIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::vector<unsigned> neighbouring_epithelial_indices;

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

		//Get cut off radius for defining neighbourhood
		double radius = GetCutOffRadius();

		// Find the indices of the elements owned by this node
		std::set<unsigned> neighbouring_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(nodeIndex, radius);

		// Iterate over these elements
		for (std::set<unsigned>::iterator elem_iter = neighbouring_indices.begin();
				elem_iter != neighbouring_indices.end();
				++elem_iter)
		{
			//Get the cell according to the index
			CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

			//Get the cell type
			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			//if the cell is not differentiated and thus an epithelial cell, we add it to the vector
			if(!p_type->IsType<DifferentiatedCellProliferativeType>())
			{
				neighbouring_epithelial_indices.push_back(*elem_iter);
			}

		}
	}

	return neighbouring_epithelial_indices;
}

/*
 * Method to get the epithelial indices that have a left and right neighbour.
 * Returns a map with pairs <index, <left, right>> for convention.
 */
std::map<unsigned, std::pair<unsigned, unsigned> > OverlappingSpheresBasedBasementMembraneForce::GetEpithelialIndicesAndTheirLeftAndRightEpithelialNeighbours(AbstractCellPopulation<2>& rCellPopulation)
{
	//Declare the map
	std::map<unsigned, std::pair<unsigned, unsigned> > epithelial_indices_and_their_left_and_right_neighbours;

	//Get the epithelial indices
	std::vector<unsigned> epithelial_indices = GetEpithelialIndices(rCellPopulation);

	/* We are going to define the basement membrane using the x-coordinate of epithelial cells.
	 * Doesn't generalise to 3D, but it will give us a confluent layer.This method only applies if
	 * we are modelling a crypt. We have to change how we define neighbours if we are modelling an organoid
	 */

	bool is_force_applied_to_crypt = GetCryptGeometryCheck();
	if (is_force_applied_to_crypt)
	{
		std::vector<std::pair<double, unsigned> > epithelial_indices_and_x_coordinates; //Define vector of pairs, so that we may sort by the x-coordinate

		for (unsigned i = 0; i < epithelial_indices.size(); i++)
		{
			unsigned epithelial_index = epithelial_indices[i]; //Get node index
			double epithelial_x_coordinate = rCellPopulation.GetNode(epithelial_index)->rGetLocation()[0]; //Get x-coordinate of node location

			//Make pair
			std::pair<double, unsigned> x_coordinate_and_index = std::make_pair(epithelial_x_coordinate, epithelial_index);

			epithelial_indices_and_x_coordinates.push_back(x_coordinate_and_index);
		}

		//Sort indices by the x-coordinate
		std::sort(epithelial_indices_and_x_coordinates.begin(), epithelial_indices_and_x_coordinates.end());

		//We now define the nodes and their left and right neighbours in the layer
		for (unsigned i = 0; i < epithelial_indices_and_x_coordinates.size(); i++) // Do not apply the force to the left-most and right-most cells
		{

			unsigned centre_node_index = epithelial_indices_and_x_coordinates[i].second; //Get index of centre node

			//Initialise left and right neighbours
			unsigned left_neighbour_index, right_neighbour_index;

			if (i == 0)
			{
				left_neighbour_index = epithelial_indices_and_x_coordinates[epithelial_indices_and_x_coordinates.size() - 1].second;
				right_neighbour_index = epithelial_indices_and_x_coordinates[i + 1].second;
			}
			else if (i == (epithelial_indices_and_x_coordinates.size() - 1) )
			{
				left_neighbour_index = epithelial_indices_and_x_coordinates[i - 1].second;
				right_neighbour_index = 0;
			}
			else
			{
				left_neighbour_index = epithelial_indices_and_x_coordinates[i - 1].second;
				right_neighbour_index = epithelial_indices_and_x_coordinates[i + 1].second;
			}

			epithelial_indices_and_their_left_and_right_neighbours[centre_node_index] = std::make_pair(left_neighbour_index, right_neighbour_index);

		}
	}
	else //If we are applying the force to an organoid, we need to sort the epithelial cells using polar coordinates (note this assumes tha the initial geomotry is that of an organoid)
	{
		std::vector<std::pair<double, unsigned> > epithelial_indices_and_angles; //Define vector of pairs, so that we may sort by the angles

		//First get centre of mass of epithelial nodes, to define a reference point for calculating the polar angle
		c_vector<double, 2> centre_of_mass = zero_vector<double>(2);

		for (unsigned i = 0; i < epithelial_indices.size(); i++)
		{
			unsigned epithelial_index = epithelial_indices[i]; //Get node index
			c_vector<double, 2> epithelial_location = rCellPopulation.GetNode(epithelial_index)->rGetLocation(); //Get x-coordinate of node location

			centre_of_mass += epithelial_location;
		}

		//Average centre of mass
		centre_of_mass /= epithelial_indices.size();

		for (unsigned k = 0; k < epithelial_indices.size(); k++)
		{
			unsigned cell_index = epithelial_indices[k]; //Get the node index

			CellPtr cell = rCellPopulation.GetCellUsingLocationIndex(cell_index); //Get the cell

			//Get the cell location
			double x = rCellPopulation.GetLocationOfCellCentre(cell)[0];
			double y = rCellPopulation.GetLocationOfCellCentre(cell)[1];

			//Make the pair of the angle and cell
			std::pair<double, unsigned> angle_cell;

			//Get point relative to the circle centre
			double rel_x = x - centre_of_mass[0];
			double rel_y = y - centre_of_mass[1];

			double circle_angle = atan(rel_y/rel_x); //Get initial angle argument

			if (rel_x<0.0) //If the point is in the second quadrant or third quadrant
			{
				circle_angle += M_PI;
			}
			else if ((rel_x>=0.0)&&(rel_y<0.0)) //Fourth quadrant
			{
				circle_angle += 2*M_PI;
			}

			angle_cell = std::make_pair(circle_angle, cell_index);

			epithelial_indices_and_angles.push_back(angle_cell); //Add the angle and node index
		}

		//Sort indices by the angle
		std::sort(epithelial_indices_and_angles.begin(), epithelial_indices_and_angles.end());

		//We now define the nodes and their left and right neighbours in the layer
		for (unsigned i = 0; i < epithelial_indices_and_angles.size(); i++)
		{
			unsigned centre_node_index = epithelial_indices_and_angles[i].second; //Get index of centre node

			//Initialise left and right neighbours
			unsigned left_neighbour_index, right_neighbour_index;

			if (i == 0) //If it is the first index, the 'left' neighbour is the last index due to periodicity
			{
				left_neighbour_index = epithelial_indices_and_angles[epithelial_indices_and_angles.size() - 1].second;
				right_neighbour_index = epithelial_indices_and_angles[1].second;
			}
			else if (i == (epithelial_indices_and_angles.size() - 1) ) //If it is the last node, the 'right' neighbour is the first index
			{
				left_neighbour_index = epithelial_indices_and_angles[epithelial_indices_and_angles.size() - 2].second;
				right_neighbour_index = epithelial_indices_and_angles[0].second;
			}
			else //Otherwise the left neighbour is the index before and the right is the index after
			{
				left_neighbour_index = epithelial_indices_and_angles[i - 1].second;
				right_neighbour_index = epithelial_indices_and_angles[i + 1].second;
			}

			epithelial_indices_and_their_left_and_right_neighbours[centre_node_index] = std::make_pair(left_neighbour_index, right_neighbour_index);
		}
	}

	return epithelial_indices_and_their_left_and_right_neighbours;
}

// Method to check if node is left or right-most cell
bool OverlappingSpheresBasedBasementMembraneForce::IsBoundaryNode(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	bool is_boundary_node = false;

	//Get the epithelial indices
	std::vector<unsigned> epithelial_indices = GetEpithelialIndices(rCellPopulation);
	unsigned num_epithelial_nodes = epithelial_indices.size();

	std::vector<std::pair<double, unsigned> > epithelial_indices_and_x_coordinates; //Define vector of pairs, so that we may sort by the x-coordinate

	for (unsigned i = 0; i < epithelial_indices.size(); i++)
	{
		unsigned epithelial_index = epithelial_indices[i]; //Get node index
		double epithelial_x_coordinate = rCellPopulation.GetNode(epithelial_index)->rGetLocation()[0]; //Get x-coordinate of node location

		//Make pair
		std::pair<double, unsigned> x_coordinate_and_index = std::make_pair(epithelial_x_coordinate, epithelial_index);

		epithelial_indices_and_x_coordinates.push_back(x_coordinate_and_index);
	}

	//Sort indices by the x-coordinate
	std::sort(epithelial_indices_and_x_coordinates.begin(), epithelial_indices_and_x_coordinates.end());

	// If the index corresponds to the first or last index, it has to be one of the boundary nodes
	if ( (nodeIndex == epithelial_indices_and_x_coordinates[0].second)||(nodeIndex == epithelial_indices_and_x_coordinates[num_epithelial_nodes - 1].second) )
	{
		is_boundary_node = true;
	}

	return is_boundary_node;
}

//Method to calculate the force due to the basement membrane on an epithelial cell
c_vector<double, 2> OverlappingSpheresBasedBasementMembraneForce::CalculateForceDueToBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned > > epithelialIndicesAndNeighbours, unsigned nodeIndex)
{
	//Get the left and right crypt boundaries and the target curvature
	double left_boundary = GetLeftCryptBoundary();
	double right_boundary = GetRightCryptBoundary();
	double target_curvature = GetTargetCurvature();

	//Get the neighbours of the node index
	std::pair<unsigned, unsigned> epithelial_neighbours = epithelialIndicesAndNeighbours[nodeIndex];
	unsigned left_node_index = epithelial_neighbours.first;
	unsigned right_node_index = epithelial_neighbours.second;

	//Get the location of all three points for the curvature calculation
	c_vector<double, 2> centre_point = rCellPopulation.GetNode(nodeIndex)->rGetLocation();
	c_vector<double, 2> left_point = rCellPopulation.GetNode(left_node_index)->rGetLocation();
	c_vector<double, 2> right_point = rCellPopulation.GetNode(right_node_index)->rGetLocation();

	double basement_membrane_parameter = GetBasementMembraneParameter(); //Get the basement membrane stiffness
	//Initialise the force vector
	c_vector<double, 2> force_due_to_basement_membrane;

	bool is_periodic_force = IsPeriodicForceApplied();

	if (!is_periodic_force)
	{
		if ( !IsBoundaryNode(rCellPopulation, nodeIndex) )
		{
			double curvature = FindParametricCurvature(rCellPopulation, left_point, centre_point, right_point);

			//Get the unit vectors from the centre points to its left and right neighbours
			c_vector<double, 2> centre_to_left = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, left_point);
			centre_to_left /= norm_2(centre_to_left); //Normalise vector

			c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, right_point);
			centre_to_right /= norm_2(centre_to_right); //Normalise vector

			/* Define direction using a formula derived (essentially solve for the vector w such that u.w = v.w, where w is the force direction
			 * and u = unit vector to the left and v = unit vector pointing to the right)
			 */
			c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_to_left, centre_to_right);
			//	PRINT_2_VARIABLES(left_point[0], left_point[1]);
			//	PRINT_2_VARIABLES(centre_point[0], centre_point[1]);
			//	PRINT_2_VARIABLES(right_point[0], right_point[1]);

			assert(norm_2(left_to_right) != 0.0);

			//Define force direction
			c_vector<double, 2> force_direction;
			force_direction(0) = left_to_right[1];
			force_direction(1) = -left_to_right[0];

			force_direction /= norm_2(left_to_right);

			/* We now ensure the vector is pointing in the appropriate direction
			 * (it will always point "down" and "outwards" initially).
			 * */
			if (target_curvature > 0.0) //If the force has overshot the target curvature, we need to reverse the force direction
			{
				if ( (curvature > 0.0)&&(curvature - target_curvature > 0.0) ) //If points look like V and the 'v' is too pointy, we send it away from the CoM
				{
					force_direction *= -1.0;
				}

			}
			else if (target_curvature < 0.0) //Similar situation but with "/\"
			{
				if ( (curvature < 0.0)&&(curvature - target_curvature < 0.0) )
				{
					force_direction *= -1.0;
				}
			}
			else //Reverse the force direction if we get a "V"
			{
				if (curvature > 0.0)
				{
					force_direction *= -1.0;
				}
			}

			//Again, the geometry of the model alters how we apply target curvature
			bool is_force_applied_to_crypt = GetCryptGeometryCheck();

			//If we are considering a crypt geometry
			if(is_force_applied_to_crypt)
			{
				// If we are looking at a boundary node, we apply no force.
				if (left_point[0] > right_point[0])
				{
					force_due_to_basement_membrane = zero_vector<double>(2);
				}
				else
				{
					//If the cell falls in the region of non-zero target curvature
					if( (centre_point[0] > left_boundary)&&(centre_point[0] < right_boundary) )
					{
						force_due_to_basement_membrane = basement_membrane_parameter*(fabs(curvature - target_curvature) )*force_direction;
					}
					else
					{
						//We take the absolute value of the local curvature as the force direction vector already accounts for which way
						//the epithelial node 'should' go
						force_due_to_basement_membrane = basement_membrane_parameter*fabs(curvature)*force_direction;
					}
				}
			}
			else //Else we are modelling organoid
			{
				force_due_to_basement_membrane = basement_membrane_parameter*(fabs(curvature - target_curvature) )*force_direction;
			}
		}
		else
		{
			force_due_to_basement_membrane = zero_vector<double>(2);
		}
	}
	else
	{
		double curvature = FindParametricCurvature(rCellPopulation, left_point, centre_point, right_point);

		//Get the unit vectors from the centre points to its left and right neighbours
		c_vector<double, 2> centre_to_left = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, left_point);
		centre_to_left /= norm_2(centre_to_left); //Normalise vector

		c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, right_point);
		centre_to_right /= norm_2(centre_to_right); //Normalise vector

		/* Define direction using a formula derived (essentially solve for the vector w such that u.w = v.w, where w is the force direction
		 * and u = unit vector to the left and v = unit vector pointing to the right)
		 */
		c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_to_left, centre_to_right);
		//	PRINT_2_VARIABLES(left_point[0], left_point[1]);
		//	PRINT_2_VARIABLES(centre_point[0], centre_point[1]);
		//	PRINT_2_VARIABLES(right_point[0], right_point[1]);

		assert(norm_2(left_to_right) != 0.0);

		//Define force direction
		c_vector<double, 2> force_direction;
		force_direction(0) = left_to_right[1];
		force_direction(1) = -left_to_right[0];

		force_direction /= norm_2(left_to_right);

		/* We now ensure the vector is pointing in the appropriate direction
		 * (it will always point "down" and "outwards" initially).
		 * */
		if (target_curvature > 0.0) //If the force has overshot the target curvature, we need to reverse the force direction
		{
			if ( (curvature > 0.0)&&(curvature - target_curvature > 0.0) ) //If points look like V and the 'v' is too pointy, we send it away from the CoM
			{
				force_direction *= -1.0;
			}

		}
		else if (target_curvature < 0.0) //Similar situation but with "/\"
		{
			if ( (curvature < 0.0)&&(curvature - target_curvature < 0.0) )
			{
				force_direction *= -1.0;
			}
		}
		else //Reverse the force direction if we get a "V"
		{
			if (curvature > 0.0)
			{
				force_direction *= -1.0;
			}
		}

		//Again, the geometry of the model alters how we apply target curvature
		bool is_force_applied_to_crypt = GetCryptGeometryCheck();

		//If we are considering a crypt geometry
		if(is_force_applied_to_crypt)
		{
			// If we are looking at a boundary node, we apply no force.
			if (left_point[0] > right_point[0])
			{
				force_due_to_basement_membrane = zero_vector<double>(2);
			}
			else
			{
				//If the cell falls in the region of non-zero target curvature
				if( (centre_point[0] > left_boundary)&&(centre_point[0] < right_boundary) )
				{
					force_due_to_basement_membrane = basement_membrane_parameter*(fabs(curvature - target_curvature) )*force_direction;
				}
				else
				{
					//We take the absolute value of the local curvature as the force direction vector already accounts for which way
					//the epithelial node 'should' go
					force_due_to_basement_membrane = basement_membrane_parameter*fabs(curvature)*force_direction;
				}
			}
		}
		else //Else we are modelling organoid
		{
			force_due_to_basement_membrane = basement_membrane_parameter*(fabs(curvature - target_curvature) )*force_direction;
		}
	}

	return force_due_to_basement_membrane;
}

//Method overriding the virtual method for AbstractForce. The crux of what really needs to be done.
void OverlappingSpheresBasedBasementMembraneForce::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{

	//Get the epithelial nodes that have a left and right nearest neighbour
	std::map<unsigned, std::pair<unsigned, unsigned> > epithelial_nodes_and_their_neighbours = GetEpithelialIndicesAndTheirLeftAndRightEpithelialNeighbours(rCellPopulation);

	/*
	 * Need to check two things:
	 * 1) the map of nodes and their neighbours defines a confluent monolayer
	 * 2) every epithelial index in the tissue population is mapped to a pair of neighbours in the map
	 */

	//Get the epithelial indices
	std::vector<unsigned> epithelial_indices = GetEpithelialIndices(rCellPopulation);

	if ( epithelial_indices.size() != epithelial_nodes_and_their_neighbours.size() )
	{
		EXCEPTION("Not every epithelial cell has been accounted for in defined monolayer");
	}


	//Iterate over each node and its neighbours
	for (std::map<unsigned, std::pair<unsigned, unsigned> >::iterator map_iter = epithelial_nodes_and_their_neighbours.begin();
			map_iter != epithelial_nodes_and_their_neighbours.end();
			map_iter++)
	{

		//Get the node and its neighbours
		unsigned centre_node_index = map_iter->first;

		//Calculate the forces exerted on the left, centre and right nodes by the basement membrane
		c_vector<double, 2> force_on_centre_node = CalculateForceDueToBasementMembrane(rCellPopulation, epithelial_nodes_and_their_neighbours, centre_node_index);

		rCellPopulation.GetNode(centre_node_index)->AddAppliedForceContribution(force_on_centre_node);

	}
}

void OverlappingSpheresBasedBasementMembraneForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n" ;
	*rParamsFile <<  "\t\t\t<TargetCurvature>" << mTargetCurvature << "</TargetCurvature> \n";
	*rParamsFile <<  "\t\t\t<LeftBoundary>"<<  mLeftBoundary << "</LeftBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<RightBoundary>"<<  mRightBoundary << "</RightBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<UsePositionDependentMembraneForce>"<<  mUsePositionDependentMembraneForce << "</UsePositionDependentMembraneForce> \n" ;
	*rParamsFile <<  "\t\t\t<MembraneForceMultiplier>"<<  mMembraneForceMultiplier << "</MembraneForceMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<ApplyPeriodicForce>" << mApplyPeriodicForce << "</ApplyPeriodicForce> \n";

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OverlappingSpheresBasedBasementMembraneForce)
