#ifndef NONPERIODICBASEMENTMEMBRANEFORCE_HPP_
#define NONPERIODICBASEMENTMEMBRANEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"

#include <cmath>
#include <list>
#include <fstream>

/**
 * A force class that defines the force due to the basement membrane.
 */

class NonPeriodicBasementMembraneForce : public AbstractForce<2>
{
    friend class TestCrossSectionModelInteractionForce;

private :

    /** Parameter that defines the stiffness of the basement membrane */
    double mBasementMembraneParameter;

    /** Target curvature for the layer of cells */
    double mTargetCurvature;

    /** Left boundary to specify heterogeneous target curvature */
    double mLeftBoundary;

    /** Right boundary to specify heterogeneous target curvature */
    double mRightBoundary;

    /** Specifies whether or not we are applying heterogeneous target curvature */
    bool mApplyForceToCrypt;

    /** Specifies whether basement membrane stiffness is heterogeneous */
    bool mUsePositionDependentMembraneForce;

    /** Specify whether target curvature is vertically heterogeneous */
    bool mApplyVerticallyDependentTargetCurvature;

    double mMembraneForceMultiplier;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mBasementMembraneParameter;
        archive & mTargetCurvature;
        archive & mLeftBoundary;
        archive & mRightBoundary;
        archive & mApplyForceToCrypt;
        archive & mApplyVerticallyDependentTargetCurvature;
        archive & mMembraneForceMultiplier;

    }

public :

    /**
     * Constructor.
     */
	NonPeriodicBasementMembraneForce();

    /**
     * Destructor.
     */
    ~NonPeriodicBasementMembraneForce();

    /* Set method for Basement Membrane Parameter
     */
    void SetBasementMembraneParameter(double basementMembraneParameter);

    /* Get method for Basement Membrane Parameter
     */
    double GetBasementMembraneParameter();

    /* Value of Target Curvature in epithelial layer */
    void SetTargetCurvature(double targetCurvature = 0.0);

    /* Get method for Target Curvature
     *
     */
    double GetTargetCurvature();

    /*
     * Set method for Left Boundary
     */
    void SetLeftCryptBoundary(double leftBoundary);

    /*
     * Get method for Left Boundary
     */
    double GetLeftCryptBoundary();

    /*
     * Set method for Right Boundary
     */
    void SetRightCryptBoundary(double rightBoundary);

    /*
     * Get method for Right Boundary
     */
    double GetRightCryptBoundary();

    bool IsForceAppliedToCrypt();

    void ApplyForceToCrypt(bool applyForceToCrypt);

    void ApplyVerticallyDependentTargetCurvature(bool applyVerticallyDependentTargetCurvature);

    bool ApplyVerticallyDependentTargetCurvature();

    /* Removing duplicated entries of a vector
     */
    void RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates);

    /* Returns a boolean for whether the element contains ghost nodes
     */
    bool DoesElementContainGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned elementIndex);

    /* Returns crypt height extremes to apply non-zero curvature to base */
    c_vector<double, 2> GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation);

    /* Finding the connected pairs of epithelial-tissue nodes
     */
    std::vector<c_vector<unsigned, 2> > GetEpithelialStromalPairs(AbstractCellPopulation<2>& rCellPopulation);

    /* Takes an epithelial node index and a tissue node index and returns the curvature of
     * the curve passing through the midpoints of the epithelial-tissue springs of the
     * common elements
     */
    double GetCurvatureFromNodePair(AbstractCellPopulation<2>& rCellPopulation, unsigned epithelialNodeIndex,
    		unsigned stromalNodeIndex);

    /*
     * Finding the curvature between three midpoints parametrically - in this case, we find the normal
     * to the vector joining the left and right midpoints, and then find the perpendicular distance of
     * the centre midpoint from the left->right vector
     */
    double GetCurvatureFromMidpoints(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> leftMidpoint,
    														c_vector<double, 2> centreMidpoint,
    														c_vector<double, 2> rightMidpoint);

    double FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
    								c_vector<double, 2> leftMidpoint,
									c_vector<double, 2> centreMidpoint,
									c_vector<double, 2> rightMidpoint);

    /* Finding the number of elements that a node belongs to, which contain only real nodes
     * and not ghost nodes
     */
    unsigned GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /* Finding the y-coordinates of the crypt orifice and crypt base
     * The first entry of the resulting vector is the orifice, the second is the base
     */
//    c_vector<double,2> GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Method to return the nodes connected to a particular node via the Delaunay
     * triangulation, excluding ghost nodes.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /* Method to return a boolean that indicates whether this node/cell has detached from the basement membrane
     */
    bool HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    /**
     * Outputs force Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(NonPeriodicBasementMembraneForce)

#endif /*NonPeriodicBASEMENTMEMBRANEFORCE_HPP_*/
