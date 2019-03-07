#ifndef _HPP_
#define OverlappingSpheresBasedBasementMembraneForce_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "RandomNumberGenerator.hpp"

#include <cmath>
#include <list>
#include <fstream>

/**
 * MODIFIED FOR PERSONAL USE BY AXEL ALMET
 * A force class that defines the force due to the basement membrane.
 */

class OverlappingSpheresBasedBasementMembraneForce : public AbstractForce<2>
{
    friend class TestCrossSectionModelInteractionForce;

private :

    /** Parameter that multiplies the curvature to give the basement membrane force */
    double mBasementMembraneParameter;

    /** Target curvature for the ring of cells (NodeBased) */
    double mTargetCurvature;

    /** x-coordinate that encloses the region in which to apply a non-zero target curvature */
    double mLeftBoundary;

    /** x-coordinate that encloses the region in which to apply a non-zero target curvature */
    double mRightBoundary;

    /** Make the basement membrane force dependent on the position of a cell up the crypt */
    bool mUsePositionDependentMembraneForce;

    /** Boolean to check whether force is applied to crypt-like epithelium or organoid */
    bool mApplyForceToCrypt;

    /** The multiplication factor for the basement membrane parameter */
    double mMembraneForceMultiplier;

    /** The cut off radius for defining neighbouring nodes */
    double mCutOffRadius;

    bool mApplyPeriodicForce;

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
        archive & mUsePositionDependentMembraneForce;
        archive & mMembraneForceMultiplier;
        archive & mCutOffRadius;
        archive & mApplyPeriodicForce;
    }

public :

    /**
     * Constructor.
     */
	OverlappingSpheresBasedBasementMembraneForce();

    /**
     * Destructor.
     */
    ~OverlappingSpheresBasedBasementMembraneForce();

    /* Set method for Basement Membrane Parameter
     */
    void SetBasementMembraneParameter(double basementMembraneParameter);

    /* Get method for Basement Membrane Parameter
     */
    double GetBasementMembraneParameter();

    /* Value of curvature at crypt base */
    void SetTargetCurvature(double targetCurvature = 0.0);

    /*
     * Get method for Target Curvature Parameters
     */
    double GetTargetCurvature();

    /*
     * Set method for left crypt boundary
     */
    void SetLeftCryptBoundary(double leftBoundary);

    /*
     * Get method for Left crypt boundary parameter
     */
    double GetLeftCryptBoundary();

    /*
     * Set method for right crypt boundary parameter
     */
    void SetRightCryptBoundary(double rightBoundary);

    /*
     * Get method for right crypt boundary
     */
    double GetRightCryptBoundary();

    /* Returns crypt height extremes to apply non-zero curvature to base */
    c_vector<double, 2> GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation);

    /* Returns crypt width extremes */
    c_vector<double, 2> GetCryptWidthExtremes(AbstractCellPopulation<2>& rCellPopulation);

    /* Return tangent vector at point based on cosine parametrisation */
    c_vector<double, 2> GetCosineBasedTangentVector(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> point);

    /*
     * Return vector of epithelial indices that are close to the considered epithelial node,
     * but based on an approximated cosine approximation
     */
    std::vector<unsigned> GetClosestNeighboursBasedOnCosineApproximation(AbstractCellPopulation<2>& rCellPopulation, unsigned epithelialIndex);

    /*
     * Return the nearest neighbour based on vector projections from the tangent vector
     * at a point along the cosine approximation of the epithelium.
     */
    unsigned GetNearestNeighboursAlongCosineApproximation(AbstractCellPopulation<2>& rCellPopulation, unsigned epithelialIndex);

    /*
     * Set method for geometry-dependent basement membrane for application, i.e.crypt or organoid
     */
    void SetCryptGeometry(bool applyForceToCrypt = true);

    /*
     * Check to determine whether or not basement membrane force is to be applied to a crypt or organoid
     */
    bool GetCryptGeometryCheck();

    /* Set method for position-dependent basement membrane force multiplier (i.e. if you want to apply a different
     * basement membrane parameter in the crypt base, or at the orifice)
     * @param usePositionDependentMembraneForce whether to multiply the basement membrane force by a factor
     * @param membraneForceMultiplier the multiplication factor for the basement membrane force
     */
    void SetPositionDependentMultiplier(bool usePositionDependentMembraneForce = false, double membraneForceMultiplier = 1.0);

    /* Get method for basement membrane force strength multiplier
     */
    double GetPositionDependentMultiplier();

    /*Get method for cut off radius */
    double GetCutOffRadius();

    /*Set method for cut off radius */
    void SetCutOffRadius(double cutOffRadius);

    /* Get method for whether we apply a periodic force */
    void ApplyPeriodicForce(bool applyPeriodicForce);

    bool IsPeriodicForceApplied();

    /* Removing duplicated entries of a vector
     */
    void RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates);

    double FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
    								c_vector<double, 2> leftPoint,
									c_vector<double, 2> centrePoint,
									c_vector<double, 2> rightPoint);

    /*
     * Method to get the indices of the monolayer
     */
    std::vector<unsigned> GetEpithelialIndices(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Method to return the nodes connected to a particular node within a defined cut-off
     * radius
     */
    std::vector<unsigned> GetNeighbouringEpithelialIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /*
     * Method to get the epithelial nodes and their left and right neighbours, if they have any
     */
    std::map<unsigned, std::pair<unsigned, unsigned> > GetEpithelialIndicesAndTheirLeftAndRightEpithelialNeighbours(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Method to check whether or not node is on the left or right boundary (and hence should have zero force).
     */
    bool IsBoundaryNode(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /*
     * Method to calculate the force due to basement membrane on an epithelial cell
     */
    c_vector<double, 2> CalculateForceDueToBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned > > epithelialIndicesAndNeighbours, unsigned nodeIndex);

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
CHASTE_CLASS_EXPORT(OverlappingSpheresBasedBasementMembraneForce)

#endif /*OVERLAPPINGSPHERESBASEDBASEMENTMEMBRANEFORCE_HPP_*/
