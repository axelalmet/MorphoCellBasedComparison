#ifndef VERTICALCOMPRESSIONFORCE_HPP_
#define VERTICALCOMPRESSIONFORCE_HPP_

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

class VerticalCompressionForce : public AbstractForce<2>
{
    friend class TestForces;

private :

    /** Force magnitude */
    double mForceMagnitude;;

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
        archive & mForceMagnitude;

    }

public :

    /**
     * Constructor.
     */
	VerticalCompressionForce();

    /**
     * Destructor.
     */
    ~VerticalCompressionForce();

    void SetForceMagnitude(double forceMagnitude);

    double GetForceMagnitude();

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
CHASTE_CLASS_EXPORT(VerticalCompressionForce)

#endif /*VERTICALCOMPRESSIONFORCE_HPP_*/
