// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

class HsmMomentsControl {
};

class HsmSourceMomentsControl : public HsmMomentsControl {
public:
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>, "Mask planes used to reject bad pixels.");

    HsmSourceMomentsControl(HsmSourceMomentsControl const & other) :
        badMaskPlanes(other.badMaskPlanes) {}

    HsmSourceMomentsControl() {
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("INTRP");
    }
};


class HsmPsfMomentsControl : public HsmMomentsControl {
public:

    HsmPsfMomentsControl() {}

    HsmPsfMomentsControl(HsmPsfMomentsControl const & other) {}
};

/// Base class to measure HSM adaptive moments
///
/// Use this to consolidate common code for HsmSourceMoments and HsmPsfMoments
class HsmMomentsAlgorithm : public base::SimpleAlgorithm {
public:
    /// A typedef to the Control object for this algorithm, defined above.
    ///  The control object contains the configuration parameters for this algorithm.
    typedef HsmMomentsControl Control;

    /// @brief Initialize with standard field names and customized documentation.
    HsmMomentsAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema, char const* doc) :
    _doc(doc),
    _centroidResultKey(
        base::CentroidResultKey::addFields(schema, name, "HSM Centroid", base::NO_UNCERTAINTY)),
    _momentsKey(
        base::ShapeResultKey::addFields(schema, name, "HSM moments", base::NO_UNCERTAINTY)
    ),
    _centroidExtractor(schema, name)
    {}

    /// Calculate moments
    template<typename PixelT>
    void calculate(
        afw::table::SourceRecord& source, // Source for recording moments
        PTR(afw::image::Image<PixelT>) const& afwImage, // Image on which to measure moments
        PTR(afw::image::Mask<afw::image::MaskPixel>) const& afwMask, // Mask for image
        afw::geom::Box2I const& bbox,     // Bounding box
        afw::geom::Point2D const& center, // Starting center for measuring moments
        afw::image::MaskPixel const badPixelMask, // Bitmask for bad pixels
        float const width            // PSF width estimate, for starting moments
        ) const;

protected:
    Control _ctrl;
    std::string _doc;
    base::CentroidResultKey _centroidResultKey;
    base::ShapeResultKey _momentsKey;
    base::FlagHandler _flagHandler;
    base::SafeCentroidExtractor _centroidExtractor;
};


/// Class to measure HSM adaptive moments of source
class HsmSourceMomentsAlgorithm : public HsmMomentsAlgorithm {
public:
    /// A typedef to the Control object for this algorithm, defined above.
    ///  The control object contains the configuration parameters for this algorithm.
    typedef HsmSourceMomentsControl Control;

    /// @brief Initialize with standard field names and customized documentation.
    HsmSourceMomentsAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmMomentsAlgorithm(ctrl, name, schema, "Source adaptive moments algorithm from HSM"), _ctrl(ctrl) {}

private:
    Control _ctrl;
    // These methods are defined private to prevent shadowing of the parent classes
    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    void fail(
        afw::table::SourceRecord & measRecord,
        meas::base::MeasurementError * error=NULL
    ) const;
};

/// Class to measure HSM adaptive moments of PSF
class HsmPsfMomentsAlgorithm : public HsmMomentsAlgorithm {
public:
    /// A typedef to the Control object for this algorithm, defined above.
    ///  The control object contains the configuration parameters for this algorithm.
    typedef HsmPsfMomentsControl Control;

    /// @brief Initialize with standard field names and customized documentation.
    HsmPsfMomentsAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmMomentsAlgorithm(ctrl, name, schema, "Psf adaptive moments algorithm from HSM"), _ctrl(ctrl) {}

private:
    Control _ctrl;
    // These methods are defined private to prevent shadowing of the parent classes
    void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    void fail(
        afw::table::SourceRecord & measRecord,
        meas::base::MeasurementError * error=NULL
    ) const;
};



}}}} // namespace lsst::meas::extensions::shapeHSM

#endif // !LSST_MEAS_EXTENSIONS_SHAPEHSM_HsmMomentsControl_h_INCLUDED
