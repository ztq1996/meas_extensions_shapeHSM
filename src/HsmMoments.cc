// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/algorithms/Measure.h"

#include "galsim/Image.h"
#include "galsim/hsm/PSFCorr.h"

#include "lsst/meas/extensions/shapeHSM/HsmAdapter.h"
#include "lsst/meas/extensions/shapeHSM/HsmMomentsControl.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

namespace {


/// Base class to measure HSM adaptive moments
///
/// Use this to consolidate common code for HsmSourceMoments and HsmPsfMoments
class HsmMoments : public algorithms::ShapeAlgorithm {
public:
    /// @brief Initialize with standard field names and customized documentation.
    template <typename ControlT>
    HsmMoments(ControlT const & ctrl, afw::table::Schema & schema, char const* doc) :
        algorithms::ShapeAlgorithm(ctrl, schema, doc),
        _centroidKeys(addCentroidFields(schema, ctrl.name + ".centroid",
                                        "centroid measured with HSM adaptive moment shape algorithm"))
        {}

    /// Calculate moments
    template <typename PixelT>
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
    afw::table::KeyTuple<afw::table::Centroid> _centroidKeys;
};


/// Class to measure HSM adaptive moments of source
class HsmSourceMoments : public HsmMoments {
public:
    /// @brief Initialize with standard field names and customized documentation.
    HsmSourceMoments(HsmSourceMomentsControl const & ctrl, afw::table::Schema & schema, char const* doc) :
        HsmMoments(ctrl, schema, doc) {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(HsmSourceMoments);

};

/// Class to measure HSM adaptive moments of PSF
class HsmPsfMoments : public HsmMoments {
public:
    /// @brief Initialize with standard field names and customized documentation.
    HsmPsfMoments(HsmPsfMomentsControl const & ctrl, afw::table::Schema & schema, char const* doc) :
        HsmMoments(ctrl, schema, doc) {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(HsmPsfMoments);

};


template <typename PixelT>
void HsmMoments::calculate(
    afw::table::SourceRecord& source,
    PTR(afw::image::Image<PixelT>) const& afwImage,
    PTR(afw::image::Mask<afw::image::MaskPixel>) const& afwMask,
    afw::geom::Box2I const& bbox,
    afw::geom::Point2D const& center,
    afw::image::MaskPixel const badPixelMask,
    float const width
) const {
    ImageConverter<PixelT> const image(afwImage, bbox);
    PTR(afw::image::Image<int>) hsmMask = convertMask(*afwMask, bbox, badPixelMask);
    ImageConverter<int> const mask(hsmMask);

    int const x0 = afwImage->getX0() - 1, y0 = afwImage->getY0() - 1;

    galsim::hsm::CppShapeData shape;
    try {
        // GalSim's HSM uses the FITS convention of 1,1 for the lower-left corner
        shape = galsim::hsm::FindAdaptiveMomView(image.getImageView(), mask.getImageView(),
                                                 width, 1.0e-6, center.getX() - x0, center.getY() - y0);
    } catch (galsim::hsm::HSMError const& e) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, e.what());
    }

    afw::geom::ellipses::DeterminantRadius const radius(shape.moments_sigma);
    afw::geom::ellipses::Distortion const ellip(shape.observed_shape.getE1(), shape.observed_shape.getE2());
    typedef afw::geom::ellipses::Separable<afw::geom::ellipses::Distortion,
                                           afw::geom::ellipses::DeterminantRadius> Ellipse;

    source.set(_centroidKeys.meas, afw::geom::Point2D(shape.moments_centroid.x + x0,
                                                      shape.moments_centroid.y + y0));
    source.set(getKeys().meas, Ellipse(ellip, radius));
    // XXX calculate errors in shape, centroid?

    source.set(getKeys().flag, false);
}



template<typename PixelT>
void HsmSourceMoments::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // bad until we are good

    HsmSourceMomentsControl const& ctrl = static_cast<HsmSourceMomentsControl const &>(getControl());
    std::vector<std::string> const & badMaskPlanes = ctrl.badMaskPlanes;
    afw::image::MaskPixel badPixelMask = exposure.getMaskedImage().getMask()->getPlaneBitMask(badMaskPlanes);

    afw::geom::Box2I bbox = source.getFootprint()->getBBox();
    if (bbox.getArea() == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "No pixels to measure.");
    }

    double const psfSigma = exposure.getPsf()->computeShape(center).getTraceRadius();

    HsmMoments::calculate(source, exposure.getMaskedImage().getImage(), exposure.getMaskedImage().getMask(),
                          bbox, center, badPixelMask, 2.5*psfSigma);
}


template<typename PixelT>
void HsmPsfMoments::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // bad until we are good

    typedef afw::image::Mask<afw::image::MaskPixel> Mask;
    PTR(afw::detection::Psf::Image) image = exposure.getPsf()->computeKernelImage(center);

    // Create a dummy mask
    PTR(Mask) mask = boost::make_shared<Mask>(image->getDimensions());
    *mask = 0;
    mask->setXY0(image->getXY0());

    double const psfSigma = exposure.getPsf()->computeShape(center).getTraceRadius();
    HsmMoments::calculate(source, image, mask, image->getBBox(afw::image::PARENT), afw::geom::Point2D(0, 0),
                          0, psfSigma);

    // Adjust centroid for size of box
    source.set(_centroidKeys.meas, source.get(_centroidKeys.meas) - afw::geom::Extent2D(image->getXY0()));
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(HsmSourceMoments);
LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(HsmPsfMoments);

} // anonymous namespace

PTR(algorithms::AlgorithmControl) HsmSourceMomentsControl::_clone() const {
    return boost::make_shared<HsmSourceMomentsControl>(*this);
}

PTR(algorithms::Algorithm) HsmSourceMomentsControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmSourceMoments>(*this, boost::ref(schema), "source adaptive moments from HSM");
}

PTR(algorithms::AlgorithmControl) HsmPsfMomentsControl::_clone() const {
    return boost::make_shared<HsmPsfMomentsControl>(*this);
}

PTR(algorithms::Algorithm) HsmPsfMomentsControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmPsfMoments>(*this, boost::ref(schema), "PSF adaptive moments from HSM");
}

}}}} // namespace lsst::meas::extensions::shapeHSM
