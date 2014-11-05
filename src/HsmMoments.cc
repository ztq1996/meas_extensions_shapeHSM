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


/// Class to measure HSM adaptive moments
class HsmMoments : public algorithms::ShapeAlgorithm {
public:
    /// @brief Initialize with standard field names and customized documentation.
    HsmMoments(HsmMomentsControl const & ctrl, afw::table::Schema & schema, char const* doc) :
        algorithms::ShapeAlgorithm(ctrl, schema, doc),
        _centroidKeys(addCentroidFields(schema, ctrl.name + ".centroid",
                                        "centroid measured with HSM adaptive moment shape algorithm"))
        {}

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(HsmMoments);

    afw::table::KeyTuple<afw::table::Centroid> _centroidKeys;
};


template<typename PixelT>
void HsmMoments::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // bad until we are good

    HsmMomentsControl const& ctrl = static_cast<HsmMomentsControl const &>(getControl());
    std::vector<std::string> const & badMaskPlanes = ctrl.badMaskPlanes;
    afw::image::MaskPixel badPixelMask = exposure.getMaskedImage().getMask()->getPlaneBitMask(badMaskPlanes);

    afw::geom::Box2I bbox = source.getFootprint()->getBBox();
    if (bbox.getArea() == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "No pixels to measure.");
    }

    ImageConverter<PixelT> const image(exposure.getMaskedImage().getImage(), bbox);

    PTR(afw::image::Image<int>) hsmMask = convertMask(*exposure.getMaskedImage().getMask(), bbox,
                                                      badPixelMask);
    ImageConverter<int> const mask(hsmMask);

    double const psfSigma = exposure.getPsf()->computeShape(center).getTraceRadius();

    galsim::hsm::CppShapeData shape;
    try {
        // GalSim's HSM uses the FITS convention of 1,1 for the lower-left corner
        shape = galsim::hsm::FindAdaptiveMomView(image.getImageView(), mask.getImageView(),
                                                 2.5*psfSigma, 1.0e-6,
                                                 center.getX() - bbox.getMinX() + 1,
                                                 center.getY() - bbox.getMinY() + 1);
    } catch (galsim::hsm::HSMError const& e) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, e.what());
    }

    afw::geom::ellipses::DeterminantRadius const radius(shape.moments_sigma);
    afw::geom::ellipses::Distortion const ellip(shape.observed_shape.getE1(), shape.observed_shape.getE2());
    typedef afw::geom::ellipses::Separable<afw::geom::ellipses::Distortion,
                                           afw::geom::ellipses::DeterminantRadius> Ellipse;

    source.set(_centroidKeys.meas, afw::geom::Point2D(shape.moments_centroid.x + bbox.getMinX() - 1.0,
                                                      shape.moments_centroid.y + bbox.getMinY() - 1.0));
    source.set(getKeys().meas, Ellipse(ellip, radius));
    // XXX calculate errors in shape, centroid?

    source.set(getKeys().flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(HsmMoments);

} // anonymous namespace

PTR(algorithms::AlgorithmControl) HsmMomentsControl::_clone() const {
    return boost::make_shared<HsmMomentsControl>(*this);
}

PTR(algorithms::Algorithm) HsmMomentsControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmMoments>(*this, boost::ref(schema), "adaptive moments from HSM");
}

}}}} // namespace lsst::meas::extensions::shapeHSM
