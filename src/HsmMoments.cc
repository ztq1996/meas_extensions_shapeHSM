// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST
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

template<typename PixelT>
void HsmMomentsAlgorithm::calculate(
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

    galsim::hsm::CppShapeData shape;
    try {
        // GalSim's HSM uses the FITS convention of 1,1 for the lower-left corner
        shape = galsim::hsm::FindAdaptiveMomView(image.getImageView(), mask.getImageView(),
                                                 width, 1.0e-6, center.getX(), center.getY());
    } catch (galsim::hsm::HSMError const& e) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, e.what());
    }

    afw::geom::ellipses::DeterminantRadius const radius(shape.moments_sigma);
    afw::geom::ellipses::Distortion const ellip(shape.observed_shape.getE1(), shape.observed_shape.getE2());
    typedef afw::geom::ellipses::Separable<afw::geom::ellipses::Distortion,
                                           afw::geom::ellipses::DeterminantRadius> Ellipse;

    base::CentroidResult centroidResult;
    centroidResult.x = shape.moments_centroid.x;
    centroidResult.y = shape.moments_centroid.y;
    source.set(_centroidResultKey, centroidResult);
    base::ShapeResult shapeResult;
    shapeResult.setShape(Ellipse(ellip, radius));
    source.set(_momentsKey, shapeResult);
    // XXX calculate errors in shape, centroid?

}

void HsmSourceMomentsAlgorithm::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {

    afw::geom::Point2D center = _centroidExtractor(source, _flagHandler);

    std::vector<std::string> const & badMaskPlanes = _ctrl.badMaskPlanes;
    afw::image::MaskPixel badPixelMask = exposure.getMaskedImage().getMask()->getPlaneBitMask(badMaskPlanes);

    afw::geom::Box2I bbox = source.getFootprint()->getBBox();
    if (bbox.getArea() == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "No pixels to measure.");
    }
    if (!bbox.contains(afw::geom::Point2I(center))) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError,
                          "Center not contained in footprint bounding box");
    }

    double const psfSigma = exposure.getPsf()->computeShape(center).getTraceRadius();

    HsmMomentsAlgorithm::calculate(source, exposure.getMaskedImage().getImage(), exposure.getMaskedImage().getMask(),
                          bbox, center, badPixelMask, 2.5*psfSigma);
}

void HsmSourceMomentsAlgorithm::fail(afw::table::SourceRecord & measRecord, base::MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}

void HsmPsfMomentsAlgorithm::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {

    afw::geom::Point2D center = _centroidExtractor(source, _flagHandler);

    typedef afw::image::Mask<afw::image::MaskPixel> Mask;
    PTR(afw::detection::Psf::Image) image = exposure.getPsf()->computeKernelImage(center);

    // Create a dummy mask
    PTR(Mask) mask = boost::make_shared<Mask>(image->getDimensions());
    *mask = 0;
    mask->setXY0(image->getXY0());

    double const psfSigma = exposure.getPsf()->computeShape(center).getTraceRadius();
    HsmMomentsAlgorithm::calculate(source, image, mask, image->getBBox(afw::image::PARENT), afw::geom::Point2D(0, 0),
                          0, psfSigma);
}

void HsmPsfMomentsAlgorithm::fail(afw::table::SourceRecord & measRecord, base::MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}
}}}} // namespace lsst::meas::extensions::shapeHSM
