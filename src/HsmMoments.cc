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
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/math/Random.h"

#include "galsim/Image.h"
#include "galsim/hsm/PSFCorr.h"

#include "lsst/meas/extensions/shapeHSM/HsmAdapter.h"
#include "lsst/meas/extensions/shapeHSM/HsmMomentsControl.h"

namespace {
lsst::meas::base::FlagDefinitionList flagDefinitions;
} // end anonymous

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

base::FlagDefinition const HsmMomentsAlgorithm::FAILURE = flagDefinitions.addFailureFlag("general failure flag, set if anything went wrong");
base::FlagDefinition const HsmMomentsAlgorithm::NO_PIXELS = flagDefinitions.add("flag_no_pixels", "no pixels to measure");
base::FlagDefinition const HsmMomentsAlgorithm::NOT_CONTAINED = flagDefinitions.add("flag_not_contained", "center not contained in footprint bounding box");
base::FlagDefinition const HsmMomentsAlgorithm::PARENT_SOURCE = flagDefinitions.add("flag_parent_source", "parent source, ignored");
base::FlagDefinition const HsmMomentsAlgorithm::GALSIM = flagDefinitions.add("flag_galsim", "GalSim failure");
base::FlagDefinition const HsmPsfMomentsDebiasedAlgorithm::EDGE = flagDefinitions.add("flag_edge", "Variance undefined outside image edge");

base::FlagDefinitionList const & HsmMomentsAlgorithm::getFlagDefinitions() {
    return flagDefinitions;
}

template<typename PixelT>
void HsmMomentsAlgorithm::calculate(
    afw::table::SourceRecord& source,
    std::shared_ptr<afw::image::Image<PixelT>> const& afwImage,
    std::shared_ptr<afw::image::Mask<afw::image::MaskPixel>> const& afwMask,
    geom::Box2I const& bbox,
    geom::Point2D const& center,
    afw::image::MaskPixel const badPixelMask,
    float const width,
    bool roundMoments,
    bool addFlux,
    bool subtractCenter
) const {
    ImageConverter<PixelT> const image(afwImage, bbox);
    std::shared_ptr<afw::image::Image<int>> hsmMask = convertMask(*afwMask, bbox, badPixelMask);
    ImageConverter<int> const mask(hsmMask, bbox);

    galsim::hsm::ShapeData shape;
    try {
        // GalSim's HSM uses the FITS convention of 1,1 for the lower-left corner
        galsim::hsm::FindAdaptiveMomView(shape,
                                         image.getImageView(), mask.getImageView(),
                                         width, 1.0e-6,
                                         galsim::Position<double>(center.getX(), center.getY()), roundMoments);
    } catch (galsim::hsm::HSMError const& e) {
        throw LSST_EXCEPT(base::MeasurementError, e.what(), GALSIM.number);
    }

    afw::geom::ellipses::DeterminantRadius const radius(shape.moments_sigma);
    afw::geom::ellipses::Distortion const ellip(shape.observed_e1, shape.observed_e2);
    typedef afw::geom::ellipses::Separable<afw::geom::ellipses::Distortion,
                                           afw::geom::ellipses::DeterminantRadius> Ellipse;

    base::CentroidResult centroidResult;
    centroidResult.x = shape.moments_centroid.x;
    centroidResult.y = shape.moments_centroid.y;
    if (subtractCenter) {
        centroidResult.x -= center.getX();
        centroidResult.y -= center.getY();
    }
    source.set(_centroidResultKey, centroidResult);
    base::ShapeResult shapeResult;
    shapeResult.setShape(Ellipse(ellip, radius));
    source.set(_momentsKey, shapeResult);
    if (addFlux) {
        assert(_fluxKey.isValid());
        source.set(_fluxKey, shape.moments_amp);
    }
    // XXX calculate errors in shape, centroid?

}

void HsmMomentsAlgorithm::fail(afw::table::SourceRecord & measRecord,
                               base::MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}


void HsmSourceMomentsAlgorithm::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {

    geom::Point2D center = _centroidExtractor(source, _flagHandler);

    std::vector<std::string> const & badMaskPlanes = _ctrl.badMaskPlanes;
    afw::image::MaskPixel badPixelMask = exposure.getMaskedImage().getMask()->getPlaneBitMask(badMaskPlanes);

    geom::Box2I bbox = source.getFootprint()->getBBox();
    if (bbox.getArea() == 0) {
        throw LSST_EXCEPT(
            base::MeasurementError,
            NO_PIXELS.doc,
            NO_PIXELS.number
        );
    }
    if (!bbox.contains(geom::Point2I(center))) {
        throw LSST_EXCEPT(
            base::MeasurementError,
            NOT_CONTAINED.doc,
            NOT_CONTAINED.number
        );
    }

    double const psfSigma = exposure.getPsf()->computeShape(center).getTraceRadius();

    HsmMomentsAlgorithm::calculate(source, exposure.getMaskedImage().getImage(),
                                   exposure.getMaskedImage().getMask(),
                                   bbox, center, badPixelMask, 2.5*psfSigma, _ctrl.roundMoments,
                                   _ctrl.addFlux);
}


// PsfMomentsAlgorithm


void HsmPsfMomentsAlgorithm::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {
    geom::Point2D center = _centroidExtractor(source, _flagHandler);

    auto psfImage = getPsfImage(center, source, exposure);
    auto psfMask = getPsfMask(center, source, exposure);
    auto badPixelMask = getBadPixelMask(exposure);

    double const sigmaGuess = exposure.getPsf()->computeShape(center).getTraceRadius();
    const geom::Point2D centroidGuess = (
        _ctrl->useSourceCentroidOffset ?
        center :
        geom::Point2D(std::floor(center.getX()+0.5), std::floor(center.getY()+0.5))
    );

    bool roundMoments = false;
    bool addFlux = false;
    bool subtractCenter = true;
    HsmMomentsAlgorithm::calculate(
        source, psfImage, psfMask, psfImage->getBBox(afw::image::PARENT),
        centroidGuess, badPixelMask, sigmaGuess, roundMoments, addFlux, subtractCenter
    );
}

std::shared_ptr<afw::detection::Psf::Image> HsmPsfMomentsAlgorithm::getPsfImage(
    geom::Point2D center,
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {
    if (_ctrl->useSourceCentroidOffset) {
        return exposure.getPsf()->computeImage(center);
    } else {
        auto psfImage = exposure.getPsf()->computeKernelImage(center);
        psfImage->setXY0(
            psfImage->getX0() + std::floor(center.getX() + 0.5),
            psfImage->getY0() + std::floor(center.getY() + 0.5)
        );
        return psfImage;
    }
}

std::shared_ptr<afw::image::Mask<afw::image::MaskPixel>> HsmPsfMomentsAlgorithm::getPsfMask(
    geom::Point2D center,
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {
    auto psf = exposure.getPsf();
    auto bbox = psf->computeBBox();

    // Add floor(center+0.5) to bbox to make it look like psf->computeImage()->getBBox()
    auto shift = geom::Extent2I(std::floor(center.getX()+0.5), std::floor(center.getY()+0.5));
    bbox.shift(shift);

    auto mask = std::make_shared<afw::image::Mask<afw::image::MaskPixel>>(bbox);
    *mask = 0;
    return mask;
}


// PsfMomentsDebiasedAlgorithm


std::shared_ptr<afw::detection::Psf::Image> HsmPsfMomentsDebiasedAlgorithm::getPsfImage(
    geom::Point2D center,
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {

    using PsfPixel = afw::detection::Psf::Pixel;

    auto thisCtrl = std::static_pointer_cast<Control::element_type>(_ctrl);

    auto psfImage = HsmPsfMomentsAlgorithm::getPsfImage(center, source, exposure);

    // Psf image crossing exposure edge is fine if we're getting the variance from metadata,
    // but not okay if we're getting the variance from the variance plane.  In both cases,
    // set the EDGE flag, but only fail hard if using variance plane.
    auto overlap = psfImage->getBBox();
    overlap.clip(exposure.getBBox());
    if (overlap != psfImage->getBBox()) {
        _flagHandler.setValue(source, EDGE.number, true);
        if (thisCtrl->noiseSource == "variance") {
            _flagHandler.setValue(source, FAILURE.number, true);
            throw LSST_EXCEPT(
                base::MeasurementError,
                "Variance undefined outside image bounds",
                EDGE.number
            );
        }
    }

    // match Psf flux to source
    auto flux = source.getPsfInstFlux();
    (*psfImage) *= flux;

    // Add Gaussian noise to image
    afw::image::Image<PsfPixel> noise(psfImage->getBBox());
    afw::math::Random rand(afw::math::Random::MT19937, source.getId() + thisCtrl->seedOffset);
    afw::math::randomGaussianImage<afw::image::Image<PsfPixel>>(&noise, rand);
    if (thisCtrl->noiseSource == "meta") {
        double bgmean;
        try {
            bgmean = exposure.getMetadata()->getAsDouble("BGMEAN");
        } catch (pex::exceptions::NotFoundError& e) {
            throw LSST_EXCEPT(base::FatalAlgorithmError, e.what());
        }
        noise *= std::sqrt(bgmean);
    } else if (thisCtrl->noiseSource == "variance") {
        afw::image::Image<float> var(
            *exposure.getMaskedImage().getVariance(),
            psfImage->getBBox(),
            afw::image::PARENT,
            true
        );
        var.sqrt();
        noise *= var;
    }
    (*psfImage) += noise;

    return psfImage;
}

std::shared_ptr<afw::image::Mask<afw::image::MaskPixel>> HsmPsfMomentsDebiasedAlgorithm::getPsfMask(
    geom::Point2D center,
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {
    // Copy part of exposure mask that overlaps psf bbox.
    // Assume rest of PSF is unmasked.
    auto psfMask = HsmPsfMomentsAlgorithm::getPsfMask(center, source, exposure);
    auto overlap = psfMask->getBBox();
    overlap.clip(exposure.getBBox());
    (*psfMask)[overlap] = (*exposure.getMaskedImage().getMask())[overlap];
    return psfMask;
}

afw::image::MaskPixel const HsmPsfMomentsDebiasedAlgorithm::getBadPixelMask(
    afw::image::Exposure<float> const & exposure
) const {
    return exposure.getMaskedImage().getMask()->getPlaneBitMask(
        std::static_pointer_cast<Control::element_type>(_ctrl)->badMaskPlanes
    );
}

}}}} // namespace lsst::meas::extensions::shapeHSM
