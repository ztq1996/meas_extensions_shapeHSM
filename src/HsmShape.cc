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

#include "lsst/afw/table/Source.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"

#include "galsim/Image.h"
#include "galsim/hsm/PSFCorr.h"

#include "lsst/meas/extensions/shapeHSM/HsmAdapter.h"
#include "lsst/meas/extensions/shapeHSM/HsmShapeControl.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

namespace {

/// Return appropriate symbol for a measurement type
char measTypeSymbol(MeasType measType)
{
    switch (measType) {
    case ELLIPTICITY:
        return 'e';
    case SHEAR:
        return 'g';
    }
    assert(false);
}

// helper function for HsmShapeAlgorithm ctor
inline afw::table::Key<double>
addEllipticityField(std::string name, char n, MeasType measType, afw::table::Schema & schema,
                    std::string doc) {
    name.push_back('_');
    name.push_back(measTypeSymbol(measType));
    name.push_back(n);
    if (n == '1') {
        doc += " (+";
    } else {
        doc += " (x";
    }
    switch (measType) {
    case ELLIPTICITY:
        doc += " component of ellipticity)";
        break;
    case SHEAR:
        doc += " component of estimated shear)";
        break;
    default:
        assert(false);
    }
    return schema.addField<double>(name, doc);
}

// helper function for HsmShapeAlgorithm ctor
inline afw::table::Key<double>
addErrorField(std::string name, MeasType measType, afw::table::Schema & schema) {
    std::string doc = "Uncertainty on ";
    doc.push_back(measType);
    doc += "1 and ";
    doc.push_back(measTypeSymbol(measType));
    doc += "2 (assumed to be the same)";
    name.push_back('_');
    name.push_back(measTypeSymbol(measType));
    return schema.addField<double>(name + "_sigma", doc);
}
} // end anonymous

HsmShapeAlgorithm::HsmShapeAlgorithm(
    Control const & ctrl,
    std::string const & name,
    std::string const & shearType,
    MeasType measType,
    std::string const & doc,
    afw::table::Schema & schema
) : _ctrl(ctrl),
    _shearType(shearType),
    _measType(measType),
    _e1Key(addEllipticityField(name, '1', measType, schema, doc)),
    _e2Key(addEllipticityField(name, '2', measType, schema, doc)),
    _sigmaKey(
        schema.addField<double>(name + "_sigma", doc + " (width)")
    ),
    _resolutionKey(
        schema.addField<double>(name + "_resolution", "resolution factor (0=unresolved, 1=resolved)")
    ),
    _doc(doc),
    _hasDeblendKey(_ctrl.deblendNChild.size() > 0),
    _centroidExtractor(schema, name)
{
    static boost::array<base::FlagDefinition,N_FLAGS> const flagDefs = {{
        {"flag", "general failure flag"},
        {"flag_no_pixels", "no pixels to measure"},
        {"flag_not_contained", "center not contained in footprint bounding box"},
        {"flag_has_deblend", "object has been deblended"},
        {"flag_ignore_parent_source", "ignoring parent source"},
    }};
    _flagHandler = base::FlagHandler::addFields(schema, name, flagDefs.begin(), flagDefs.end());

    if (_hasDeblendKey) {
        _deblendKey = schema[ctrl.deblendNChild];
    }
}

void HsmShapeAlgorithm::measure(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const & exposure
) const {
    afw::geom::Point2D center = _centroidExtractor(source, _flagHandler);
    if (_hasDeblendKey && source.get(_deblendKey) > 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, "Ignoring parent source");
    }

    std::vector<std::string> const & badMaskPlanes = _ctrl.badMaskPlanes;

    afw::image::MaskPixel badPixelMask = exposure.getMaskedImage().getMask()->getPlaneBitMask(badMaskPlanes);

    afw::geom::Box2I bbox = source.getFootprint()->getBBox();
    if (bbox.getArea() == 0) {
        throw LSST_EXCEPT(
            base::MeasurementError,
            _flagHandler.getDefinition(NO_PIXELS).doc,
            NO_PIXELS
        );
    }
    if (!bbox.contains(afw::geom::Point2I(center))) {
        throw LSST_EXCEPT(
            base::MeasurementError,
            _flagHandler.getDefinition(NOT_CONTAINED).doc,
            NOT_CONTAINED
        );
    }

    PTR(afw::detection::Psf::Image) psf = exposure.getPsf()->computeImage(center);
    psf->setXY0(0, 0);
    ImageConverter<float> const image(exposure.getMaskedImage().getImage(), bbox);
    ImageConverter<afw::detection::Psf::Image::Pixel> const psfImage(psf);

    typedef afw::image::Image<int> ImageI;
    PTR(afw::image::Exposure<float>::MaskedImageT::Mask) afwMask = exposure.getMaskedImage().getMask();
    PTR(ImageI) hsmMask = convertMask(*afwMask, bbox, badPixelMask);
    ImageConverter<int> const mask(hsmMask, bbox);

    PTR(ImageI) dummyMask = boost::make_shared<ImageI>(psf->getDimensions());
    *dummyMask = 1;
    ImageConverter<int> const psfMask(dummyMask);

    afw::math::StatisticsControl sctrl;
    sctrl.setAndMask(badPixelMask);
    afw::image::Exposure<float>::MaskedImageT::Variance const variance(*exposure.getMaskedImage().getVariance(),
                                                                       bbox, afw::image::PARENT);
    afw::math::Statistics stat = afw::math::makeStatistics(variance, *afwMask, afw::math::MEDIAN, sctrl);
    double const skyvar = sqrt(stat.getValue(afw::math::MEDIAN));

    double const psfSigma = exposure.getPsf()->computeShape(center).getTraceRadius();

    galsim::hsm::CppShapeData shape, psfShape;

    try {
        shape = galsim::hsm::EstimateShearView(image.getImageView(), psfImage.getImageView(),
                                               mask.getImageView(), skyvar, _shearType.c_str(), "FIT",
                                               2.5*psfSigma, psfSigma, 1.0e-6, center.getX(), center.getY());
    } catch (galsim::hsm::HSMError const& e) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, e.what());
    }

    assert(shape.meas_type[0] == measTypeSymbol(_measType));
    switch (*shape.meas_type.c_str()) {
    case 'e':
        source.set(_e1Key, shape.corrected_e1);
        source.set(_e2Key, shape.corrected_e2);
        source.set(_sigmaKey, 2.0*shape.corrected_shape_err);
        break;
    case 'g':
        source.set(_e1Key, shape.corrected_g1);
        source.set(_e2Key, shape.corrected_g2);
        source.set(_sigmaKey, shape.corrected_shape_err);
        break;
    default:
        assert(false);
    }
    source.set(_resolutionKey, shape.resolution_factor);
    _flagHandler.setValue(source, FAILURE, shape.correction_status != 0);
}

void HsmShapeAlgorithm::fail(afw::table::SourceRecord & measRecord, base::MeasurementError * error) const {
    _flagHandler.handleFailure(measRecord, error);
}
}}}} // namespace lsst::meas::extensions::shapeHSM
