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
#include "lsst/meas/algorithms/Measure.h"

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
    name.push_back('.');
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
addErrorField(std::string const & name, MeasType measType, afw::table::Schema & schema) {
    std::string doc = "Uncertainty on ";
    doc.push_back(measType);
    doc += "1 and ";
    doc.push_back(measTypeSymbol(measType));
    doc += "2 (assumed to be the same)";
    return schema.addField<double>(name + ".err", doc);
}

/*
 *  HSM shape algorithm class - we use one class for all algorithms; all the specialized
 *  work is done by the control class _makeAlgorithm implementations.
 *
 *  We don't inherit from ShapeAlgorithm because none of these compute any errors,
 *  and we don't want to add fields that are always NaNs.  We probably need to think
 *  about how to handle the "is-a" ShapeAlgorithm test a little better, but we don't
 *  rely on it anywhere presently.
 */
class HsmShape : public algorithms::Algorithm {
public:

    /// @brief Initialize with standard field names and customized documentation.
    HsmShape(
        HsmShapeControl const & ctrl,
        std::string const & shearType,
        MeasType measType,
        afw::table::Schema & schema,
        std::string const & doc
    );

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(HsmShape);

    std::string _shearType;
    MeasType _measType;
    afw::table::Key<double> _e1Key;
    afw::table::Key<double> _e2Key;
    afw::table::Key<double> _sigmaKey;
    afw::table::Key<double> _resolutionKey;
    afw::table::Key<afw::table::Flag> _flagKey;
    bool _hasDeblendKey;
    afw::table::Key<int> _deblendKey;
};

HsmShape::HsmShape(
    HsmShapeControl const & ctrl, 
    std::string const & shearType,
    MeasType measType,
    afw::table::Schema & schema,
    std::string const & doc
) : Algorithm(ctrl),
    _shearType(shearType),
    _measType(measType),
    _e1Key(addEllipticityField(ctrl.name, '1', measType, schema, doc)),
    _e2Key(addEllipticityField(ctrl.name, '2', measType, schema, doc)),
    _sigmaKey(
        schema.addField<double>(ctrl.name + ".sigma", doc + " (width)")
    ),
    _resolutionKey(
        schema.addField<double>(ctrl.name + ".resolution", "resolution factor (0=unresolved, 1=resolved)")
    ),
    _flagKey(schema.addField<afw::table::Flag>(ctrl.name + ".flags", "set if measurement failed in any way")),
    _hasDeblendKey(ctrl.deblendNChild.size() > 0)
{
    if (_hasDeblendKey) {
        _deblendKey = schema[ctrl.deblendNChild];
    }
}

template<typename PixelT>
void HsmShape::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    typedef afw::image::Exposure<PixelT> ExposureT;

    source.set(_flagKey, true); // bad until we are good

    if (_hasDeblendKey && source.get(_deblendKey) > 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, "Ignoring parent source");
    }

    HsmShapeControl const& ctrl = static_cast<HsmShapeControl const &>(getControl());
    std::vector<std::string> const & badMaskPlanes = ctrl.badMaskPlanes;

    afw::image::MaskPixel badPixelMask = exposure.getMaskedImage().getMask()->getPlaneBitMask(badMaskPlanes);

    afw::geom::Box2I bbox = source.getFootprint()->getBBox();
    if (bbox.getArea() == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "No pixels to measure.");
    }

    PTR(afw::detection::Psf::Image) psf = exposure.getPsf()->computeImage(center);
    psf->setXY0(0, 0);
    ImageConverter<PixelT> const image(exposure.getMaskedImage().getImage(), bbox);
    ImageConverter<afw::detection::Psf::Image::Pixel> const psfImage(psf);

    typedef afw::image::Image<int> ImageI;
    PTR(typename ExposureT::MaskedImageT::Mask) afwMask = exposure.getMaskedImage().getMask();
    PTR(ImageI) hsmMask = convertMask(*afwMask, bbox, badPixelMask);
    ImageConverter<int> const mask(hsmMask);

    PTR(ImageI) dummyMask = boost::make_shared<ImageI>(psf->getDimensions());
    *dummyMask = 1;
    ImageConverter<int> const psfMask(dummyMask);

    // Calculate the sky variance
    afw::math::StatisticsControl sctrl;
    sctrl.setAndMask(badPixelMask);
    typename ExposureT::MaskedImageT::Variance const variance(*exposure.getMaskedImage().getVariance(), bbox);
    afw::math::Statistics stat = afw::math::makeStatistics(variance, *afwMask, afw::math::MEDIAN, sctrl);
    double const skyvar = sqrt(stat.getValue(afw::math::MEDIAN));

    double const psfSigma = exposure.getPsf()->computeShape(center).getTraceRadius();

    galsim::hsm::CppShapeData shape, psfShape;

    try {
        shape = galsim::hsm::EstimateShearView(image.getImageView(), psfImage.getImageView(),
                                               mask.getImageView(), skyvar, _shearType.c_str(), "FIT",
                                               2.5*psfSigma, psfSigma, 1.0e-6,
                                               center.getX() - bbox.getMinX() + 1,
                                               center.getY() - bbox.getMinY() + 1);
    } catch (galsim::hsm::HSMError const& e) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeError, e.what());
    }

    assert(shape.meas_type[0] == measTypeSymbol(_measType));
    switch (*shape.meas_type.c_str()) {
    case 'e':
        source.set(_e1Key, shape.corrected_e1);
        source.set(_e2Key, shape.corrected_e2);
        source.set(_sigmaKey, 0.5*shape.corrected_shape_err);
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
    source.set(_flagKey, shape.correction_status != 0);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(HsmShape);

} // anonymous namespace

PTR(algorithms::AlgorithmControl) HsmShapeBjControl::_clone() const {
    return boost::make_shared<HsmShapeBjControl>(*this);
}

PTR(algorithms::AlgorithmControl) HsmShapeLinearControl::_clone() const {
    return boost::make_shared<HsmShapeLinearControl>(*this);
}

PTR(algorithms::AlgorithmControl) HsmShapeKsbControl::_clone() const {
    return boost::make_shared<HsmShapeKsbControl>(*this);
}

PTR(algorithms::AlgorithmControl) HsmShapeRegaussControl::_clone() const {
    return boost::make_shared<HsmShapeRegaussControl>(*this);
}

PTR(algorithms::Algorithm) HsmShapeBjControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmShape>(
        *this, "BJ", ELLIPTICITY, boost::ref(schema),
        "PSF-corrected ellipticity using Bernstein & Jarvis (2002) method"
    );
}

PTR(algorithms::Algorithm) HsmShapeLinearControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmShape>(
        *this, "LINEAR", ELLIPTICITY, boost::ref(schema),
        "PSF-corrected ellipticity using Hirata & Seljak (2003) 'linear' method"
    );
}

PTR(algorithms::Algorithm) HsmShapeKsbControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmShape>(
        *this, "KSB", SHEAR, boost::ref(schema),
        "PSF-corrected shear using KSB method"
    );
}

PTR(algorithms::Algorithm) HsmShapeRegaussControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmShape>(
        *this, "REGAUSS", ELLIPTICITY, boost::ref(schema),
        "PSF-corrected ellipticity using Hirata & Seljak (2003) 'regaussianization' method"
    );
}

}}}} // namespace lsst::meas::extensions::shapeHSM
