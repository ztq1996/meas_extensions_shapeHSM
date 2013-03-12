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
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/Statistics.h"

#include "lsst/meas/extensions/shapeHSM/HsmShapeControl.h"
#include "hsm/PSFCorr.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

namespace {

// helper function for HsmShapeAlgorithm ctor
inline afw::table::Key<double>
addEllipticityField(std::string name, char n, char measType, afw::table::Schema & schema, std::string doc) {
    name.push_back('.');
    name.push_back(measType);
    name.push_back(n);
    if (n == '1') {
        doc += " (+";
    } else {
        doc += " (x";
    }
    if (measType == 'e') {
        doc += " component of ellipticity)";
    } else {
        doc += " component of estimated shear)";
    }
    return schema.addField<double>(name, doc);
}

// helper function for HsmShapeAlgorithm ctor
inline afw::table::Key<double>
addErrorField(std::string const & name, char measType, afw::table::Schema & schema) {
    std::string doc = "Uncertainty on ";
    doc.push_back(measType);
    doc += "1 and ";
    doc.push_back(measType);
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
        char measType,
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
    afw::table::Centroid::MeasKey _centroidKey;
    afw::table::Shape::MeasKey _momentsKey;
#if 0
    afw::table::Shape::MeasKey _psfMomentsKey;
#endif
    afw::table::Key<double> _e1Key;
    afw::table::Key<double> _e2Key;
    afw::table::Key<double> _errKey;
    afw::table::Key<double> _resolutionKey;
    afw::table::Key<afw::table::Flag> _flagKey;
};

HsmShape::HsmShape(
    HsmShapeControl const & ctrl, 
    std::string const & shearType,
    char measType,
    afw::table::Schema & schema,
    std::string const & doc
) : Algorithm(ctrl),
    _shearType(shearType),
    _centroidKey(
        schema.addField<afw::table::Centroid::MeasTag>(ctrl.name + ".centroid", doc + " (centroid)")
    ),
    _momentsKey(
        schema.addField<afw::table::Shape::MeasTag>(ctrl.name + ".moments", doc + " (uncorrected moments)")
    ),
#if 0
    _psfMomentsKey(schema.addField<afw::table::Shape::MeasTag>(ctrl.name + ".psf", doc + " (PSF moments)")),
#endif
    _e1Key(addEllipticityField(ctrl.name, '1', measType, schema, doc)),
    _e2Key(addEllipticityField(ctrl.name, '2', measType, schema, doc)),
    _errKey(addErrorField(ctrl.name, measType, schema)),
    _resolutionKey(
        schema.addField<double>(ctrl.name + ".resolution", "resolution factor (0=unresolved, 1=resolved)")
    ),
    _flagKey(schema.addField<afw::table::Flag>(ctrl.name + ".flags", "set if measurement failed in any way"))
{}

template<typename PixelT>
void HsmShape::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_flagKey, true); // bad until we are good
    std::vector<std::string> const & badMaskPlanes 
        = static_cast<HsmShapeControl const &>(getControl()).badMaskPlanes;

    afw::image::MaskPixel badPixelMask 
        = exposure.getMaskedImage().getMask()->getPlaneBitMask(badMaskPlanes);

    afw::geom::Box2I bbox(source.getFootprint()->getBBox());

    afw::image::MaskedImage<PixelT> subMaskedImage(exposure.getMaskedImage(), bbox, afw::image::PARENT);
    afw::image::Image<int> goodPixels(bbox);
    for (int iY = 0; iY < bbox.getWidth(); ++iY) {
        afw::image::Mask<>::x_iterator pIn = subMaskedImage.getMask()->row_begin(iY);
        afw::image::Image<int>::x_iterator pOut = goodPixels.row_begin(iY);
        afw::image::Image<int>::x_iterator const pOutEnd = goodPixels.row_end(iY);
        for (; pOut != pOutEnd; ++pOut, ++pIn) {
            *pOut = (badPixelMask & *pIn) ? 0 : 1;
        }
    }

    afw::math::StatisticsControl sctrl;
    sctrl.setAndMask(badPixelMask);
    afw::math::Statistics stat = afw::math::makeStatistics(*subMaskedImage.getVariance(),
                                                           *subMaskedImage.getMask(),
                                                           afw::math::MEDIAN, sctrl);
    float skyvar = std::sqrt(stat.getValue(afw::math::MEDIAN));

    try {

        galsim::hsm::CppHSMShapeData data = galsim::hsm::EstimateShearHSMView(
            galsim::makeImageView(*subMaskedImage.getImage()),
            galsim::makeImageView(*exposure.getPsf()->computeImage(center)),
            galsim::makeImageView(goodPixels),
            skyvar, _shearType.c_str(), 0xe,
            5.0, 3.0, // initial sigmas guesses; should feed previous shape estimate into these
            1E-6,     // precision used for convergence criteria
            center.getX(),
            center.getY()
        );

        afw::geom::ellipses::Quadrupole moments(
            afw::geom::ellipses::SeparableDistortionDeterminantRadius(
                data.observed_shape.distortion,
                data.moments_sigma
            )
        );
        
        source.set(_centroidKey, data.moments_centroid.getPoint());
        source.set(_momentsKey, moments);
#if 0
        source.set(_psfMomentsKey, shearEst.getPsfMoments());
#endif
        source.set(_resolutionKey, data.resolution_factor);
        
        if (data.meas_type == "e") {
            source.set(_errKey, 0.5 * data.corrected_shape_err);
            source.set(_e1Key, data.corrected_e1);
            source.set(_e2Key, data.corrected_e2);
        } else if (data.meas_type == "g") {
            source.set(_errKey, data.corrected_shape_err);
            source.set(_e1Key, data.corrected_g1);
            source.set(_e2Key, data.corrected_g2);
        }
        source.set(_flagKey, data.moments_status | data.correction_status); 
    } catch (galsim::hsm::HSMError & err) {
        source.set(_flagKey, 1);
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            err.what()
        );
    }
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

PTR(algorithms::AlgorithmControl) HsmShapeShapeletControl::_clone() const {
    return boost::make_shared<HsmShapeShapeletControl>(*this);
}

PTR(algorithms::Algorithm) HsmShapeBjControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmShape>(
        *this, "BJ", 'e', boost::ref(schema),
        "PSF-corrected shear using Bernstein & Jarvis (2002) method"
    );
}

PTR(algorithms::Algorithm) HsmShapeLinearControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmShape>(
        *this, "LINEAR", 'e', boost::ref(schema),
        "PSF-corrected shear using Hirata & Seljak (2003) 'linear' method"
    );
}

PTR(algorithms::Algorithm) HsmShapeKsbControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmShape>(
        *this, "KSB", 'g', boost::ref(schema),
        "PSF-corrected shear using KSB method"
    );
}

PTR(algorithms::Algorithm) HsmShapeRegaussControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<HsmShape>(
        *this, "REGAUSS", 'e', boost::ref(schema),
        "PSF-corrected shear using Hirata & Seljak (2003) 'regaussianization' method"
    );
}

PTR(algorithms::Algorithm) HsmShapeShapeletControl::_makeAlgorithm(
    afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
) const {
    std::string shearType = (boost::format("SHAPELET%d,%d") % maxOrderPsf % maxOrderGalaxy).str();
    return boost::make_shared<HsmShape>(
        *this, shearType, 'g', boost::ref(schema),
        "PSF-corrected shear using Hirata & Seljak (2003) 'shapelet' method"
    );
}

}}}} // namespace lsst::meas::extensions::shapeHSM
