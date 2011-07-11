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

/**
 * \file
 * Measure adaptive moments.
 *
 * Originally provided by Phil Fischer, based on code from Tim McKay's group.  Error calculations by Dave
 * Johnston.  Major reworking by RHL for SDSS, and now a major rewrite for LSST
 */

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/meas/extensions/shapeHSM/HsmShape.h"
#include "lsst/meas/extensions/shapeHSM/HsmShapeAdapter.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace measAlg  = lsst::meas::algorithms;

namespace extendShapeHsm = lsst::meas::extensions::shapeHSM;

LSST_REGISTER_SERIALIZER(extendShapeHsm::HsmShape)

int extendShapeHsm::HsmShape::_max_order_psf = 8;   // used only for HSM_SHAPELET
int extendShapeHsm::HsmShape::_max_order_gal = 8;   // used only for HSM_SHAPELET
std::vector<std::string> extendShapeHsm::HsmShape::_badMaskPlanes = std::vector<std::string>();
short extendShapeHsm::HsmShape::_status = -1; 

// Use templates to allow the different algorithms to registered with the same boilerplate code.
// The adapter just needs to be called with a different string for the corresponding algorithm
// We don't want to rewrite the same code for each one.  By registering (see below) the names
// for templates instantiated with these classes, we can pull out the string (char* rather).
namespace {
    class Bj {
    public:
        const char *getName() { return std::string("BJ").c_str(); }
    };
    class Linear {
    public:
        const char *getName() { return std::string("LINEAR").c_str(); }
    };
    class Ksb {
    public:
        const char *getName() { return std::string("KSB").c_str(); }
    };
    class Regauss {
    public:
        const char *getName() { return std::string("REGAUSS").c_str(); }
    };
    class Shapelet {
    public:
        const char *getName() { return std::string("SHAPELET").c_str(); }
    };
}


/**
 * @brief Given an image and a pixel position, return a Shape using the HSM algorithm
 */
template<typename ExposureT, typename AlgorithmT>
afwDetection::Shape::Ptr extendShapeHsm::HsmShape::doMeasure(
                                                          CONST_PTR(ExposureT) exposure,
                                                          CONST_PTR(afwDetection::Peak) peak,
                                                          CONST_PTR(afwDetection::Source) source
                                                         ) {
    if (!peak) {
        double const NaN = std::numeric_limits<double>::quiet_NaN();
        return boost::shared_ptr<extendShapeHsm::HsmShape>(new extendShapeHsm::HsmShape(NaN, NaN, NaN, NaN,
                                                                                        NaN, NaN,
                                                                                        NaN, NaN, NaN, NaN));
    }

    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename ExposureT::MaskedImageT::Image ImageT;
    
    afwImage::MaskPixel badPixelMask =
        exposure->getMaskedImage().getMask()->getPlaneBitMask(_badMaskPlanes);
    typename extendShapeHsm::HsmShapeAdapter<ExposureT>
        shearEst(exposure, peak, source, badPixelMask);

    // some fussiness here.
    // we usually just want alg.getName() to pass to the ShearEstimator.measure() method.
    // but HSM_SHAPELET also needs two integers appended to the char*, so we must handle that.
    AlgorithmT alg;
    int const maxChar = 32;
    char algName[maxChar];
    strncpy(algName, alg.getName(), maxChar);
    
    if (strncmp(algName, "SHAPELET", maxChar) == 0) {
        
        // limit order sizes
        // partly to prevent exceeding maxChar characters, but mainly to avoid absurdly high values
        // for the algorithm
        if ((_max_order_psf < 100) && (_max_order_gal < 100)) {
            sprintf(algName, "%s%d,%d", alg.getName(), _max_order_psf, _max_order_gal);
        } else {
            throw LSST_EXCEPT(pexExceptions::InvalidParameterException,
                              (boost::format("max_order_psf (%d) and max_order_gal (%d) must be <100") %
                               _max_order_psf % _max_order_gal).str());
        }
    }
    _status = shearEst.measure(algName);   
 
    float x = shearEst.getX();
    float xErr = 0.0;
    float y = shearEst.getY();
    float yErr = 0.0;
    
    float ixx = shearEst.getIxx();
    float ixy = shearEst.getIxy();
    float iyy = shearEst.getIyy();
    float ixxErr = 0;
    float ixyErr = 0;
    float iyyErr = 0;

    // flags: SHAPE_MAXITER, SHAPE_UNWEIGHTED, SHAPE_UNWEIGHTED_BAD
    //shape->setFlags(shape->getFlags() | Flags::SHAPE_UNWEIGHTED);
    extendShapeHsm::HsmShape::Ptr shape =
        extendShapeHsm::HsmShape::Ptr(new extendShapeHsm::HsmShape(x, xErr, y, yErr,
                                                                   ixx, ixxErr,
                                                                   ixy, ixyErr,
                                                                   iyy, iyyErr));

    shape->set<SHAPE_STATUS>(_status);
    shape->set<SIGMA>(shearEst.getSigma());
    
    // if we have ellipticities 'e' put them in e1/e2
    // if we have shear 'g', put them in shear1/shear2
    char meas_type = shearEst.getMeasType();
    if (meas_type == 'e') {
        shape->set<E1>(shearEst.getE1());
        shape->set<E1_ERR>(0.5 * shearEst.getShearSig());
        shape->set<E2>(shearEst.getE2());
        shape->set<E2_ERR>(0.5 * shearEst.getShearSig());
    } else if (meas_type == 'g') {
        shape->set<SHEAR1>(shearEst.getE1());
        shape->set<SHEAR1_ERR>(shearEst.getShearSig());
        shape->set<SHEAR2>(shearEst.getE2());
        shape->set<SHEAR2_ERR>(shearEst.getShearSig());
    }

    shape->set<RESOLUTION>(shearEst.getResolution());
    
    shape->set<PSF_IXX>(shearEst.getPsfIxx());
    shape->set<PSF_IXY>(shearEst.getPsfIxy());
    shape->set<PSF_IYY>(shearEst.getPsfIyy());

    return shape;
}

/*
 * Declare the existence of the algorithm to MeasureShape
 */
#define ALG_INSTANTIATE(TYPE, LABEL, ALG)                                       \
    measAlg::MeasureShape<afwImage::Exposure<TYPE> >::declare(LABEL,  \
        &extendShapeHsm::HsmShape::doMeasure<afwImage::Exposure<TYPE>, ALG>, \
        &extendShapeHsm::HsmShape::doConfigure \
        )

#define INSTANTIATE(TYPE)                       \
    ALG_INSTANTIATE(TYPE, "HSM_BJ", Bj),        \
        ALG_INSTANTIATE(TYPE, "HSM_LINEAR", Linear),     \
        ALG_INSTANTIATE(TYPE, "HSM_KSB", Ksb),            \
        ALG_INSTANTIATE(TYPE, "HSM_REGAUSS", Regauss), \
        ALG_INSTANTIATE(TYPE, "HSM_SHAPELET", Shapelet)

volatile bool isInstance[] = {
    INSTANTIATE(int),                             
    INSTANTIATE(float),           
    INSTANTIATE(double)
};                                          



