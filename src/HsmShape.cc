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

#include "lsst/meas/extensions/shapeHSM/HsmShapeAdapter.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace measAlg  = lsst::meas::algorithms;

namespace extendShapeHsm = lsst::meas::extensions::shapeHSM;

namespace lsst { namespace meas { namespace algorithms {

template<typename ExposureT>
class HsmShapeBase : public Algorithm<afwDet::Shape, ExposureT>
{
public:
    typedef Algorithm<afwDet::Shape, ExposureT> AlgorithmT;
    typedef boost::shared_ptr<HsmShapeBase> Ptr;
    typedef boost::shared_ptr<HsmShapeBase const> ConstPtr;

    /// Ctor
    HsmShapeBase(std::vector<std::string> const& badMaskPlanes=std::vector<std::string>()) :
        AlgorithmT(), _badMaskPlanes(badMaskPlanes) {}

    virtual void configure(pexPolicy::Policy const& policy) {
        if (policy.exists("badmaskplanes")) {
            std::vector<std::string> planes = policy.getStringArray("badmaskplanes");
            _badMaskPlanes.clear();
            for (std::vector<std::string>::const_iterator it = planes.begin(); it != planes.end(); ++it) {
                _badMaskPlanes.push_back(*it);
            }
        }
    }

    virtual PTR(afwDet::Shape) measureNull(void) const {
        double const NaN = std::numeric_limits<double>::quiet_NaN();
        return boost::shared_ptr<afwDet::Shape>(new afwDet::Shape(NaN, NaN, NaN, NaN, NaN, 
                                                                  NaN, NaN, NaN, NaN, NaN));
    }

    // We usually just want alg.getName() to pass to the
    // ShearEstimator.measure() method, but HSM_SHAPELET also needs two integers
    // appended to the string, so we must handle that.
    virtual std::string getName() const = 0;
    virtual std::string getAdapterName() const { return getName().substr(4); }

    virtual PTR(afwDet::Shape) measureSingle(afwDet::Source const&, afwDet::Source const&,
                                             ExposurePatch<ExposureT> const&) const;
protected:
    std::vector<std::string> _badMaskPlanes;
};

/// CRTP to clone the derived classes that are all the same
template<typename ExposureT, typename HsmShapeDerived>
class HsmShapeCloner : public HsmShapeBase<ExposureT> {
public:
    HsmShapeCloner(std::vector<std::string> const& badMaskPlanes) : HsmShapeBase<ExposureT>(badMaskPlanes) {}
    virtual PTR(typename HsmShapeBase<ExposureT>::AlgorithmT) clone() const {
        return boost::make_shared<HsmShapeDerived>(this->_badMaskPlanes);
    }
};

template<typename ExposureT> 
class HsmShapeBj : public HsmShapeCloner<ExposureT, HsmShapeBj<ExposureT> > {
public:
    HsmShapeBj(std::vector<std::string> const& badMaskPlanes=std::vector<std::string>()) : 
        HsmShapeCloner<ExposureT, HsmShapeBj<ExposureT> >(badMaskPlanes) {}
    virtual std::string getName() const { return "HSM_BJ"; }
};

template<typename ExposureT>
class HsmShapeLinear : public HsmShapeCloner<ExposureT, HsmShapeLinear<ExposureT> > {
public:
    HsmShapeLinear(std::vector<std::string> const& badMaskPlanes=std::vector<std::string>()) : 
        HsmShapeCloner<ExposureT, HsmShapeLinear<ExposureT> >(badMaskPlanes) {}
    virtual std::string getName() const { return "HSM_LINEAR"; }
};

template<typename ExposureT>
class HsmShapeKsb : public HsmShapeCloner<ExposureT, HsmShapeKsb<ExposureT> > {
public:
    HsmShapeKsb(std::vector<std::string> const& badMaskPlanes=std::vector<std::string>()) : 
        HsmShapeCloner<ExposureT, HsmShapeKsb<ExposureT> >(badMaskPlanes) {}
    virtual std::string getName() const { return "HSM_KSB"; }
};

template<typename ExposureT>
class HsmShapeRegauss : public HsmShapeCloner<ExposureT, HsmShapeRegauss<ExposureT> > {
public:
    HsmShapeRegauss(std::vector<std::string> const& badMaskPlanes=std::vector<std::string>()) : 
        HsmShapeCloner<ExposureT, HsmShapeRegauss<ExposureT> >(badMaskPlanes) {}
    virtual std::string getName() const { return "HSM_REGAUSS"; }
};

template<typename ExposureT>
class HsmShapeShapelet : public HsmShapeBase<ExposureT> {
public:
    HsmShapeShapelet(int max_order_psf=8, int max_order_gal=8, 
                     std::vector<std::string> badMaskPlanes=std::vector<std::string>()) :
        HsmShapeBase<ExposureT>(badMaskPlanes), _max_order_psf(max_order_psf), 
        _max_order_gal(max_order_gal) {}

    virtual std::string getName() const { return "HSM_SHAPELET"; }
    virtual std::string getAdapterName() const {
        return getName().substr(4) + (boost::format("%d,%d") % _max_order_psf % _max_order_gal).str();
    }

    virtual void configure(pexPolicy::Policy const& policy) {
        if (policy.exists("max_order_psf")) {
            _max_order_psf = policy.getInt("max_order_psf");
        }
        if (policy.exists("max_order_gal")) {
            _max_order_gal = policy.getInt("max_order_gal");
        }
        HsmShapeBase<ExposureT>::configure(policy);
    }

    virtual PTR(typename HsmShapeBase<ExposureT>::AlgorithmT) clone() const {
        return boost::make_shared<HsmShapeShapelet>(_max_order_psf, _max_order_gal, 
                                                    this->_badMaskPlanes);
    }

private:
    int _max_order_psf;
    int _max_order_gal;
};


/**
 * @brief Given an image and a pixel position, return a Shape using the HSM algorithm
 */
template<typename ExposureT>
PTR(afwDet::Shape) HsmShapeBase<ExposureT>::measureSingle(
    afwDet::Source const& target,
    afwDet::Source const& source,
    ExposurePatch<ExposureT> const& patch
    ) const
{
    CONST_PTR(ExposureT) exposure = patch.getExposure();
    CONST_PTR(afwDet::Footprint) foot = patch.getFootprint();

    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename ExposureT::MaskedImageT::Image ImageT;
    
    afwImage::MaskPixel badPixelMask =
        exposure->getMaskedImage().getMask()->getPlaneBitMask(_badMaskPlanes);
    extendShapeHsm::HsmShapeAdapter<ExposureT> shearEst(exposure, patch.getCenter(), *foot, badPixelMask);
    
    short status = shearEst.measure(getAdapterName());   
 
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
    PTR(afwDet::Shape) shape = PTR(afwDet::Shape)(new afwDet::Shape(x, xErr, y, yErr, ixx, ixxErr, 
                                                                    ixy, ixyErr, iyy, iyyErr));

    shape->set<afwDet::Shape::SHAPE_STATUS, short>(status);
    shape->set<afwDet::Shape::SIGMA, double>(shearEst.getSigma());
    
    // if we have ellipticities 'e' put them in e1/e2
    // if we have shear 'g', put them in shear1/shear2
    char meas_type = shearEst.getMeasType();
    if (meas_type == 'e') {
        shape->set<afwDet::Shape::E1, double>(shearEst.getE1());
        shape->set<afwDet::Shape::E1_ERR, double>(0.5 * shearEst.getShearSig());
        shape->set<afwDet::Shape::E2, double>(shearEst.getE2());
        shape->set<afwDet::Shape::E2_ERR, double>(0.5 * shearEst.getShearSig());
    } else if (meas_type == 'g') {
        shape->set<afwDet::Shape::SHEAR1, double>(shearEst.getE1());
        shape->set<afwDet::Shape::SHEAR1_ERR, double>(shearEst.getShearSig());
        shape->set<afwDet::Shape::SHEAR2, double>(shearEst.getE2());
        shape->set<afwDet::Shape::SHEAR2_ERR, double>(shearEst.getShearSig());
    }

    shape->set<afwDet::Shape::RESOLUTION, double>(shearEst.getResolution());
    
    shape->set<afwDet::Shape::PSF_IXX, double>(shearEst.getPsfIxx());
    shape->set<afwDet::Shape::PSF_IXY, double>(shearEst.getPsfIxy());
    shape->set<afwDet::Shape::PSF_IYY, double>(shearEst.getPsfIyy());

    return shape;
}

LSST_DECLARE_ALGORITHM(HsmShapeBj, afwDet::Shape);
LSST_DECLARE_ALGORITHM(HsmShapeLinear, afwDet::Shape);
LSST_DECLARE_ALGORITHM(HsmShapeKsb, afwDet::Shape);
LSST_DECLARE_ALGORITHM(HsmShapeRegauss, afwDet::Shape);
LSST_DECLARE_ALGORITHM(HsmShapeShapelet, afwDet::Shape);


}}} // namespace
