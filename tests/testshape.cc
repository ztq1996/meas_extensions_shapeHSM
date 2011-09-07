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
 
//
// test a perfect Gaussian PSF and measure aperture photometry at different radii
//
#include <iostream>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/ImageAlgorithm.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/afw/coord/Coord.h"

#include "lsst/meas/extensions/shapeHSM/HsmShape.h"


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HsmShape

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

using namespace std;
namespace pexPolicy = lsst::pex::policy;
namespace measAlgorithms = lsst::meas::algorithms;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace afwGeom = lsst::afw::geom;
namespace afwCoord = lsst::afw::coord;

typedef afwImage::Exposure<float, short unsigned int, float> ExposureT;
typedef ExposureT::MaskedImageT MImage;


/**
 */

BOOST_AUTO_TEST_CASE(HsmMoments) {

    float bkgd = 100.0;
    int xwidth = 100;
    int ywidth = 100;
    MImage mimg(afwGeom::Extent2I(xwidth, ywidth));
    *mimg.getImage() = bkgd;

    int x = 30;
    int y = 40;
    float I0 = 1000.0;

    // make a masked image and an exposure
    (*mimg.getImage())(x,y) = I0 + bkgd;
    *mimg.getMask() = 0;
    *mimg.getVariance() = sqrt(I0);
    ExposureT::Ptr exposure(new ExposureT(mimg));

    int kwid = 35;
    float sigma = 1.0;
    afwDetection::Psf::Ptr psf = afwDetection::createPsf("SingleGaussian", kwid, kwid, sigma);
    exposure->setPsf(psf);
    exposure->setWcs(*afwImage::makeWcs(afwCoord::makeCoord(afwCoord::ICRS, 0.0 * afwGeom::degrees, 0.0 * afwGeom::degrees),
                                        afwGeom::Point2D(1.0, 1.0),
                                        0.2/3600.0, 0.0, 0.0, 0.2/3600.0));
    
    // Add a Gaussian to the image
    float sigma_xx = std::pow(1.5, 2);
    float sigma_xy = 0;
    float sigma_yy = std::pow(2.5, 2);
    int ksize = 15;


    // compute the known values of the moments
    float phi = 0.0;  // rotate +ve this far (in radians)
    float c = cos(phi);
    float s = sin(phi);

    float sum, sumxx, sumxy, sumyy;
    sum = sumxx = sumxy = sumyy = 0.0;

    for (int dx = -(ksize/2); dx < ksize/2 + 1; dx++) {
        for (int dy = -(ksize/2); dy < ksize/2 + 1; dy++) {
            float u = c*dx + s*dy;
            float v = s*dx + c*dy;
            float I = I0*exp(-0.5*(u*u/sigma_xx + v*v/sigma_yy));
            MImage::iterator ptr = mimg.at(x+dx, y+dy);
            ptr.image() = bkgd + I;
            sum += I;
            sumxx += I*dx*dx;
            sumxy += I*dx*dy;
            sumyy += I*dy*dy;
        }
    }
    sumxx /= sum;
    sumxy /= sum;
    sumyy /= sum;
        
    printf("%g %g %g\n", sumxx, sumxy, sumyy);


    // measure the shape with our algorithms
    std::vector<std::string> algList;
    algList.push_back("HSM_BJ");
    algList.push_back("HSM_LINEAR");
    algList.push_back("HSM_KSB");
    algList.push_back("HSM_REGAUSS");
    algList.push_back("HSM_SHAPELET");

    PTR(MImage::Image) img = exposure->getMaskedImage().getImage();
    *img -= bkgd;

    for (std::vector<std::string>::iterator it = algList.begin(); it != algList.end(); ++it) {

        std::string alg = *it;
        
        measAlgorithms::MeasureShape<ExposureT>::Ptr measureShape =
            measAlgorithms::makeMeasureShape<ExposureT::Ptr>();
        measureShape->addAlgorithm(alg);
        pexPolicy::Policy policy;
        policy.set(alg+".background", bkgd);
        policy.set(alg+".max_order_psf", 8);
        policy.set(alg+".max_order_gal", 8);
        policy.set(alg+".badmaskplanes", "BAD");
        policy.add(alg+".badmaskplanes", "SAT");
        
        measureShape->configure(policy);
        CONST_PTR(afwDetection::Peak) peak = boost::make_shared<afwDetection::Peak>(x, y);
        measureShape->setImage(exposure);

        afwDetection::Shape::Ptr shape = measureShape->measure(peak)->find(alg);

        // compare to known values
        float Ixx = shape->getIxx();
        float Iyy = shape->getIyy();
        float Ixy = shape->getIxy();
        float A2 = 0.5*(Ixx + Iyy) + sqrt( pow(0.5*(Ixx - Iyy), 2) + pow(Ixy, 2) );
        float B2 = 0.5*(Ixx + Iyy) - sqrt( pow(0.5*(Ixx - Iyy), 2) + pow(Ixy, 2) );

        printf("Algorithm: %s\n", alg.c_str());
        printf("x:     %d %.2f\n", x, shape->getX());
        printf("y:     %d %.2f\n", y, shape->getY());
        printf("I_xx:  %.5f %.5f\n", sigma_xx, Ixx);
        printf("I_xy:  %.5f %.5f\n", sigma_xy, Ixy);
        printf("I_yy:  %.5f %.5f\n", sigma_yy, Iyy);
        printf("A2, B2 = %.5f, %.5f\n", A2, B2);      
    
        // 1/100 pixel limit for x,y
        BOOST_CHECK(fabs(static_cast<double>(x) - shape->getX()) < 1e-2);
        BOOST_CHECK(fabs(static_cast<double>(y) - shape->getY()) < 1e-2);

        // higher limit for ixx/xy/yy
        BOOST_CHECK(fabs(Ixx - sigma_xx) < 0.02*(1.0+sigma_xx));
        BOOST_CHECK(fabs(Ixy - sigma_xy) < 0.02*(1.0+sigma_xy));
        BOOST_CHECK(fabs(Iyy - sigma_yy) < 0.03*(1.0+sigma_yy));
        
    }
}
    
