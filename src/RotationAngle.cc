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
 
// -*- LSST-C++ -*-
#include <cmath>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"


#include <stdio.h>

namespace afwDetection = lsst::afw::detection;
namespace afwCoord = lsst::afw::coord;
namespace afwImage = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/**
 * @brief A class that calculates the rotation angles (angles between x axis of CCD and East,North)
 */
class RotationAngle : public afwDetection::Astrometry
{
protected:
    enum { EAST = NVALUE, NORTH };
public:
    typedef boost::shared_ptr<RotationAngle> Ptr;
    typedef boost::shared_ptr<RotationAngle const> ConstPtr;

    /// Ctor
    RotationAngle(double east, double north) {
        init();
        // Everything has to be set, even to a meaningless value
        set<X>(std::numeric_limits<double>::quiet_NaN());
        set<X_ERR>(std::numeric_limits<double>::quiet_NaN());
        set<Y>(std::numeric_limits<double>::quiet_NaN());
        set<Y_ERR>(std::numeric_limits<double>::quiet_NaN());

        set<EAST>(east);
        set<NORTH>(north);
    }
    
    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema) {
        schema->add(afwDetection::SchemaEntry("east", EAST, afwDetection::Schema::DOUBLE, 1, "radians"));
        schema->add(afwDetection::SchemaEntry("north", NORTH, afwDetection::Schema::DOUBLE, 1, "radians"));
    }

    template<typename ExposureT>
    static Astrometry::Ptr doMeasure(CONST_PTR(ExposureT),
                                     CONST_PTR(afwDetection::Peak),
                                     CONST_PTR(afwDetection::Source)
                                    );

    double getEast() const { return afwDetection::Astrometry::get<EAST, double>(); }
    double getNorth() const { return afwDetection::Astrometry::get<NORTH, double>(); }
};

template<typename ExposureT>
afwDetection::Astrometry::Ptr RotationAngle::doMeasure(CONST_PTR(ExposureT) image,
                                                       CONST_PTR(afwDetection::Peak) peak,
                                                       CONST_PTR(afwDetection::Source) source
    )
{
    double east = std::numeric_limits<double>::quiet_NaN();
    double north = std::numeric_limits<double>::quiet_NaN();
    if (peak) {
        afwImage::Wcs::ConstPtr wcs = image->getWcs();
        if (wcs) {
            afwGeom::Point2D const pix(peak->getFx(), peak->getFy());
            afwGeom::AffineTransform const lin = wcs->linearizePixelToSky(pix, afwCoord::RADIANS);
            afwGeom::AffineTransform::ParameterVector const param = lin.getParameterVector();
            east = ::atan2(param[lin.XY], param[lin.XX]);
            north = ::atan2(param[lin.YY], param[lin.YX]);
        }
    }   
    return boost::make_shared<RotationAngle>(east, north);
}

/*
 * Declare the existence of the algorithm to MeasureAstrometry
 *
 * @cond
 */
#define INSTANTIATE(TYPE) \
    MeasureAstrometry<afwImage::Exposure<TYPE> >::declare("ROTANGLE", \
        &RotationAngle::doMeasure<afwImage::Exposure<TYPE> > \
        )

volatile bool isInstance[] = {
    INSTANTIATE(int),
    INSTANTIATE(float),
    INSTANTIATE(double),
};

// \endcond

}}}}
