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

#if !defined(LSST_MEAS_EXTENSIONS_SHAPEHSM_H)
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_H

/**
 * \file
 * Measure adaptive moments using the Hirata/Seljac/Mandelbaum code.
 */
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/detection/Shape.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace shapeHSM {
    
/**
 * @brief A class that knows how to calculate adaptive moment shape measurements
 */
class HsmShape : public lsst::afw::detection::Shape {
public:
    typedef boost::shared_ptr<HsmShape> Ptr;
    typedef boost::shared_ptr<HsmShape const> ConstPtr;
    
    /// constructor
    HsmShape(double x, double xErr, double y, double yErr,
             double ixx, double ixxErr, double ixy, double ixyErr, double iyy, double iyyErr) :
        lsst::afw::detection::Shape(x, xErr, y, yErr, ixx, ixxErr, ixy, ixyErr, iyy, iyyErr) {}
    
    /// Add desired fields to the schema
    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema) {
        Shape::defineSchema(schema);
    }

    template<typename ExposureT, typename AlgorithmT>
    static Shape::Ptr doMeasure(
                                CONST_PTR(ExposureT) im,
                                CONST_PTR(lsst::afw::detection::Peak) peak,
                                CONST_PTR(lsst::afw::detection::Source) source
                               );
    
    static bool doConfigure(lsst::pex::policy::Policy const& policy) {

        if (policy.exists("background")) {
            _background = policy.getDouble("background");
        }

        // a little sloppy ... these should only exist for HSM_SHAPELET
        // but it hurts no one to have them handled when defined by other HSM shapes.
        if (policy.exists("max_order_psf")) {
            _max_order_psf = policy.getInt("max_order_psf");
        }
        if (policy.exists("max_order_gal")) {
            _max_order_gal = policy.getInt("max_order_gal");
        }

        if (policy.exists("badmaskplanes")) {
            std::vector<std::string> planes = policy.getStringArray("badmaskplanes");
            _badMaskPlanes.clear();
            for (std::vector<std::string>::const_iterator it = planes.begin(); it != planes.end(); ++it) {
                _badMaskPlanes.push_back(*it);
            }
        }
        
        return true;
    }
    
private:
    static double _background;
    static int _max_order_psf;
    static int _max_order_gal;
    static std::vector<std::string> _badMaskPlanes;
    static short _status;
    HsmShape(void) : lsst::afw::detection::Shape() { }
    LSST_SERIALIZE_PARENT(lsst::afw::detection::Shape);
};

}}}}

#endif
