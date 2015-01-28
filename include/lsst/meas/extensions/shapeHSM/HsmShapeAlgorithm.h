// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef LSST_MEAS_BASE_HsmShape_h_INCLUDED
#define LSST_MEAS_BASE_HsmShape_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/ShapeUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

/**
 *  @brief A C++ control class to handle HsmShapeAlgorithm's configuration
 *
 *  In C++, we define Control objects to handle configuration information.  Using the LSST_CONTROL_FIELD
 *  macro and lsst.pex.config.wrap.makeConfigClass, we can turn these into more full-featured Config classes
 *  in Python.  While the user will usually interact with the Config class, the plugin wrapper system will
 *  turn Config instances into Control instances when passing them to C++.
 *
 *  This should logically be an inner class, but Swig doesn't know how to parse those.
 */
class HsmShapeControl {
public:

    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>,
                       "Mask planes that indicate pixels that should be excluded from the fit");

    /**
     *  @brief Default constructor
     *
     *  All control classes should define a default constructor that sets all fields to their default values.
     */
    HsmShapeControl() : _name("") {
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("INTRP");
    }

    HsmShapeControl(std::string const & name) : _name(name) {
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("INTRP");
    }

private:
    std::string const & _name;
};

class HsmShapeBjControl : public HsmShapeControl {
public:
    HsmShapeBjControl() : HsmShapeControl("extensions_shapeHSM_HsmShapeBj") {}
};

class HsmShapeLinearControl : public HsmShapeControl {
public:
    HsmShapeLinearControl() : HsmShapeControl("extensions_shapeHSM_HsmShapeLinear") {}
};

class HsmShapeKsbControl : public HsmShapeControl {
public:
    HsmShapeKsbControl() : HsmShapeControl("extensions_shapeHSM_HsmShapeKsb") {}
};

class HsmShapeRegaussControl : public HsmShapeControl {
public:
    HsmShapeRegaussControl() : HsmShapeControl("extensions_shapeHSM_HsmShapeRegauss") {}
};

class HsmShapeShapeletControl : public HsmShapeControl {
public:
    LSST_CONTROL_FIELD(maxOrderPsf, int, "Maximum shapelet order for PSF");
    LSST_CONTROL_FIELD(maxOrderGalaxy, int, "Maximum shapelet order for galaxy");

    HsmShapeShapeletControl() : HsmShapeControl("extensions_shapeHSM_HsmShapeShapelet"), maxOrderPsf(8), maxOrderGalaxy(8) {}
};

/**
 *  @brief A measurement algorithm that estimates flux using a linear least-squares fit with the Psf model
 *
 *  The HsmShape algorithm is extremely simple: we do a least-squares fit of the Psf model (evaluated
 *  at a given position) to the data.  For point sources, this provides the optimal flux measurement
 *  in the limit where the Psf model is correct.  We do not use per-pixel weights in the fit, as this
 *  results in bright stars being fit with a different effective profile than faint stairs.
 */
class HsmShapeAlgorithm : public base::SimpleAlgorithm {
public:

    enum {
        FAILURE=base::FlagHandler::FAILURE,
        N_FLAGS
    };

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef HsmShapeControl Control;

    HsmShapeAlgorithm(Control const & ctrl, std::string const & name, std::string const & shearType, char measType,
                      std::string const & doc, afw::table::Schema & schema);

private:

    // These are private so they doesn't shadow the other overloads in base classes;
    // we can still call it via the public method on the base class.  We could have
    // used a using declaration instead, but Swig had trouble with that here.

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        base::MeasurementError * error=NULL
    ) const;

    Control _ctrl;
    std::string _shearType;
    char _measType;
    std::string const & _doc;

    base::CentroidResultKey _centroidResultKey;
    base::ShapeResultKey _momentsKey;
    base::ShapeResultKey _psfMomentsKey;
    afw::table::Key<double> _e1Key;
    afw::table::Key<double> _e2Key;
    afw::table::Key<double> _errKey;
    afw::table::Key<double> _sigmaKey;
    afw::table::Key<double> _resolutionKey;
    base::FlagHandler _flagHandler;
    base::SafeCentroidExtractor _centroidExtractor;
};

class HsmShapeBjAlgorithm : public HsmShapeAlgorithm {
    typedef HsmShapeBjControl Control;
    HsmShapeBjAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmShapeAlgorithm( ctrl, name + "Bj", "BJ", 'e', "PSF-corrected shear using Bernstein & Jarvis (2002) method", schema) {}
};

class HsmShapeLinearAlgorithm : public HsmShapeAlgorithm {
    typedef HsmShapeLinearControl Control;
    HsmShapeLinearAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmShapeAlgorithm(ctrl, name + "Linear", "LINEAR", 'e', "PSF-corrected shear using Hirata & Seljak (2003) 'linear' method", schema) {}
};

class HsmShapeKsbAlgorithm : public HsmShapeAlgorithm {
    typedef HsmShapeKsbControl Control;
    HsmShapeKsbAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmShapeAlgorithm(ctrl, name + "Ksb", "KSB", 'g', "PSF-corrected shear using KSB method", schema) {}
};

class HsmShapeRegaussAlgorithm : public HsmShapeAlgorithm {
    typedef HsmShapeRegaussControl Control;
    HsmShapeRegaussAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema) :
        HsmShapeAlgorithm(ctrl, name + "Regauss", "REGAUSS", 'e', "PSF-corrected shear using Hirata & Seljak (2003) 'regaussianization' method", schema) {}
};

class HsmShapeShapeletAlgorithm : public HsmShapeAlgorithm {
    typedef HsmShapeShapeletControl Control;
    HsmShapeShapeletAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);
};

}}}} // namespace lsst::meas::extensions::shapeHSM

#endif // !LSST_MEAS_BASE_HsmShape_h_INCLUDED
