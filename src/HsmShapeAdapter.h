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
#if !defined(LSST_MEAS_EXTENSIONS_SHAPEHSM_ADAPTER_H)
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_ADAPTER_H

#include "psfcorr.h"

#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace shapeHSM {

/*
 * @class HsmShapeAdapter
 *
 * @brief Wrapper/adapter for Hirata/Seljac/Mandelbaum adaptive moment shape code.
 *
 */
    
template<typename ExposureT>
class HsmShapeAdapter {
public:

    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename ExposureT::MaskedImageT::Image ImageT;
    typedef typename ExposureT::MaskedImageT::Mask MaskT;
    typedef typename lsst::afw::detection::Psf::Image  PsfImageT;

    /*
     * constructor
     */
    HsmShapeAdapter(
                    CONST_PTR(ExposureT) exposure,
                    CONST_PTR(lsst::afw::detection::Peak) peak,
                    CONST_PTR(lsst::afw::detection::Source) source=PTR(lsst::afw::detection::Source)(),
                    lsst::afw::image::MaskPixel badPixelMask=0x0,
                    double background=0
                   );
    
    ~HsmShapeAdapter();

    /*
     * @brief Call the measurement routine for the specified algorithm
     * options include: BJ, KSB, REGAUSS, LINEAR, SHAPELET
     * @return status computed by HSM routine (0 = success, not 0 = failure)
     */
    unsigned measure(std::string shearType);

    /*
     * Accessors for the HSM outputs
     */
    double getX()     { return _galaxyData.x0 + _bbox.getMinX(); }
    double getY()     { return _galaxyData.y0 + _bbox.getMinY(); }
    double getE1()    { return _galaxyData.e1; }
    double getE1Err() { return getShearSig(); }
    double getE2()    { return _galaxyData.e2; }
    double getE2Err() { return getShearSig(); }
    
    double getFlux()  { return _galaxyData.flux; }
    double getResponsivity() { return _galaxyData.responsivity; }
    double getResolution()   { return _galaxyData.resolution; }
    char   getMeasType()     { return _galaxyData.meas_type; }

    double getIxx()   { return _galaxyData.mxx; }
    double getIyy()   { return _galaxyData.myy; }
    double getIxy()   { return _galaxyData.mxy; }
    double getSigma() { return _galaxyData.sigma; }

    double getPsfIxx()   { return _psfData.mxx; }
    double getPsfIyy()   { return _psfData.myy; }
    double getPsfIxy()   { return _psfData.mxy; }
    double getPsfSigma() { return _psfData.sigma; }

    double getShearSig() {
        return sqrt( 4.0*M_PI*_skyvar ) * _galaxyData.sigma / (_galaxyData.resolution * _galaxyData.flux);
    }
    
private:
    CONST_PTR(ExposureT) _exposure;
    CONST_PTR(lsst::afw::detection::Peak) _peak;
    CONST_PTR(lsst::afw::detection::Source) _source;
    lsst::afw::geom::Box2I _bbox;
    lsst::afw::image::MaskPixel _badPixelMask;
    RECT_IMAGE _atlasImage;
    RECT_IMAGE _psfImage;
    OBJECT_DATA _galaxyData;
    OBJECT_DATA _psfData;
    double _shearSig;
    double _skyvar;
};

}}}}

#endif
    
