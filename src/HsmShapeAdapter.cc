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
 * Measure adaptive moments using Hirata, Seljac, Mandelbaum (HSM) code
 *
 * Originally provided by Rachel Mandelbaum
 * This adapter/wrapper written by S. Bickerton
 */

#include "HsmShapeAdapter.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/PSF.h"

namespace afwMath        = lsst::afw::math;
namespace afwGeom        = lsst::afw::geom;
namespace afwImage       = lsst::afw::image;
namespace afwDetection   = lsst::afw::detection;
namespace extendShapeHsm = lsst::meas::extensions::shapeHSM;
namespace measAlg        = lsst::meas::algorithms;

/**
 * @brief Constructor for HsmShapeAdapter
 *
 */
template<typename ExposureT>
extendShapeHsm::HsmShapeAdapter<ExposureT>::HsmShapeAdapter(
    CONST_PTR(ExposureT) exposure,          ///< exposure containing the pixels to measure
    CONST_PTR(afwDetection::Peak) peak,     ///< location of the object to measure
    CONST_PTR(afwDetection::Source) source, ///< source object containing Footprint of the object
    afwImage::MaskPixel badPixelMask        ///< specify mask bits for pixels to *ignore*
                                                        ) :
    _exposure(exposure), _peak(peak), _source(source), _badPixelMask(badPixelMask) {


    // decide on the bbox, default to the whole image if the source or footprint aren't specified
    if (_source && _source->getFootprint()) {
        _bbox = _source->getFootprint()->getBBox();
    } else {
        CONST_PTR(ImageT) imTmp = exposure->getMaskedImage().getImage();
        _bbox = afwGeom::Box2I(afwGeom::Point2I(imTmp->getXY0()), 
                               afwGeom::Extent2I(imTmp->getWidth(), imTmp->getHeight()));
    }

    // make a shallow image copy in the bbox, allocate the image structure and copy into it
    CONST_PTR(MaskedImageT) img = 
        boost::make_shared<MaskedImageT>(exposure->getMaskedImage(), _bbox, afwImage::LOCAL, false);
    
    allocate_rect_image(&_atlasImage, 0, img->getWidth() - 1, 0, img->getHeight() - 1);
    for (int iY=0; iY < img->getHeight(); ++iY) {
        int iX = 0;
        for (typename ImageT::x_iterator ptr = img->getImage()->row_begin(iY); 
             ptr != img->getImage()->row_end(iY); ++ptr, ++iX) {
            _atlasImage.image[iX][iY] = static_cast<double>(*ptr);
        }
    }

    afwMath::StatisticsControl sctrl;
    sctrl.setAndMask(_badPixelMask);
    afwMath::Statistics stat = afwMath::makeStatistics(*img->getVariance(), *img->getMask(),
                                                       afwMath::MEDIAN, sctrl);
    _skyvar = sqrt(stat.getValue(afwMath::MEDIAN));
    
    // shallow copy the mask and use the user-provided badPixelMask to set the mask structure
    CONST_PTR(MaskT) msk =
        boost::make_shared<MaskT>(*exposure->getMaskedImage().getMask(), _bbox, afwImage::LOCAL, false);
    
    for (int iY = 0; iY < msk->getHeight(); ++iY) {
        int iX = 0;
        for (typename MaskT::x_iterator ptr = msk->row_begin(iY); ptr != msk->row_end(iY); ++ptr, ++iX) {
            _atlasImage.mask[iX][iY] = (_badPixelMask & *ptr) ? 0 : 1;
        }
    }
    
    // init the galaxyData
    double x = _peak->getFx() - _bbox.getMinX();
    double y = _peak->getFy() - _bbox.getMinY();
    _galaxyData.x0 = x;
    _galaxyData.y0 = y;

    
    // get a local image of the psf, allocate the psf structure and copy the psf into it
    typename PsfImageT::Ptr psf = exposure->getPsf()->computeImage(afwGeom::Point2D(x, y));
    allocate_rect_image(&_psfImage, 0, psf->getWidth() - 1, 0, psf->getHeight() - 1);
    for (int iY = 0; iY < psf->getHeight(); ++iY) {
        int iX = 0;
        for (typename PsfImageT::x_iterator ptr = psf->row_begin(iY); ptr != psf->row_end(iY); ++ptr, ++iX) {
            _psfImage.image[iX][iY] = static_cast<double>(*ptr);
        }
    }


    // init the psfData structure
    _psfData.x0 = 0.5 * (_psfImage.xmin + _psfImage.xmax);
    _psfData.y0 = 0.5 * (_psfImage.ymin + _psfImage.ymax);
    measAlg::PsfAttributes psfAttrib(exposure->getPsf(),
                                     static_cast<int>(_psfData.x0), static_cast<int>(_psfData.y0));
    _psfData.sigma = psfAttrib.computeGaussianWidth(measAlg::PsfAttributes::ADAPTIVE_MOMENT);
    _galaxyData.sigma = 2.5*_psfData.sigma;
    
}


/**
 * @brief Destructor to free memory allocated for the temp copies of the image and psf.
 *
 */
template<typename ExposureT>
extendShapeHsm::HsmShapeAdapter<ExposureT>::~HsmShapeAdapter() {
    deallocate_rect_image(&_atlasImage);    
    deallocate_rect_image(&_psfImage);
}


/**
 * @brief Forward the calculations through to the HSM code's general_shear_estimator
 *
 */
template<typename ExposureT>
unsigned extendShapeHsm::HsmShapeAdapter<ExposureT>::measure(
                                                             std::string shearType ///< algorithm to use
                                                            ) {
    try {
        char *cShearType = (char*) shearType.c_str();
        return  general_shear_estimator(&_atlasImage, &_psfImage,
                                        &_galaxyData, &_psfData,
                                        cShearType, 0xe);
    } catch (lsst::pex::exceptions::RuntimeErrorException& e) {
        LSST_EXCEPT_ADD(e, (boost::format("Using shape algorithm %s") % shearType).str());
        throw e;
    } 
}



#define INSTANTIATE(T) \
    template class extendShapeHsm::HsmShapeAdapter<afwImage::Exposure<T> >;


INSTANTIATE(int);
INSTANTIATE(float);
INSTANTIATE(double);
