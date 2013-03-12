// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_Image_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_Image_h_INCLUDED

/*
 * This is an adapter file that provides a minimal reimplementation
 * of the GalSim Image objects using LSST objects, defined here
 * so that we can use PSFCorr.[h,cc] with minimal modification.
 */

#include "boost/make_shared.hpp"

#include "lsst/afw/image/Image.h"

#include "Bounds.h"

namespace galsim {

template <typename T>
struct ConstImageView {

    explicit ConstImageView(lsst::afw::geom::Point2I const & xy0_, ndarray::Array<T,2,1> const & array_) :
        xy0(xy0_), array(array_)
    {}

    T const * getData() const { return array.getData(); }
    
    int getStride() const { return array.template getStride<0>(); }

    T const & operator()(int x, int y) const { return array[ndarray::makeVector(y, x)]; }
    
    int getXMin() const { return xy0.getX(); }
    int getYMin() const { return xy0.getY(); }
    int getXMax() const { return xy0.getX() + array.template getSize<1>() - 1; }
    int getYMax() const { return xy0.getY() + array.template getSize<0>() - 1; }
    
    Bounds<int> getBounds() const { return Bounds<int>(getXMin(), getXMax(), getYMin(), getYMax()); }

    ConstImageView<T> view() const { return ConstImageView<T>(*this); }
    
    lsst::afw::geom::Point2I xy0;
    ndarray::Array<T,2,1> array;
};

template <typename T>
struct ImageView : public ConstImageView<T> {

    explicit ImageView(lsst::afw::geom::Point2I const & xy0_, ndarray::Array<T,2,1> const & array_) :
        ConstImageView<T>(xy0_, array_)
    {}

    T * getData() const { return this->array.getData(); }

    int getStride() const { return this->array.template getStride<0>(); }

    T & operator()(int x, int y) const { return this->array[ndarray::makeVector(y, x)]; }
    
    ImageView<T> const & operator*=(T v) const {
        this->array.deep() *= v;
        return *this;
    }

    void fill(T v) { this->array.deep() = v; }

    ImageView<T> view() const { return ImageView<T>(*this); }
};

template <typename T>
struct Image : public ImageView<T> {

    explicit Image(Bounds<int> const & bounds) :
        ImageView<T>(bounds.box.getMin(), ndarray::allocate(bounds.box.getWidth(), bounds.box.getHeight()))
    {}

    Image(Image<T> const & other) : ImageView<T>(other.xy0, ndarray::copy(other.array)) {}

    Image(ConstImageView<T> const & other) : ImageView<T>(other.xy0, ndarray::copy(other.array)) {}

    template <typename U>
    Image(ConstImageView<U> const & other) :
        ImageView<T>(other.xy0, ndarray::allocate(other.array.getShape()))
    {
        this->array.deep() = other.array;
    }

};

template <typename T>
inline ImageView<T> makeImageView(lsst::afw::image::Image<T> image) {
    return ImageView<T>(image.getXY0(), image.getArray());
}

} // namespace galsim

#endif // !LSST_MEAS_EXTENSIONS_SHAPEHSM_Image_h_INCLUDED
