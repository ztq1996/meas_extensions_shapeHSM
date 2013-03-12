// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_Bounds_h_INCLUDED
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_Bounds_h_INCLUDED

/*
 * This is an adapter file that provides a minimal reimplementation
 * of the GalSim Bounds object using LSST objects, defined here
 * so that we can use PSFCorr.[h,cc] with minimal modification.
 */

#include "lsst/afw/geom/Box.h"

#define dbg if (false) std::cerr

namespace galsim {

template <typename T> struct Bounds;

template <typename T>
struct Position {

    explicit Position(T x_, T y_) : x(x_), y(y_) {}

    lsst::afw::geom::Point<T,2> getPoint() const { return lsst::afw::geom::Point<T,2>(x, y); }

    T x;
    T y;
};

template <>
struct Bounds<int> {

    explicit Bounds() : box() {}

    explicit Bounds(int xmin, int xmax, int ymin, int ymax) :
        box(lsst::afw::geom::Point2I(xmin, ymin), lsst::afw::geom::Point2I(xmax, ymax))
    {}

    explicit Bounds(lsst::afw::geom::Box2I const & box_) : box(box_) {}

    int getXMin() const { return box.getMinX(); }
    int getYMin() const { return box.getMinY(); }
    int getXMax() const { return box.getMaxX(); }
    int getYMax() const { return box.getMaxX(); }

    lsst::afw::geom::Box2I box;
};

} // namespace galsim

#endif // !LSST_MEAS_EXTENSIONS_SHAPEHSM_Bounds_h_INCLUDED
