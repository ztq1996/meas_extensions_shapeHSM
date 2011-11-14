// -*- lsst-c++ -*-

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

/* **EDIT** fooLib (1 place) */
%define hsmLib_DOCSTRING
"
Various swigged-up C++ classes for testing
"
%enddef

/* **EDIT** fooLib (3 places) */
%feature("autodoc", "1");
%module(package="lsst.meas.extensions.shapeHSM.hsmLib", docstring=hsmLib_DOCSTRING) hsmLib

%pythonnondynamic;
%naturalvar;  // use const reference typemaps

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()

%import "lsst/afw/detection/detectionLib.i"

/* **EDIT** all remaining lines to include your header and handle shared pointer to your class */
%{
#include "lsst/meas/extensions/shapeHSM/HsmShapeAdapter.h"
%}

SWIG_SHARED_PTR_DERIVED(HsmShapePtr, lsst::afw::detection::Shape, lsst::meas::extensions::shapeHSM::HsmShape);

%include "lsst/meas/extensions/shapeHSM/HsmShapeAdapter.h"


%define %declareShape(PIXTYPE, SUFFIX)
%template(HsmShapeAdapter ## SUFFIX) lsst::meas::extensions::shapeHSM::HsmShapeAdapter<lsst::afw::image::Exposure<PIXTYPE> >;
%enddef

%declareShape(double, D)
%declareShape(float, F)
%declareShape(int, I)
