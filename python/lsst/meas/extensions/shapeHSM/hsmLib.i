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

%{
#include "lsst/pex/logging.h"
#include "lsst/afw/image.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/algorithms.h"
%}

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()

%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/meas/algorithms/algorithmsLib.i"

%{
#include "lsst/meas/extensions/shapeHSM/HsmShapeControl.h"
#include "lsst/meas/extensions/shapeHSM/HsmMomentsControl.h"
%}

%shared_ptr(lsst::meas::extensions::shapeHSM::HsmShapeControl);
%shared_ptr(lsst::meas::extensions::shapeHSM::HsmShapeBjControl);
%shared_ptr(lsst::meas::extensions::shapeHSM::HsmShapeLinearControl);
%shared_ptr(lsst::meas::extensions::shapeHSM::HsmShapeKsbControl);
%shared_ptr(lsst::meas::extensions::shapeHSM::HsmShapeRegaussControl);
%shared_ptr(lsst::meas::extensions::shapeHSM::HsmSourceMomentsControl);
%shared_ptr(lsst::meas::extensions::shapeHSM::HsmPsfMomentsControl);

%include "lsst/meas/extensions/shapeHSM/HsmShapeControl.h"
%include "lsst/meas/extensions/shapeHSM/HsmMomentsControl.h"

