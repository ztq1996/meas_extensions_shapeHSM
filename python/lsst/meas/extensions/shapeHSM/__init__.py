# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

from .hsmLib import *
from .version import *

import lsst.meas.base
lsst.meas.base.wrapSimpleAlgorithm(HsmShapeBjAlgorithm, name="ext_shapeHSM_HsmShapeBj",
    Control=HsmShapeBjControl, executionOrder=3.0)
lsst.meas.base.wrapSimpleAlgorithm(HsmShapeLinearAlgorithm, name="ext_shapeHSM_HsmShapeLinear",
    Control=HsmShapeLinearControl, executionOrder=3.0)
lsst.meas.base.wrapSimpleAlgorithm(HsmShapeKsbAlgorithm, name="ext_shapeHSM_HsmShapeKsb",
    Control=HsmShapeKsbControl, executionOrder=3.0)
lsst.meas.base.wrapSimpleAlgorithm(HsmShapeRegaussAlgorithm, name="ext_shapeHSM_HsmShapeRegauss",
    Control=HsmShapeRegaussControl, executionOrder=3.0)
lsst.meas.base.wrapSimpleAlgorithm(HsmSourceMomentsAlgorithm, name="ext_shapeHSM_HsmSourceMoments",
    Control=HsmSourceMomentsControl, executionOrder=3.0)
lsst.meas.base.wrapSimpleAlgorithm(HsmPsfMomentsAlgorithm, name="ext_shapeHSM_HsmPsfMoments",
    Control=HsmPsfMomentsControl, executionOrder=3.0)
del lsst # cleanup namespace
