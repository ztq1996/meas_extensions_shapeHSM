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

from lsst.meas.base import *
wrapSimpleAlgorithm(HsmShapeBjAlgorithm, needsMetadata=False, Control=HsmShapeBjControl, executionOrder=3.0)
wrapSimpleAlgorithm(HsmShapeLinearAlgorithm, needsMetadata=False, Control=HsmShapeLinearControl, executionOrder=3.0)
wrapSimpleAlgorithm(HsmShapeKsbAlgorithm, needsMetadata=False, Control=HsmShapeKsbControl, executionOrder=3.0)
wrapSimpleAlgorithm(HsmShapeRegaussAlgorithm, needsMetadata=False, Control=HsmShapeRegaussControl, executionOrder=3.0)
wrapSimpleAlgorithm(HsmShapeShapeletAlgorithm, needsMetadata=False, Control=HsmShapeShapeletControl, executionOrder=3.0)
del lsst # cleanup namespace
