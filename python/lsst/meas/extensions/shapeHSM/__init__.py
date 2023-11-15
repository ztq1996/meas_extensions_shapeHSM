#
# LSST Data Management System
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
# See the COPYRIGHT file
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#

"""lsst.meas.extensions.shapeHSM
"""
from lsst.meas.base import BasePlugin, wrapSimpleAlgorithm
from ._hsm_higher_moments import *
# from .hsmMomentsControl import *
# from .hsmShapeControl import *

from .version import *

# wrapSimpleAlgorithm(HsmShapeBjAlgorithm, name="ext_shapeHSM_HsmShapeBj",
#                     Control=HsmShapeBjControl, executionOrder=BasePlugin.SHAPE_ORDER)
# wrapSimpleAlgorithm(HsmShapeLinearAlgorithm, name="ext_shapeHSM_HsmShapeLinear",
#                     Control=HsmShapeLinearControl, executionOrder=BasePlugin.SHAPE_ORDER)
# wrapSimpleAlgorithm(HsmShapeKsbAlgorithm, name="ext_shapeHSM_HsmShapeKsb",
#                     Control=HsmShapeKsbControl, executionOrder=BasePlugin.SHAPE_ORDER)
# wrapSimpleAlgorithm(HsmShapeRegaussAlgorithm, name="ext_shapeHSM_HsmShapeRegauss",
#                     Control=HsmShapeRegaussControl, executionOrder=BasePlugin.SHAPE_ORDER)
# wrapSimpleAlgorithm(HsmSourceMomentsAlgorithm, name="ext_shapeHSM_HsmSourceMoments",
#                     Control=HsmSourceMomentsControl, executionOrder=BasePlugin.SHAPE_ORDER)
# wrapSimpleAlgorithm(HsmSourceMomentsRoundAlgorithm, name="ext_shapeHSM_HsmSourceMomentsRound",
#                     Control=HsmSourceMomentsRoundControl, executionOrder=BasePlugin.SHAPE_ORDER)
# wrapSimpleAlgorithm(HsmPsfMomentsAlgorithm, name="ext_shapeHSM_HsmPsfMoments",
#                     Control=HsmPsfMomentsControl, executionOrder=BasePlugin.SHAPE_ORDER)
# wrapSimpleAlgorithm(HsmPsfMomentsDebiasedAlgorithm, name="ext_shapeHSM_HsmPsfMomentsDebiased",
#                     Control=HsmPsfMomentsDebiasedControl, executionOrder=BasePlugin.FLUX_ORDER+1)
