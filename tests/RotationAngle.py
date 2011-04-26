#!/usr/bin/env python

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

import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests

import lsst.meas.extensions.shapeHSM

try:
    type(verbose)
except NameError:
    verbose = 0

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class RotationAngleTestCase(unittest.TestCase):
    """A test case for rotation angle"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testMeasure(self):
        """Test that we can instantiate and play with RotationAngle"""

        scale = 5.0e-5                  # Pixel scale: 0.18 arcsec/pixel
        for angle in (0, 45, 90, 120, 180, 275):
            for expFactory in (afwImage.ExposureF, afwImage.ExposureI):
                angle = math.radians(angle)
                cdMatrix = [[scale * math.cos(angle), scale * math.sin(angle)],
                            [scale * -math.sin(angle), scale * math.cos(angle)]]
                
                wcs = afwImage.Wcs(afwGeom.Point2D(0, 0), afwGeom.Point2D(50, 50), cdMatrix)
                exp = expFactory(afwGeom.ExtentI(100, 100), wcs)

                astrom = algorithms.makeMeasureAstrometry(exp)
                astrom.addAlgorithm("ROTANGLE")

                x, y = 10, 20
                ast = astrom.measure(afwDetection.Peak(x, y)).find()

                east, north = ast.get("east"), ast.get("north")

                eastTrue = angle
                northTrue = angle + math.radians(90)
                while east < 0: east += 2.0 * math.pi
                while north < 0: north += 2.0 * math.pi
                while northTrue > 2.0 * math.pi: northTrue -= 2.0 * math.pi

                self.assertAlmostEqual(east, eastTrue)
                self.assertAlmostEqual(north, northTrue)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(RotationAngleTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
