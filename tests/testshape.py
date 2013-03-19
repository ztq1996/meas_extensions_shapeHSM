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

import re, os, sys
import glob
import math
import unittest

import lsst.pex.policy as pexPolicy
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests
import lsst.afw.detection as afwDetection
import lsst.afw.table as afwTable

import lsst.afw.display.ds9 as ds9

import lsst.meas.extensions.shapeHSM as hsm

try:
    type(verbose)
except NameError:
    verbose = 0
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class HsmMomentTestCase(unittest.TestCase):
    """A test case for the HSM shape moments (un-psf-corrected)."""

    def setUp(self):
        
        self.bkgd = 100.0
        self.x, self.y = 30, 40
        self.I0 = 1000.0

        # make a masked image and an exposure
        im = afwImage.ImageF(afwGeom.Extent2I(100, 100))
        im.set(self.bkgd)
        im.set(self.x, self.y, self.I0 + self.bkgd)
        msk = afwImage.MaskU(im.getDimensions()); msk.set(0)
        var = afwImage.ImageF(im.getDimensions()); var.set(math.sqrt(self.I0))
        self.mimg = afwImage.MaskedImageF(im, msk, var)
        del im; del msk; del var
        self.exposure = afwImage.makeExposure(self.mimg)
        
        kwid = 35
        sigma = 2.5
        self.psf = algorithms.SingleGaussianPsf(kwid, kwid, sigma)

        self.exposure.setPsf(self.psf)
        
        # Add a Gaussian to the image
        self.sigma_xx, self.sigma_xy, self.sigma_yy, self.ksize = math.pow(1.5, 2), 0, math.pow(2.5, 2), 15

        phi = math.radians(0.0)      # rotate +ve this far
        c, s = math.cos(phi), math.sin(phi)

        self.sum, self.sumxx, self.sumxy, self.sumyy = 4*[0.0]
        for dx in range(-(self.ksize/2), self.ksize/2 + 1):
            for dy in range(-(self.ksize/2), self.ksize/2 + 1):
                u, v = c*dx + s*dy,  s*dx - c*dy
                I = self.I0*math.exp(-0.5*(u*u/self.sigma_xx + v*v/self.sigma_yy))
                self.mimg.getImage().set(self.x + dx, self.y + dy, self.bkgd + I)

                self.sum += I
                self.sumxx += I*dx*dx
                self.sumxy += I*dx*dy
                self.sumyy += I*dy*dy

        self.sumxx /= self.sum
        self.sumxy /= self.sum
        self.sumyy /= self.sum


        self.policy = pexPolicy.Policy(pexPolicy.PolicyString(
                """#<?cfg paf policy?>
            HSM_SHAPELET.max_order_psf: 8
            HSM_SHAPELET.max_order_gal: 8
            """
                ))
        
        print "%g %g %g" % (self.sumxx, self.sumxy, self.sumyy)


        
    def tearDown(self):
        del self.mimg
        del self.exposure
        del self.psf
        del self.policy


        
    def testHsmMoments(self):
        """Test that the un-psf-corrected moments are correct"""

        # measure shape with our algorithm
        algorithmNames = ["shape.hsm.bj", "shape.hsm.linear", "shape.hsm.ksb",
                          "shape.hsm.regauss", "shape.hsm.shapelet"]

        mim = self.exposure.getMaskedImage()
        im  = mim.getImage()
        im -= self.bkgd
        
        for algName in algorithmNames:
            schema = afwTable.SourceTable.makeMinimalSchema()
            msConfig = algorithms.SourceMeasurementConfig()
            msConfig.algorithms.names = [algName]
            shapeFinder = msConfig.makeMeasureSources(schema)

            if display:
                ds9.mtv(self.mimg)

            center = afwGeom.Point2D(self.x, self.y)
            table = afwTable.SourceTable.make(schema)
            source = table.makeRecord()
            source.setFootprint(afwDetection.Footprint(self.exposure.getBBox()))

            shapeFinder.apply(source, self.exposure, center)

            p = source.get(algName + ".centroid")
            s = source.get(algName + ".moments")

            # see how we did
            Ixx, Iyy, Ixy = s.getIxx(), s.getIyy(), s.getIxy()
            A2 = 0.5*(Ixx + Iyy) + math.sqrt( (0.5*(Ixx - Iyy))**2 + Ixy**2 )
            B2 = 0.5*(Ixx + Iyy) - math.sqrt( (0.5*(Ixx - Iyy))**2 + Ixy**2 )

            print "I_xx:  %.5f %.5f" % (Ixx, self.sigma_xx)
            print "I_xy:  %.5f %.5f" % (Ixy, self.sigma_xy)
            print "I_yy:  %.5f %.5f" % (Iyy, self.sigma_yy)
            print "A2, B2 = %.5f, %.5f" % (A2, B2)         

            # 1/100 pixel limit
            self.assertTrue(abs(self.x - p.getX()) < 1e-2, "%g v. %g" % (self.x, p.getX()))
            self.assertTrue(abs(self.y - p.getY()) < 1e-2, "%g v. %g" % (self.y, p.getY()))

            # need higher limit for these (they're just placeholder shapes computed, after all)
            self.assertTrue(abs(s.getIxx() - self.sigma_xx) < 2e-2*(1 + self.sigma_xx),
                            "%g v. %g" % (self.sigma_xx, s.getIxx()))
            self.assertTrue(abs(s.getIxy() - self.sigma_xy) < 2e-2*(1 + self.sigma_xy),
                            "%g v. %g" % (self.sigma_xy, s.getIxy()))
            self.assertTrue(abs(s.getIyy() - self.sigma_yy) < 3e-2*(1 + self.sigma_yy),
                            "%g v. %g" % (self.sigma_yy, s.getIyy()))


        
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(HsmMomentTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
