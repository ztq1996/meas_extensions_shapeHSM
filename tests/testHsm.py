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
import numpy
import unittest

import lsst.pex.policy as pexPolicy
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests
import lsst.afw.detection as afwDetection
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.display.ds9 as ds9

import lsst.meas.extensions.shapeHSM

try:
    type(verbose)
except NameError:
    verbose = 0
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ShapeTestCase(unittest.TestCase):
    """A test case for shape measurement"""

    def setUp(self):

        # load the known values
        self.dataDir = os.path.join(os.getenv('MEAS_EXTENSIONS_SHAPEHSM_DIR'), "tests", "data")
        
        algNames = []

        knownValuesFile = os.path.join(self.dataDir, "knownValuesTestHsm.dat")
        self.knownValues = {}
        columnNums = {}
        columnNames = [""]*20
        fp = open(knownValuesFile, 'r')
        for line in fp.readlines():
            if len(line.strip()) == 0:
                continue
            
            fields = line.split()
            if re.search("^#", line):
                columnNums[fields[2]] = int(fields[1])
                columnNames[int(fields[1])] = fields[2]
                continue

            imageid = fields[columnNums['imageid']]
            algorithm = fields[columnNums['algorithm']]
            algNames.append(algorithm)
            name = str(imageid) + '-' + algorithm
            if not self.knownValues.has_key(name):
                self.knownValues[name] = {}

            for i in range(2, len(fields)):
                field = fields[i]
                fname = columnNames[i]
                self.knownValues[name][fname] = field
        fp.close()
                
        self.algNames = list(set(algNames))
        self.bkgd = 1000.0 # standard for atlas image
        
        policyFile = os.path.join(os.getenv("MEAS_EXTENSIONS_SHAPEHSM_DIR"), "policy", "hsmShape.paf")
        self.policy = pexPolicy.Policy(policyFile).get("shape")
        
    def tearDown(self):
        del self.policy

    def testHsmShape(self):
        """Test that we can instantiate and play with a measureShape"""

        nFail = 0
        msg = ""

        for name, known in self.knownValues.items():

            imageid, algName = name.split("-")

            #if not name == "6-REGAUSS":
            #    continue
            
            # load the test image
            imgFile = os.path.join(self.dataDir, "image."+imageid+".fits")
            img = afwImage.ImageF(imgFile)
            img -= self.bkgd
            nx, ny = img.getWidth(), img.getHeight()
            msk = afwImage.MaskU(afwGeom.Extent2I(nx, ny), 0x0)
            var = afwImage.ImageF(imgFile)
            mimg = afwImage.MaskedImageF(img, msk, var)

            # code won't run if we mask the bad pixels ... worrisome.
            if False:
                for i in range(mimg.getWidth()):
                    for j in range(mimg.getHeight()):
                        mpix = mimg.get(i, j)
                        print mpix
                        if abs(mpix[0]) < 1.0e-8:
                            mimg.set(i, j, (mpix[0], msk.getPlaneBitMask("BAD"), mpix[2]))

            exposure = afwImage.makeExposure(mimg)
            exposure.setWcs(afwImage.makeWcs(afwCoord.makeCoord(afwCoord.ICRS, 0. * afwGeom.degrees, 0. * afwGeom.degrees),
                                             afwGeom.Point2D(1.0,1.0),
                                             1.0/(2.53*3600.0), 0.0, 0.0, 1.0/(2.53*3600.0)))

            
            # load the corresponding test psf
            psfFile = os.path.join(self.dataDir, "psf."+imageid+".fits")
            psfImg = afwImage.ImageD(psfFile)
            psfImg -= self.bkgd
            
            kernel = afwMath.FixedKernel(psfImg)
            kernelPsf = algorithms.KernelPsf(kernel)
            exposure.setPsf(kernelPsf)

            
            # perform the shape measurement
            schema = afwTable.SourceTable.makeMinimalSchema()
            msConfig = algorithms.SourceMeasurementConfig()
            algorithmName = "shape.hsm." + algName.lower()
            msConfig.algorithms.names = [algorithmName]
            shapeFinder = msConfig.makeMeasureSources(schema)

            x, y = float(known['x']), float(known['y'])
            x2, y2 = int(x+0.5), int(y+0.5)

            center = afwGeom.Point2D(x2, y2)
            table = afwTable.SourceTable.make(schema)
            source = table.makeRecord()
            source.setFootprint(afwDetection.Footprint(exposure.getBBox()))

            shapeFinder.apply(source, exposure, center)
            
            ##########################################
            # see how we did

            s = source.get(algorithmName + ".moments")
            p = source.get(algorithmName + ".centroid")

            if algName in ("KSB", "SHAPELET"):
                e1 = source.get(algorithmName + ".g1")
                e2 = source.get(algorithmName + ".g2")
                shearsig = source.get(algorithmName + ".err")
            else:
                e1 = source.get(algorithmName + ".e1")
                e2 = source.get(algorithmName + ".e2")
                shearsig = 2.0 * source.get(algorithmName + ".err")
            sigma = source.get(algorithmName + ".sigma")
            resolution = source.get(algorithmName + ".resolution")
            flags = source.get(algorithmName + ".flags")
                
            limit = 1.0e-5
            tests = [
                # label      known-value                 measured              tolerance
                ["sigma",      float(known['sigma']),       sigma,             limit],
                ["e1",         float(known['e1']),          e1,                limit],
                ["e2",         float(known['e2']),          e2,                limit],
                ["resolution", float(known['resolution']),  resolution,        limit],
                
                # shearsig won't match exactly
                # we're using skyvar=mean(var) instead of measured value ... expected a difference
                ["shearsig", float(known["shearsig"]),      shearsig,          0.065],
                ["shapeStatus", 0,                          flags,             0],
                ]

            
            for test in tests:
                label, know, hsm, limit = test
                if know != 0:
                    err = (hsm-know)/know
                else:
                    err = hsm-know
                msgTmp = "%-12s %s  %5s:   %6.6f %6.6f  (val-known)/known = %.3g\n" % (algName, imageid,
                                                                                       label, know, hsm, err)
                if not numpy.isfinite(err) or abs(err) > limit:
                    msg += msgTmp
                    nFail += 1

        self.assertEqual(nFail, 0, "\n"+msg)
        

        
        
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ShapeTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    utilsTests.run(suite(), exit)
 
if __name__ == "__main__":
    run(True)
