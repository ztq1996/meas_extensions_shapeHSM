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
import numpy as np
import unittest
import itertools

import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as algorithms
import lsst.utils.tests as utilsTests
import lsst.afw.detection as afwDetection
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.coord as afwCoord
import lsst.afw.display.ds9 as ds9

import lsst.meas.extensions.shapeHSM

try:
    type(verbose)
except NameError:
    verbose = 0
    display = False

SIZE_DECIMALS = 2 # Number of decimals for equality in sizes
SHAPE_DECIMALS = 3 # Number of decimals for equality in shapes

### The following values are pulled directly from GalSim's test_hsm.py:
file_indices = [0, 2, 4, 6, 8]
x_centroid = [35.888, 19.44, 8.74, 20.193, 57.94]
y_centroid = [19.845, 25.047, 11.92, 38.93, 27.73]
sky_var = [35.01188, 35.93418, 35.15456, 35.11146, 35.16454]
correction_methods = ["KSB", "BJ", "LINEAR", "REGAUSS"]
# Note: expected results give shear for KSB and distortion for others, but the results below have
# converted KSB expected results to distortion for the sake of consistency
e1_expected = np.array([
        [0.467603106752, 0.381211727, 0.398856937, 0.401755571],
        [0.28618443944, 0.199222784, 0.233883543, 0.234257525],
        [0.271533794146, 0.158049396, 0.183517068, 0.184893412],
        [-0.293754156071, -0.457024541, 0.123946584, -0.609233462],
        [0.557720893779, 0.374143023, 0.714147448, 0.435404409] ])
e2_expected = np.array([
        [-0.867225166489, -0.734855778, -0.777027588, -0.774684891],
        [-0.469354341577, -0.395520479, -0.502540961, -0.464466257],
        [-0.519775291311, -0.471589061, -0.574750641, -0.529664935],
        [0.345688365839, -0.342047099, 0.120603755, -0.44609129428863525],
        [0.525728304099, 0.370691830, 0.702724807, 0.433999442] ])
resolution_expected = np.array([
        [0.796144249, 0.835624917, 0.835624917, 0.827796187],
        [0.685023735, 0.699602704, 0.699602704, 0.659457638],
        [0.634736458, 0.651040481, 0.651040481, 0.614663396],
        [0.477027015, 0.477210752, 0.477210752, 0.423157447],
        [0.595205998, 0.611824797, 0.611824797, 0.563582092] ])
sigma_e_expected = np.array([
        [0.016924826, 0.014637648, 0.014637648, 0.014465546],
        [0.075769504, 0.073602324, 0.073602324, 0.064414520],
        [0.110253112, 0.106222900, 0.106222900, 0.099357106],
        [0.185276702, 0.184300955, 0.184300955, 0.173478300],
        [0.073020065, 0.070270966, 0.070270966, 0.061856263] ])
### End of GalSim's values

# These values calculated using GalSim's HSM as part of GalSim
moments_expected = np.array([ # sigma, e1, e2
    [2.24490427971, 0.336240686301, -0.627372910656],
    [1.9031778574, 0.150566105384, -0.245272792302],
    [1.77790760994, 0.112286123389, -0.286203939641],
    [1.45464873314, -0.155597168978, -0.102008266223],
    [1.63144648075, 0.22886961923, 0.228813588897],
    ])
centroid_expected = np.array([ # x, y
    [36.218247328, 20.5678722157],
    [20.325744838, 25.4176650386],
    [9.54257706283, 12.6134786199],
    [20.6407850048, 39.5864802706],
    [58.5008586442, 28.2850942049],
    ])

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ShapeTestCase(unittest.TestCase):
    """A test case for shape measurement"""

    def setUp(self):

        # load the known values
        self.dataDir = os.path.join(os.getenv('MEAS_EXTENSIONS_SHAPEHSM_DIR'), "tests", "data")
        self.bkgd = 1000.0 # standard for atlas image

    def runMeasurement(self, algorithmName, imageid, x, y):
        """Run the measurement algorithm on an image"""
        # load the test image
        imgFile = os.path.join(self.dataDir, "image.%d.fits" % imageid)
        img = afwImage.ImageF(imgFile)
        img -= self.bkgd
        nx, ny = img.getWidth(), img.getHeight()
        msk = afwImage.MaskU(afwGeom.Extent2I(nx, ny), 0x0)
        var = afwImage.ImageF(imgFile)
        mimg = afwImage.MaskedImageF(img, msk, var)
        msk.getArray()[:] = np.where(np.fabs(img.getArray()) < 1.0e-8, msk.getPlaneBitMask("BAD"), 0)

        exposure = afwImage.makeExposure(mimg)
        exposure.setWcs(afwImage.makeWcs(afwCoord.makeCoord(afwCoord.ICRS, 0. * afwGeom.degrees, 0. * afwGeom.degrees),
                                         afwGeom.Point2D(1.0,1.0),
                                         1.0/(2.53*3600.0), 0.0, 0.0, 1.0/(2.53*3600.0)))


        # load the corresponding test psf
        psfFile = os.path.join(self.dataDir, "psf.%d.fits" % imageid)
        psfImg = afwImage.ImageD(psfFile)
        psfImg -= self.bkgd

        kernel = afwMath.FixedKernel(psfImg)
        kernelPsf = algorithms.KernelPsf(kernel)
        exposure.setPsf(kernelPsf)


        # perform the shape measurement
        schema = afwTable.SourceTable.makeMinimalSchema()
        msConfig = algorithms.SourceMeasurementConfig()
        msConfig.algorithms.names = [algorithmName]
        shapeFinder = msConfig.makeMeasureSources(schema)

        # Note: It is essential to remove the floating point part of the position for the
        # Algorithm._apply.  Otherwise, when the PSF is realised it will have been warped
        # to account for the sub-pixel offset and we won't get *exactly* this PSF.
        center = afwGeom.Point2D(int(x), int(y))
        table = afwTable.SourceTable.make(schema)
        source = table.makeRecord()
        source.setFootprint(afwDetection.Footprint(exposure.getBBox()))

        shapeFinder.apply(source, exposure, center)

        return source

    def testHsmShape(self):
        """Test that we can instantiate and play with a measureShape"""

        nFail = 0
        msg = ""

        for (algNum, algName), (i, imageid) in itertools.product(enumerate(correction_methods),
                                                                 enumerate(file_indices)):
            algorithmName = "shape.hsm." + algName.lower()
            source = self.runMeasurement(algorithmName, imageid, x_centroid[i], y_centroid[i])
            
            ##########################################
            # see how we did

            if algName in ("KSB"):
                # Need to convert g1,g2 --> e1,e2 because GalSim has done that
                # for the expected values ("for consistency")
                g1 = source.get(algorithmName + ".g1")
                g2 = source.get(algorithmName + ".g2")
                scale = 2.0/(1.0 + g1**2 + g2**2)
                e1 = g1*scale
                e2 = g2*scale
            else:
                e1 = source.get(algorithmName + ".e1")
                e2 = source.get(algorithmName + ".e2")
            sigma = source.get(algorithmName + ".sigma")
            resolution = source.get(algorithmName + ".resolution")
            flags = source.get(algorithmName + ".flags")
                
            tests = [
                # label        known-value                            measured              tolerance
                ["e1",         float(e1_expected[i][algNum]),         e1,             0.5*10**-SHAPE_DECIMALS],
                ["e2",         float(e2_expected[i][algNum]),         e2,             0.5*10**-SHAPE_DECIMALS],
                ["resolution", float(resolution_expected[i][algNum]), resolution,     0.5*10**-SIZE_DECIMALS],
                
                # sigma won't match exactly because
                # we're using skyvar=mean(var) instead of measured value ... expected a difference
                ["sigma",      float(sigma_e_expected[i][algNum]),    sigma,             0.07],
                ["shapeStatus", 0,                          flags,             0],
                ]

            
            for test in tests:
                label, know, hsm, limit = test
                err = hsm - know
                msgTmp = "%-12s %s  %5s:   %6.6f %6.6f  (val-known) = %.3g\n" % (algName, imageid,
                                                                                 label, know, hsm, err)
                if not np.isfinite(err) or abs(err) > limit:
                    msg += msgTmp
                    nFail += 1

        self.assertEqual(nFail, 0, "\n"+msg)
        
    def testHsmMoments(self):
        for (i, imageid) in enumerate(file_indices):
            source = self.runMeasurement("shape.hsm.moments", imageid, x_centroid[i], y_centroid[i])
            moments = source.get("shape.hsm.moments")
            centroid = source.get("shape.hsm.moments.centroid")

            self.assertAlmostEqual(centroid[0], centroid_expected[i][0], 3)
            self.assertAlmostEqual(centroid[1], centroid_expected[i][1], 3)

            expected = afwEll.Quadrupole(afwEll.SeparableDistortionDeterminantRadius(
                moments_expected[i][1], moments_expected[i][2], moments_expected[i][0]))

            self.assertAlmostEqual(moments.getIxx(), expected.getIxx(), SHAPE_DECIMALS)
            self.assertAlmostEqual(moments.getIxy(), expected.getIxy(), SHAPE_DECIMALS)
            self.assertAlmostEqual(moments.getIyy(), expected.getIyy(), SHAPE_DECIMALS)


        
        
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
