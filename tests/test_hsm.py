# This file is part of meas_extensions_shapeHSM.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import itertools
import os
import unittest

import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.geom as geom
import lsst.meas.algorithms as algorithms
import lsst.meas.base as base
import lsst.meas.base.tests
import lsst.meas.extensions.shapeHSM as shapeHSM
import lsst.utils.tests
import numpy as np
from lsst.daf.base import PropertySet
import lsst.pex.config as pexConfig
import galsim

SIZE_DECIMALS = 2  # Number of decimals for equality in sizes
SHAPE_DECIMALS = 3  # Number of decimals for equality in shapes

# The following values are pulled directly from GalSim's test_hsm.py:
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
    [0.557720893779, 0.374143023, 0.714147448, 0.435404409]])
e2_expected = np.array([
    [-0.867225166489, -0.734855778, -0.777027588, -0.774684891],
    [-0.469354341577, -0.395520479, -0.502540961, -0.464466257],
    [-0.519775291311, -0.471589061, -0.574750641, -0.529664935],
    [0.345688365839, -0.342047099, 0.120603755, -0.44609129428863525],
    [0.525728304099, 0.370691830, 0.702724807, 0.433999442]])
resolution_expected = np.array([
    [0.796144249, 0.835624917, 0.835624917, 0.827796187],
    [0.685023735, 0.699602704, 0.699602704, 0.659457638],
    [0.634736458, 0.651040481, 0.651040481, 0.614663396],
    [0.477027015, 0.477210752, 0.477210752, 0.423157447],
    [0.595205998, 0.611824797, 0.611824797, 0.563582092]])
sigma_e_expected = np.array([
    [0.016924826, 0.014637648, 0.014637648, 0.014465546],
    [0.075769504, 0.073602324, 0.073602324, 0.064414520],
    [0.110253112, 0.106222900, 0.106222900, 0.099357106],
    [0.185276702, 0.184300955, 0.184300955, 0.173478300],
    [0.073020065, 0.070270966, 0.070270966, 0.061856263]])
# End of GalSim's values

# These values calculated using GalSim's HSM as part of GalSim
galsim_e1 = np.array([
    [0.399292618036, 0.381213068962, 0.398856908083, 0.401749581099],
    [0.155929282308, 0.199228107929, 0.233882278204, 0.234371587634],
    [0.150018423796, 0.158052951097, 0.183515056968, 0.184561833739],
    [-2.6984937191, -0.457033962011, 0.123932465911, -0.60886412859],
    [0.33959621191, 0.374140143394, 0.713756918907, 0.43560180068],
])
galsim_e2 = np.array([
    [-0.74053555727, -0.734855830669, -0.777024209499, -0.774700462818],
    [-0.25573053956, -0.395517915487, -0.50251352787, -0.464388132095],
    [-0.287168383598, -0.471584022045, -0.574719130993, -0.5296921134],
    [3.1754450798, -0.342054128647, 0.120592080057, -0.446093201637],
    [0.320115834475, 0.370669454336, 0.702303349972, 0.433968126774],
])
galsim_resolution = np.array([
    [0.79614430666, 0.835625052452, 0.835625052452, 0.827822327614],
    [0.685023903847, 0.699601829052, 0.699601829052, 0.659438848495],
    [0.634736537933, 0.651039719582, 0.651039719582, 0.614759743214],
    [0.477026551962, 0.47721144557, 0.47721144557, 0.423227936029],
    [0.595205545425, 0.611821532249, 0.611821532249, 0.563564240932],
])
galsim_err = np.array([
    [0.0169247947633, 0.0146376201883, 0.0146376201883, 0.0144661813974],
    [0.0757696777582, 0.0736026018858, 0.0736026018858, 0.0644160583615],
    [0.110252402723, 0.106222368777, 0.106222368777, 0.0993555411696],
    [0.185278102756, 0.184301897883, 0.184301897883, 0.17346136272],
    [0.0730196461082, 0.0702708885074, 0.0702708885074, 0.0618583671749],
])

moments_expected = np.array([  # sigma, e1, e2
    [2.24490427971, 0.336240686301, -0.627372910656],
    [1.9031778574, 0.150566105384, -0.245272792302],
    [1.77790760994, 0.112286123389, -0.286203939641],
    [1.45464873314, -0.155597168978, -0.102008266223],
    [1.63144648075, 0.22886961923, 0.228813588897],
])
centroid_expected = np.array([  # x, y
    [36.218247328, 20.5678722157],
    [20.325744838, 25.4176650386],
    [9.54257706283, 12.6134786199],
    [20.6407850048, 39.5864802706],
    [58.5008586442, 28.2850942049],
])

round_moments_expected = np.array([  # sigma, e1, e2, flux, x, y
    [2.40270376205, 0.197810277343, -0.372329413891, 3740.22436523, 36.4032272633, 20.4847916447],
    [1.89714717865, 0.046496052295, -0.0987404286861, 776.709594727, 20.2893584046, 25.4230368047],
    [1.77995181084, 0.0416346564889, -0.143147706985, 534.59197998, 9.51994111869, 12.6250775205],
    [1.46549296379, -0.0831127092242, -0.0628845766187, 348.294403076, 20.6242279632, 39.5941625731],
    [1.64031589031, 0.0867517963052, 0.0940798297524, 793.374450684, 58.4728765002, 28.2686937854],
])


def makePluginAndCat(alg, name, control=None, metadata=False, centroid=None, psfflux=None, addFlux=False):
    if control is None:
        control = alg.ConfigClass()
    if addFlux:
        control.addFlux = True
    schema = afwTable.SourceTable.makeMinimalSchema()
    if centroid:
        lsst.afw.table.Point2DKey.addFields(schema, centroid, "centroid", "pixel")
        schema.getAliasMap().set("slot_Centroid", centroid)
    if psfflux:
        base.PsfFluxAlgorithm(base.PsfFluxControl(), psfflux, schema)
        schema.getAliasMap().set("slot_PsfFlux", psfflux)
    if metadata:
        plugin = alg(control, name, schema, PropertySet())
    else:
        plugin = alg(control, name, schema)
    cat = afwTable.SourceCatalog(schema)
    if centroid:
        cat.defineCentroid(centroid)
    return plugin, cat


class MomentsTestCase(unittest.TestCase):
    """A test case for shape measurement"""

    def setUp(self):
        # load the known values
        self.dataDir = os.path.join(os.getenv("MEAS_EXTENSIONS_SHAPEHSM_DIR"), "tests", "data")
        self.bkgd = 1000.0  # standard for atlas image
        self.offset = geom.Extent2I(1234, 1234)
        self.xy0 = geom.Point2I(5678, 9876)

    def tearDown(self):
        del self.offset
        del self.xy0

    def runMeasurement(self, algorithmName, imageid, x, y, v, addFlux=False):
        """Run the measurement algorithm on an image"""
        # load the test image
        imgFile = os.path.join(self.dataDir, "image.%d.fits" % imageid)
        img = afwImage.ImageF(imgFile)
        img -= self.bkgd
        nx, ny = img.getWidth(), img.getHeight()
        msk = afwImage.Mask(geom.Extent2I(nx, ny), 0x0)
        var = afwImage.ImageF(geom.Extent2I(nx, ny), v)
        mimg = afwImage.MaskedImageF(img, msk, var)
        msk.getArray()[:] = np.where(np.fabs(img.getArray()) < 1.0e-8, msk.getPlaneBitMask("BAD"), 0)

        # Put it in a bigger image, in case it matters
        big = afwImage.MaskedImageF(self.offset + mimg.getDimensions())
        big.getImage().set(0)
        big.getMask().set(0)
        big.getVariance().set(v)
        subBig = afwImage.MaskedImageF(big, geom.Box2I(big.getXY0() + self.offset, mimg.getDimensions()))
        subBig.assign(mimg)
        mimg = big
        mimg.setXY0(self.xy0)

        exposure = afwImage.makeExposure(mimg)
        cdMatrix = np.array([1.0 / (2.53 * 3600.0), 0.0, 0.0, 1.0 / (2.53 * 3600.0)])
        cdMatrix.shape = (2, 2)
        exposure.setWcs(
            afwGeom.makeSkyWcs(
                crpix=geom.Point2D(1.0, 1.0), crval=geom.SpherePoint(0, 0, geom.degrees), cdMatrix=cdMatrix
            )
        )

        # load the corresponding test psf
        psfFile = os.path.join(self.dataDir, "psf.%d.fits" % imageid)
        psfImg = afwImage.ImageD(psfFile)
        psfImg -= self.bkgd

        kernel = afwMath.FixedKernel(psfImg)
        kernelPsf = algorithms.KernelPsf(kernel)
        exposure.setPsf(kernelPsf)

        # perform the shape measurement
        msConfig = base.SingleFrameMeasurementConfig()
        msConfig.plugins.names |= [algorithmName]
        control = msConfig.plugins[algorithmName]
        alg = base.SingleFramePlugin.registry[algorithmName].PluginClass
        # NOTE: It is essential to remove the floating point part of the position for the
        # Algorithm._apply.  Otherwise, when the PSF is realised it will have been warped
        # to account for the sub-pixel offset and we won't get *exactly* this PSF.
        plugin, table = makePluginAndCat(
            alg, algorithmName, control, centroid="centroid", metadata=True, addFlux=addFlux
        )
        center = geom.Point2D(int(x), int(y)) + geom.Extent2D(self.offset + geom.Extent2I(self.xy0))
        source = table.makeRecord()
        source.set("centroid_x", center.getX())
        source.set("centroid_y", center.getY())
        source.setFootprint(afwDetection.Footprint(afwGeom.SpanSet(exposure.getBBox(afwImage.PARENT))))
        plugin.measure(source, exposure)

        return source

    def testHsmSourceMoments(self):
        for i, imageid in enumerate(file_indices):
            source = self.runMeasurement(
                "ext_shapeHSM_HsmSourceMoments", imageid, x_centroid[i], y_centroid[i], sky_var[i]
            )
            x = source.get("ext_shapeHSM_HsmSourceMoments_x")
            y = source.get("ext_shapeHSM_HsmSourceMoments_y")
            xx = source.get("ext_shapeHSM_HsmSourceMoments_xx")
            yy = source.get("ext_shapeHSM_HsmSourceMoments_yy")
            xy = source.get("ext_shapeHSM_HsmSourceMoments_xy")

            # Centroids from GalSim use the FITS lower-left corner of 1,1
            offset = self.xy0 + self.offset
            self.assertAlmostEqual(x - offset.getX(), centroid_expected[i][0] - 1, 3)
            self.assertAlmostEqual(y - offset.getY(), centroid_expected[i][1] - 1, 3)

            expected = afwEll.Quadrupole(
                afwEll.SeparableDistortionDeterminantRadius(
                    moments_expected[i][1], moments_expected[i][2], moments_expected[i][0]
                )
            )

            self.assertAlmostEqual(xx, expected.getIxx(), SHAPE_DECIMALS)
            self.assertAlmostEqual(xy, expected.getIxy(), SHAPE_DECIMALS)
            self.assertAlmostEqual(yy, expected.getIyy(), SHAPE_DECIMALS)

    def testHsmSourceMomentsRound(self):
        for i, imageid in enumerate(file_indices):
            source = self.runMeasurement(
                "ext_shapeHSM_HsmSourceMomentsRound",
                imageid,
                x_centroid[i],
                y_centroid[i],
                sky_var[i],
                addFlux=True,
            )
            x = source.get("ext_shapeHSM_HsmSourceMomentsRound_x")
            y = source.get("ext_shapeHSM_HsmSourceMomentsRound_y")
            xx = source.get("ext_shapeHSM_HsmSourceMomentsRound_xx")
            yy = source.get("ext_shapeHSM_HsmSourceMomentsRound_yy")
            xy = source.get("ext_shapeHSM_HsmSourceMomentsRound_xy")
            flux = source.get("ext_shapeHSM_HsmSourceMomentsRound_Flux")

            # Centroids from GalSim use the FITS lower-left corner of 1,1
            offset = self.xy0 + self.offset
            self.assertAlmostEqual(x - offset.getX(), round_moments_expected[i][4] - 1, 3)
            self.assertAlmostEqual(y - offset.getY(), round_moments_expected[i][5] - 1, 3)

            expected = afwEll.Quadrupole(
                afwEll.SeparableDistortionDeterminantRadius(
                    round_moments_expected[i][1], round_moments_expected[i][2], round_moments_expected[i][0]
                )
            )
            self.assertAlmostEqual(xx, expected.getIxx(), SHAPE_DECIMALS)
            self.assertAlmostEqual(xy, expected.getIxy(), SHAPE_DECIMALS)
            self.assertAlmostEqual(yy, expected.getIyy(), SHAPE_DECIMALS)

            self.assertAlmostEqual(flux, round_moments_expected[i][3], SHAPE_DECIMALS)

    def testHsmSourceMomentsVsSdssShape(self):
        # Initialize a config and activate the plugins.
        sfmConfig = base.SingleFrameMeasurementConfig()
        sfmConfig.plugins.names |= ["ext_shapeHSM_HsmSourceMoments", "base_SdssShape"]

        # Create a minimal schema (columns).
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()

        # Instantiate the task.
        sfmTask = base.SingleFrameMeasurementTask(config=sfmConfig, schema=schema)

        # Create a simple, test dataset.
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
        dataset = lsst.meas.base.tests.TestDataset(bbox)

        # First source is a point.
        dataset.addSource(100000.0, lsst.geom.Point2D(49.5, 49.5))

        # Second source is a galaxy.
        dataset.addSource(300000.0, lsst.geom.Point2D(76.3, 79.2), afwGeom.Quadrupole(2.0, 3.0, 0.5))

        # Third source is also a galaxy.
        dataset.addSource(250000.0, lsst.geom.Point2D(28.9, 41.35), afwGeom.Quadrupole(1.8, 3.5, 0.4))

        # Get the exposure and catalog.
        exposure, catalog = dataset.realize(10.0, sfmTask.schema, randomSeed=0)

        # Run the measurement task.
        sfmTask.run(catalog, exposure)
        cat = catalog.asAstropy()

        # Get the moments from the catalog.
        xSdss, ySdss = cat["base_SdssShape_x"], cat["base_SdssShape_y"]
        xxSdss, xySdss, yySdss = cat["base_SdssShape_xx"], cat["base_SdssShape_xy"], cat["base_SdssShape_yy"]
        xHsm, yHsm = cat["ext_shapeHSM_HsmSourceMoments_x"], cat["ext_shapeHSM_HsmSourceMoments_y"]
        xxHsm, xyHsm, yyHsm = (
            cat["ext_shapeHSM_HsmSourceMoments_xx"],
            cat["ext_shapeHSM_HsmSourceMoments_xy"],
            cat["ext_shapeHSM_HsmSourceMoments_yy"],
        )

        # Loop over the sources and check that the moments are the same.
        for i in range(3):
            self.assertAlmostEqual(xSdss[i], xHsm[i], 2)
            self.assertAlmostEqual(ySdss[i], yHsm[i], 2)
            self.assertAlmostEqual(xxSdss[i], xxHsm[i], SHAPE_DECIMALS)
            self.assertAlmostEqual(xySdss[i], xyHsm[i], SHAPE_DECIMALS)
            self.assertAlmostEqual(yySdss[i], yyHsm[i], SHAPE_DECIMALS)


class ShapeTestCase(unittest.TestCase):
    """A test case for shape measurement"""

    def setUp(self):

        # load the known values
        self.dataDir = os.path.join(os.getenv('MEAS_EXTENSIONS_SHAPEHSM_DIR'), "tests", "data")
        self.bkgd = 1000.0  # standard for atlas image
        self.offset = geom.Extent2I(1234, 1234)
        self.xy0 = geom.Point2I(5678, 9876)

    def tearDown(self):
        del self.offset
        del self.xy0

    def runMeasurement(self, algorithmName, imageid, x, y, v):
        """Run the measurement algorithm on an image"""
        # load the test image
        imgFile = os.path.join(self.dataDir, "image.%d.fits" % imageid)
        img = afwImage.ImageF(imgFile)
        img -= self.bkgd
        nx, ny = img.getWidth(), img.getHeight()
        msk = afwImage.Mask(geom.Extent2I(nx, ny), 0x0)
        var = afwImage.ImageF(geom.Extent2I(nx, ny), v)
        mimg = afwImage.MaskedImageF(img, msk, var)
        msk.getArray()[:] = np.where(np.fabs(img.getArray()) < 1.0e-8, msk.getPlaneBitMask("BAD"), 0)

        # Put it in a bigger image, in case it matters
        big = afwImage.MaskedImageF(self.offset + mimg.getDimensions())
        big.getImage().set(0)
        big.getMask().set(0)
        big.getVariance().set(v)
        subBig = afwImage.MaskedImageF(big, geom.Box2I(big.getXY0() + self.offset, mimg.getDimensions()))
        subBig.assign(mimg)
        mimg = big
        mimg.setXY0(self.xy0)

        exposure = afwImage.makeExposure(mimg)
        cdMatrix = np.array([1.0/(2.53*3600.0), 0.0, 0.0, 1.0/(2.53*3600.0)])
        cdMatrix.shape = (2, 2)
        exposure.setWcs(afwGeom.makeSkyWcs(crpix=geom.Point2D(1.0, 1.0),
                                           crval=geom.SpherePoint(0, 0, geom.degrees),
                                           cdMatrix=cdMatrix))

        # load the corresponding test psf
        psfFile = os.path.join(self.dataDir, "psf.%d.fits" % imageid)
        psfImg = afwImage.ImageD(psfFile)
        psfImg -= self.bkgd

        kernel = afwMath.FixedKernel(psfImg)
        kernelPsf = algorithms.KernelPsf(kernel)
        exposure.setPsf(kernelPsf)

        # perform the shape measurement
        msConfig = base.SingleFrameMeasurementConfig()
        alg = base.SingleFramePlugin.registry[algorithmName].PluginClass.AlgClass
        control = base.SingleFramePlugin.registry[algorithmName].PluginClass.ConfigClass().makeControl()
        msConfig.algorithms.names = [algorithmName]
        # Note: It is essential to remove the floating point part of the position for the
        # Algorithm._apply.  Otherwise, when the PSF is realised it will have been warped
        # to account for the sub-pixel offset and we won't get *exactly* this PSF.
        plugin, table = makePluginAndCat(alg, algorithmName, control, centroid="centroid")
        center = geom.Point2D(int(x), int(y)) + geom.Extent2D(self.offset + geom.Extent2I(self.xy0))
        source = table.makeRecord()
        source.set("centroid_x", center.getX())
        source.set("centroid_y", center.getY())
        source.setFootprint(afwDetection.Footprint(afwGeom.SpanSet(exposure.getBBox(afwImage.PARENT))))
        plugin.measure(source, exposure)

        return source

    def testHsmShape(self):
        """Test that we can instantiate and play with a measureShape"""

        nFail = 0
        msg = ""

        for (algNum, algName), (i, imageid) in itertools.product(enumerate(correction_methods),
                                                                 enumerate(file_indices)):
            algorithmName = "ext_shapeHSM_HsmShape" + algName[0:1].upper() + algName[1:].lower()

            source = self.runMeasurement(algorithmName, imageid, x_centroid[i], y_centroid[i], sky_var[i])

            ##########################################
            # see how we did
            if algName in ("KSB"):
                # Need to convert g1,g2 --> e1,e2 because GalSim has done that
                # for the expected values ("for consistency")
                g1 = source.get(algorithmName + "_g1")
                g2 = source.get(algorithmName + "_g2")
                scale = 2.0/(1.0 + g1**2 + g2**2)
                e1 = g1*scale
                e2 = g2*scale
                sigma = source.get(algorithmName + "_sigma")
            else:
                e1 = source.get(algorithmName + "_e1")
                e2 = source.get(algorithmName + "_e2")
                sigma = 0.5*source.get(algorithmName + "_sigma")
            resolution = source.get(algorithmName + "_resolution")
            flags = source.get(algorithmName + "_flag")

            tests = [
                # label        known-value                            measured              tolerance
                ["e1", float(e1_expected[i][algNum]), e1, 0.5*10**-SHAPE_DECIMALS],
                ["e2", float(e2_expected[i][algNum]), e2, 0.5*10**-SHAPE_DECIMALS],
                ["resolution", float(resolution_expected[i][algNum]), resolution, 0.5*10**-SIZE_DECIMALS],

                # sigma won't match exactly because
                # we're using skyvar=mean(var) instead of measured value ... expected a difference
                ["sigma", float(sigma_e_expected[i][algNum]), sigma, 0.07],
                ["shapeStatus", 0, flags, 0],
            ]

            for test in tests:
                label, know, hsm, limit = test
                err = hsm - know
                msgTmp = "%-12s %s  %5s:   %6.6f %6.6f  (val-known) = %.3g\n" % (algName, imageid,
                                                                                 label, know, hsm, err)
                if not np.isfinite(err) or abs(err) > limit:
                    msg += msgTmp
                    nFail += 1

            self.assertAlmostEqual(g1 if algName in ("KSB") else e1, galsim_e1[i][algNum], SHAPE_DECIMALS)
            self.assertAlmostEqual(g2 if algName in ("KSB") else e2, galsim_e2[i][algNum], SHAPE_DECIMALS)
            self.assertAlmostEqual(resolution, galsim_resolution[i][algNum], SIZE_DECIMALS)
            self.assertAlmostEqual(sigma, galsim_err[i][algNum], delta=0.07)

        self.assertEqual(nFail, 0, "\n"+msg)


class PyGaussianPsf(afwDetection.Psf):
    # Like afwDetection.GaussianPsf, but handles computeImage exactly instead of
    # via interpolation.  This is a subminimal implementation.  It works for the
    # tests here but isn't fully functional as a Psf class.

    def __init__(self, width, height, sigma, varyBBox=False, wrongBBox=False):
        afwDetection.Psf.__init__(self, isFixed=not varyBBox)
        self.dimensions = geom.Extent2I(width, height)
        self.sigma = sigma
        self.varyBBox = varyBBox  # To address DM-29863
        self.wrongBBox = wrongBBox  # To address DM-30426

    def _doComputeKernelImage(self, position=None, color=None):
        bbox = self.computeBBox(position, color)
        img = afwImage.Image(bbox, dtype=np.float64)
        x, y = np.ogrid[bbox.minY:bbox.maxY+1, bbox.minX:bbox.maxX+1]
        rsqr = x**2 + y**2
        img.array[:] = np.exp(-0.5*rsqr/self.sigma**2)
        img.array /= np.sum(img.array)
        return img

    def _doComputeImage(self, position=None, color=None):
        bbox = self.computeBBox(position, color)
        if self.wrongBBox:
            # For DM-30426:
            # Purposely make computeImage.getBBox() and computeBBox()
            # inconsistent.  Old shapeHSM code attempted to infer the former
            # from the latter, but was unreliable.  New code infers the former
            # directly, so this inconsistency no longer breaks things.
            bbox.shift(geom.Extent2I(1, 1))
        img = afwImage.Image(bbox, dtype=np.float64)
        y, x = np.ogrid[float(bbox.minY):bbox.maxY+1, bbox.minX:bbox.maxX+1]
        x -= (position.x - np.floor(position.x+0.5))
        y -= (position.y - np.floor(position.y+0.5))
        rsqr = x**2 + y**2
        img.array[:] = np.exp(-0.5*rsqr/self.sigma**2)
        img.array /= np.sum(img.array)
        img.setXY0(geom.Point2I(
            img.getX0() + np.floor(position.x+0.5),
            img.getY0() + np.floor(position.y+0.5)
        ))
        return img

    def _doComputeBBox(self, position=None, color=None):
        # Variable size bbox for addressing DM-29863
        dims = self.dimensions
        if self.varyBBox:
            if position.x > 20.0:
                dims = dims + geom.Extent2I(2, 2)
        return geom.Box2I(geom.Point2I(-dims/2), dims)

    def _doComputeShape(self, position=None, color=None):
        return afwGeom.ellipses.Quadrupole(self.sigma**2, self.sigma**2, 0.0)


class PsfMomentsTestCase(unittest.TestCase):
    """A test case for PSF moments measurement"""

    @staticmethod
    def computeDirectPsfMomentsFromGalSim(
        psf,
        center,
        useSourceCentroidOffset=False
    ):
        """Directly from GalSim."""
        psfBBox = psf.computeImageBBox(center)
        psfSigma = psf.computeShape(center).getTraceRadius()
        if useSourceCentroidOffset:
            psfImage = psf.computeImage(center)
            centroid = center
        else:
            psfImage = psf.computeKernelImage(center)
            psfImage.setXY0(psfBBox.getMin())
            centroid = geom.Point2D(psfBBox.getMin() + psfBBox.getDimensions() // 2)
        bbox = psfImage.getBBox(afwImage.PARENT)
        bounds = galsim.bounds.BoundsI(bbox.getMinX(), bbox.getMaxX(), bbox.getMinY(), bbox.getMaxY())
        image = galsim.Image(psfImage.array, bounds=bounds, copy=False)
        guessCentroid = galsim.PositionD(centroid.x, centroid.y)
        shape = galsim.hsm.FindAdaptiveMom(
            image,
            weight=None,
            badpix=None,
            guess_sig=psfSigma,
            precision=1e-6,
            guess_centroid=guessCentroid,
            strict=True,
            round_moments=False,
            hsmparams=None,
        )
        ellipse = lsst.afw.geom.ellipses.SeparableDistortionDeterminantRadius(
            e1=shape.observed_shape.e1,
            e2=shape.observed_shape.e2,
            radius=shape.moments_sigma,
            normalize=True,  # Fail if |e|>1.
        )
        quad = lsst.afw.geom.ellipses.Quadrupole(ellipse)
        ixx = quad.getIxx()
        iyy = quad.getIyy()
        ixy = quad.getIxy()
        return ixx, iyy, ixy

    @lsst.utils.tests.methodParameters(
        # Make Cartesian product of settings to feed to methodParameters
        **dict(list(zip(
            (kwargs := dict(
                # Increasing the width beyond 4.5 leads to noticeable
                # truncation of the PSF, i.e. a PSF that is too large for the
                # box. While this truncated state leads to incorrect
                # measurements, it is necessary for testing purposes to
                # evaluate the behavior under these extreme conditions.
                width=(2.0, 3.0, 4.0, 10.0, 40.0, 100.0),
                useSourceCentroidOffset=(True, False),
                varyBBox=(True, False),
                wrongBBox=(True, False),
                center=(
                    (23.0, 34.0),  # various offsets that might cause trouble
                    (23.5, 34.0),
                    (23.5, 34.5),
                    (23.15, 34.25),
                    (22.81, 34.01),
                    (22.81, 33.99),
                    (1.2, 1.3),  # psfImage extends outside exposure; that's okay
                    (-100.0, -100.0),
                    (-100.5, -100.0),
                    (-100.5, -100.5),
                )
            )).keys(),
            zip(*itertools.product(*kwargs.values()))
        )))
    )
    def testHsmPsfMoments(
        self, width, useSourceCentroidOffset, varyBBox, wrongBBox, center
    ):
        psf = PyGaussianPsf(
            35, 35, width,
            varyBBox=varyBBox,
            wrongBBox=wrongBBox
        )
        exposure = afwImage.ExposureF(45, 56)
        exposure.getMaskedImage().set(1.0, 0, 1.0)
        exposure.setPsf(psf)

        # perform the moment measurement
        algorithmName = "ext_shapeHSM_HsmPsfMoments"
        msConfig = base.SingleFrameMeasurementConfig()
        msConfig.algorithms.names = [algorithmName]
        control = msConfig.plugins[algorithmName]
        alg = base.SingleFramePlugin.registry[algorithmName].PluginClass
        self.assertFalse(control.useSourceCentroidOffset)
        control.useSourceCentroidOffset = useSourceCentroidOffset
        plugin, cat = makePluginAndCat(
            alg, algorithmName,
            centroid="centroid",
            control=control,
            metadata=True,
        )
        source = cat.addNew()
        source.set("centroid_x", center[0])
        source.set("centroid_y", center[1])
        offset = geom.Point2I(*center)
        tmpSpans = afwGeom.SpanSet.fromShape(int(width), offset=offset)
        source.setFootprint(afwDetection.Footprint(tmpSpans))
        plugin.measure(source, exposure)
        x = source.get("ext_shapeHSM_HsmPsfMoments_x")
        y = source.get("ext_shapeHSM_HsmPsfMoments_y")
        xx = source.get("ext_shapeHSM_HsmPsfMoments_xx")
        yy = source.get("ext_shapeHSM_HsmPsfMoments_yy")
        xy = source.get("ext_shapeHSM_HsmPsfMoments_xy")

        self.assertFalse(source.get("ext_shapeHSM_HsmPsfMoments_flag"))
        self.assertFalse(source.get("ext_shapeHSM_HsmPsfMoments_flag_no_pixels"))
        self.assertFalse(source.get("ext_shapeHSM_HsmPsfMoments_flag_not_contained"))
        self.assertFalse(source.get("ext_shapeHSM_HsmPsfMoments_flag_parent_source"))

        if width < 4.5:
            # i.e., as long as the PSF is not truncated for our 35x35 box.
            self.assertAlmostEqual(x, 0.0, 3)
            self.assertAlmostEqual(y, 0.0, 3)
            expected = afwEll.Quadrupole(afwEll.Axes(width, width, 0.0))
            self.assertAlmostEqual(xx, expected.getIxx(), SHAPE_DECIMALS)
            self.assertAlmostEqual(xy, expected.getIxy(), SHAPE_DECIMALS)
            self.assertAlmostEqual(yy, expected.getIyy(), SHAPE_DECIMALS)

        # Test schema documentation
        for fieldName in cat.schema.extract("*HsmPsfMoments_[xy]"):
            self.assertEqual(cat.schema[fieldName].asField().getDoc(),
                             "Centroid of the PSF via the HSM shape algorithm")
        for fieldName in cat.schema.extract("*HsmPsfMoments_[xy][xy]*"):
            self.assertEqual(cat.schema[fieldName].asField().getDoc(),
                             "Adaptive moments of the PSF via the HSM shape algorithm")

        # Test that the moments are identical to those obtained directly by
        # GalSim. For `width` > 4.5 where the truncation becomes significant,
        # the answer might not be 'correct' but should remain 'consistent'.
        xxDirect, yyDirect, xyDirect = self.computeDirectPsfMomentsFromGalSim(
            psf,
            geom.Point2D(*center),
            useSourceCentroidOffset=useSourceCentroidOffset,
        )
        self.assertEqual(xx, xxDirect)
        self.assertEqual(yy, yyDirect)
        self.assertEqual(xy, xyDirect)

    @lsst.utils.tests.methodParameters(
        # Make Cartesian product of settings to feed to methodParameters
        **dict(list(zip(
            (kwargs := dict(
                width=(2.0, 3.0, 4.0),
                useSourceCentroidOffset=(True, False),
                varyBBox=(True, False),
                wrongBBox=(True, False),
                center=(
                    (23.0, 34.0),  # various offsets that might cause trouble
                    (23.5, 34.0),
                    (23.5, 34.5),
                    (23.15, 34.25),
                    (22.81, 34.01),
                    (22.81, 33.99),
                )
            )).keys(),
            zip(*itertools.product(*kwargs.values()))
        )))
    )
    def testHsmPsfMomentsDebiased(
        self, width, useSourceCentroidOffset, varyBBox, wrongBBox, center
    ):
        # As a note, it's really hard to actually unit test whether we've
        # succesfully "debiased" these measurements.  That would require a
        # many-object comparison of moments with and without noise.  So we just
        # test similar to the biased moments above.
        var = 1.2
        # As we reduce the flux, our deviation from the expected value
        # increases, so decrease tolerance.
        for flux, decimals in [
            (1e6, 3),
            (1e4, 1),
            (1e3, 0),
        ]:
            psf = PyGaussianPsf(
                35, 35, width,
                varyBBox=varyBBox,
                wrongBBox=wrongBBox
            )
            exposure = afwImage.ExposureF(45, 56)
            exposure.getMaskedImage().set(1.0, 0, var)
            exposure.setPsf(psf)

            algorithmName = "ext_shapeHSM_HsmPsfMomentsDebiased"
            alg = base.SingleFramePlugin.registry[algorithmName].PluginClass

            # perform the shape measurement
            control = lsst.meas.extensions.shapeHSM.HsmPsfMomentsDebiasedConfig()
            self.assertTrue(control.useSourceCentroidOffset)
            self.assertEqual(control.noiseSource, "variance")
            control.useSourceCentroidOffset = useSourceCentroidOffset
            plugin, cat = makePluginAndCat(
                alg,
                algorithmName,
                centroid="centroid",
                psfflux="base_PsfFlux",
                control=control,
                metadata=True,
            )
            source = cat.addNew()
            source.set("centroid_x", center[0])
            source.set("centroid_y", center[1])
            offset = geom.Point2I(*center)
            source.set("base_PsfFlux_instFlux", flux)
            tmpSpans = afwGeom.SpanSet.fromShape(int(width), offset=offset)
            source.setFootprint(afwDetection.Footprint(tmpSpans))

            plugin.measure(source, exposure)
            x = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_x")
            y = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_y")
            xx = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_xx")
            yy = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_yy")
            xy = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_xy")
            for flag in [
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag",
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag_no_pixels",
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag_not_contained",
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag_parent_source",
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag_edge"
            ]:
                self.assertFalse(source.get(flag))

            expected = afwEll.Quadrupole(afwEll.Axes(width, width, 0.0))

            self.assertAlmostEqual(x, 0.0, decimals)
            self.assertAlmostEqual(y, 0.0, decimals)

            T = expected.getIxx() + expected.getIyy()
            self.assertAlmostEqual((xx-expected.getIxx())/T, 0.0, decimals)
            self.assertAlmostEqual((xy-expected.getIxy())/T, 0.0, decimals)
            self.assertAlmostEqual((yy-expected.getIyy())/T, 0.0, decimals)

            # Repeat using noiseSource='meta'.  Should get nearly the same
            # results if BGMEAN is set to `var` above.
            exposure2 = afwImage.ExposureF(45, 56)
            # set the variance plane to something else to ensure we're
            # ignoring it
            exposure2.getMaskedImage().set(1.0, 0, 2*var+1.1)
            exposure2.setPsf(psf)
            exposure2.getMetadata().set("BGMEAN", var)

            control2 = shapeHSM.HsmPsfMomentsDebiasedConfig()
            control2.noiseSource = "meta"
            control2.useSourceCentroidOffset = useSourceCentroidOffset
            plugin2, cat2 = makePluginAndCat(
                alg,
                algorithmName,
                centroid="centroid",
                psfflux="base_PsfFlux",
                control=control2,
                metadata=True,
            )
            source2 = cat2.addNew()
            source2.set("centroid_x", center[0])
            source2.set("centroid_y", center[1])
            offset2 = geom.Point2I(*center)
            source2.set("base_PsfFlux_instFlux", flux)
            tmpSpans2 = afwGeom.SpanSet.fromShape(int(width), offset=offset2)
            source2.setFootprint(afwDetection.Footprint(tmpSpans2))

            plugin2.measure(source2, exposure2)
            x2 = source2.get("ext_shapeHSM_HsmPsfMomentsDebiased_x")
            y2 = source2.get("ext_shapeHSM_HsmPsfMomentsDebiased_y")
            xx2 = source2.get("ext_shapeHSM_HsmPsfMomentsDebiased_xx")
            yy2 = source2.get("ext_shapeHSM_HsmPsfMomentsDebiased_yy")
            xy2 = source2.get("ext_shapeHSM_HsmPsfMomentsDebiased_xy")
            for flag in [
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag",
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag_no_pixels",
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag_not_contained",
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag_parent_source",
                "ext_shapeHSM_HsmPsfMomentsDebiased_flag_edge"
            ]:
                self.assertFalse(source.get(flag))

            # Would be identically equal, but variance input via "BGMEAN" is
            # consumed in c++ as a double, where variance from the variance
            # plane is a c++ float.
            self.assertAlmostEqual(x, x2, 8)
            self.assertAlmostEqual(y, y2, 8)
            self.assertAlmostEqual(xx, xx2, 5)
            self.assertAlmostEqual(xy, xy2, 5)
            self.assertAlmostEqual(yy, yy2, 5)

        # Test schema documentation
        for fieldName in cat.schema.extract("*HsmPsfMoments_[xy]"):
            self.assertEqual(cat.schema[fieldName].asField().getDoc(),
                             "Debiased centroid of the PSF via the HSM shape algorithm")
        for fieldName in cat.schema.extract("*HsmPsfMoments_[xy][xy]*"):
            self.assertEqual(cat.schema[fieldName].asField().getDoc(),
                             "Debiased adaptive moments of the PSF via the HSM shape algorithm")

    testHsmPsfMomentsDebiasedEdgeArgs = dict(
        width=(2.0, 3.0, 4.0),
        useSourceCentroidOffset=(True, False),
        center=(
            (1.2, 1.3),
            (33.2, 50.1)
        )
    )

    @lsst.utils.tests.methodParameters(
        # Make Cartesian product of settings to feed to methodParameters
        **dict(list(zip(
            (kwargs := dict(
                width=(2.0, 3.0, 4.0),
                useSourceCentroidOffset=(True, False),
                center=[
                    (1.2, 1.3),
                    (33.2, 50.1)
                ]
            )).keys(),
            zip(*itertools.product(*kwargs.values()))
        )))
    )
    def testHsmPsfMomentsDebiasedEdge(self, width, useSourceCentroidOffset, center):
        # As we reduce the flux, our deviation from the expected value
        # increases, so decrease tolerance.
        var = 1.2
        for flux, decimals in [
            (1e6, 3),
            (1e4, 2),
            (1e3, 1),
        ]:
            psf = PyGaussianPsf(35, 35, width)
            exposure = afwImage.ExposureF(45, 56)
            exposure.getMaskedImage().set(1.0, 0, 2*var+1.1)
            exposure.setPsf(psf)

            algorithmName = "ext_shapeHSM_HsmPsfMomentsDebiased"
            alg = base.SingleFramePlugin.registry[algorithmName].PluginClass

            # perform the shape measurement
            control = shapeHSM.HsmPsfMomentsDebiasedConfig()
            control.useSourceCentroidOffset = useSourceCentroidOffset
            self.assertEqual(control.noiseSource, "variance")
            plugin, cat = makePluginAndCat(
                alg,
                algorithmName,
                centroid="centroid",
                psfflux="base_PsfFlux",
                control=control,
                metadata=True,
            )
            source = cat.addNew()
            source.set("centroid_x", center[0])
            source.set("centroid_y", center[1])
            offset = geom.Point2I(*center)
            source.set("base_PsfFlux_instFlux", flux)
            tmpSpans = afwGeom.SpanSet.fromShape(int(width), offset=offset)
            source.setFootprint(afwDetection.Footprint(tmpSpans))

            # Edge fails when setting noise from var plane
            with self.assertRaises(base.MeasurementError):
                plugin.measure(source, exposure)

            # Succeeds when noise is from meta
            exposure.getMetadata().set("BGMEAN", var)
            control.noiseSource = "meta"
            plugin, cat = makePluginAndCat(
                alg,
                algorithmName,
                centroid="centroid",
                psfflux="base_PsfFlux",
                control=control,
                metadata=True,
            )
            source = cat.addNew()
            source.set("centroid_x", center[0])
            source.set("centroid_y", center[1])
            offset = geom.Point2I(*center)
            source.set("base_PsfFlux_instFlux", flux)
            tmpSpans = afwGeom.SpanSet.fromShape(int(width), offset=offset)
            source.setFootprint(afwDetection.Footprint(tmpSpans))
            plugin.measure(source, exposure)

            x = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_x")
            y = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_y")
            xx = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_xx")
            yy = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_yy")
            xy = source.get("ext_shapeHSM_HsmPsfMomentsDebiased_xy")
            self.assertFalse(source.get("ext_shapeHSM_HsmPsfMomentsDebiased_flag"))
            self.assertFalse(source.get("ext_shapeHSM_HsmPsfMomentsDebiased_flag_no_pixels"))
            self.assertFalse(source.get("ext_shapeHSM_HsmPsfMomentsDebiased_flag_not_contained"))
            self.assertFalse(source.get("ext_shapeHSM_HsmPsfMomentsDebiased_flag_parent_source"))
            # but _does_ set EDGE flag in this case
            self.assertTrue(source.get("ext_shapeHSM_HsmPsfMomentsDebiased_flag_edge"))

            expected = afwEll.Quadrupole(afwEll.Axes(width, width, 0.0))

            self.assertAlmostEqual(x, 0.0, decimals)
            self.assertAlmostEqual(y, 0.0, decimals)

            T = expected.getIxx() + expected.getIyy()
            self.assertAlmostEqual((xx-expected.getIxx())/T, 0.0, decimals)
            self.assertAlmostEqual((xy-expected.getIxy())/T, 0.0, decimals)
            self.assertAlmostEqual((yy-expected.getIyy())/T, 0.0, decimals)

            # But fails hard if meta doesn't contain BGMEAN
            exposure.getMetadata().remove("BGMEAN")
            plugin, cat = makePluginAndCat(
                alg,
                algorithmName,
                centroid="centroid",
                psfflux="base_PsfFlux",
                control=control,
                metadata=True,
            )
            source = cat.addNew()
            source.set("centroid_x", center[0])
            source.set("centroid_y", center[1])
            offset = geom.Point2I(*center)
            source.set("base_PsfFlux_instFlux", flux)
            tmpSpans = afwGeom.SpanSet.fromShape(int(width), offset=offset)
            source.setFootprint(afwDetection.Footprint(tmpSpans))
            with self.assertRaises(base.FatalAlgorithmError):
                plugin.measure(source, exposure)

    def testHsmPsfMomentsDebiasedBadNoiseSource(self):
        control = shapeHSM.HsmPsfMomentsDebiasedConfig()
        with self.assertRaises(pexConfig.FieldValidationError):
            control.noiseSource = "ACM"


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
