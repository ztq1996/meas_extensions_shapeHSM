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


import sys
sys.path = ['/pscratch/sd/z/ztq1996/DM/meas_extensions_shapeHSM/python'] + sys.path 

from lsst.meas.base import SingleFrameMeasurementConfig, SingleFrameMeasurementTask
import lsst.meas.extensions.shapeHSM as shapeHSM
import lsst.meas.base.tests
import lsst.afw.geom

import os
import numpy as np
import unittest
import itertools


MOMENTS_DECIMAL = 3  # Number of decimals for equality in moments



class HigherMomentsTestCase(unittest.TestCase):
    """A test case for shape measurement"""


    def runMeasurement(self, testname):
        """Create an exposure and run measurement on the source and the PSF"""
        
        
        
        # Initialize a config and activate the plugin
        sfmConfig = SingleFrameMeasurementConfig()
        sfmConfig.plugins.names |= ["ext_shapeHSM_HigherOrderMomentsSource", "ext_shapeHSM_HigherOrderMomentsPSF"]

        # Create a minimal schema (columns)
        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()
        
        # Create a task
        sfmTask = SingleFrameMeasurementTask(config=sfmConfig, schema=schema)
        
        if testname == "point_source+gaussian_PSF":
            # Create a simple, fake dataset
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
            dataset = lsst.meas.base.tests.TestDataset(bbox)
            # Create a point source with Gaussian PSF
            dataset.addSource(100000.0, lsst.geom.Point2D(49.5, 49.5))

            # Get the exposure and catalog.
            exposure, catalog = dataset.realize(0.0, sfmTask.schema, randomSeed=0)

        sfmTask.run(catalog, exposure)
        

        return catalog.asAstropy()

    def testHsmSourceMoments(self):
        """Test that we can instantiate and play with a measureShape"""
        
        nFail = 0
        msg = ""
        
        test1 = "point_source+gaussian_PSF"
        catalog1 = self.runMeasurement(test1)
        
        M_source_40 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_40"][0]
        M_source_31 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_31"][0]
        M_source_22 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_22"][0]
        M_source_13 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_13"][0]
        M_source_04 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_04"][0]
        
        M_source_30 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_30"][0]
        M_source_21 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_21"][0]
        M_source_12 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_12"][0]
        M_source_03 = catalog1["ext_shapeHSM_HigherOrderMomentsSource_03"][0]
        
        self.assertAlmostEqual(M_source_40, 0.75, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_source_31, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_source_22, 0.25, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_source_13, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_source_04, 0.75, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_source_30, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_source_21, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_source_12, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_source_03, 0.0, MOMENTS_DECIMAL)
        
        M_psf_40 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_40"][0]
        M_psf_31 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_31"][0]
        M_psf_22 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_22"][0]
        M_psf_13 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_13"][0]
        M_psf_04 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_04"][0]
        
        M_psf_30 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_30"][0]
        M_psf_21 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_21"][0]
        M_psf_12 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_12"][0]
        M_psf_03 = catalog1["ext_shapeHSM_HigherOrderMomentsPSF_03"][0]
        
        self.assertAlmostEqual(M_psf_40, 0.75, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_psf_31, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_psf_22, 0.25, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_psf_13, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_psf_04, 0.75, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_psf_30, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_psf_21, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_psf_12, 0.0, MOMENTS_DECIMAL)
        self.assertAlmostEqual(M_psf_03, 0.0, MOMENTS_DECIMAL)


        self.assertEqual(nFail, 0, "\n"+msg)

