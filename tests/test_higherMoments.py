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
import lsst.utils.tests as tests
import itertools


MOMENTS_DECIMAL = 3  # Number of decimals for equality in moments



class HigherMomentsTestCase(tests.TestCase):
    """A test case for shape measurement"""


    def runMeasurement(self, testname):
        """Create an exposure and run measurement on the source and the PSF"""
        
        
        

        
        if testname == "point_source+gaussian_PSF":
            # Initialize a config and activate the plugin
            sfmConfig = SingleFrameMeasurementConfig()
            sfmConfig.plugins.names |= ["ext_shapeHSM_HigherOrderMomentsSource", "ext_shapeHSM_HigherOrderMomentsPSF"]
            sfmConfig.plugins['ext_shapeHSM_HigherOrderMomentsPSF'].orderMax = 6
            sfmConfig.plugins['ext_shapeHSM_HigherOrderMomentsSource'].orderMax = 6
            # Create a minimal schema (columns)
            schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()

            # Create a task
            sfmTask = SingleFrameMeasurementTask(config=sfmConfig, schema=schema)
            # Create a simple, fake dataset
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
            dataset = lsst.meas.base.tests.TestDataset(bbox)
            # Create a point source with Gaussian PSF
            dataset.addSource(100000.0, lsst.geom.Point2D(49.5, 49.5))
            
            # Create a galaxy with Gaussian PSF
            dataset.addSource(300000.0, lsst.geom.Point2D(76.3, 79.2), lsst.afw.geom.Quadrupole(2.0, 3.0, 0.5))

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
        
        for row in catalog1:
            
            M_source_30 = row["ext_shapeHSM_HigherOrderMomentsSource_30"]
            M_source_21 = row["ext_shapeHSM_HigherOrderMomentsSource_21"]
            M_source_12 = row["ext_shapeHSM_HigherOrderMomentsSource_12"]
            M_source_03 = row["ext_shapeHSM_HigherOrderMomentsSource_03"]
        
            M_source_40 = row["ext_shapeHSM_HigherOrderMomentsSource_40"]
            M_source_31 = row["ext_shapeHSM_HigherOrderMomentsSource_31"]
            M_source_22 = row["ext_shapeHSM_HigherOrderMomentsSource_22"]
            M_source_13 = row["ext_shapeHSM_HigherOrderMomentsSource_13"]
            M_source_04 = row["ext_shapeHSM_HigherOrderMomentsSource_04"]
            
            M_source_50 = row["ext_shapeHSM_HigherOrderMomentsSource_50"]
            M_source_41 = row["ext_shapeHSM_HigherOrderMomentsSource_41"]
            M_source_32 = row["ext_shapeHSM_HigherOrderMomentsSource_32"]
            M_source_23 = row["ext_shapeHSM_HigherOrderMomentsSource_23"]
            M_source_14 = row["ext_shapeHSM_HigherOrderMomentsSource_14"]
            M_source_05 = row["ext_shapeHSM_HigherOrderMomentsSource_05"]
            
            M_source_60 = row["ext_shapeHSM_HigherOrderMomentsSource_60"]
            M_source_51 = row["ext_shapeHSM_HigherOrderMomentsSource_51"]
            M_source_42 = row["ext_shapeHSM_HigherOrderMomentsSource_42"]
            M_source_33 = row["ext_shapeHSM_HigherOrderMomentsSource_33"]
            M_source_24 = row["ext_shapeHSM_HigherOrderMomentsSource_24"]
            M_source_15 = row["ext_shapeHSM_HigherOrderMomentsSource_15"]
            M_source_06 = row["ext_shapeHSM_HigherOrderMomentsSource_06"]
            
            
            self.assertAlmostEqual(M_source_30, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_21, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_12, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_03, 0.0, MOMENTS_DECIMAL)
            
            
            self.assertAlmostEqual(M_source_40, 0.75, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_31, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_22, 0.25, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_13, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_04, 0.75, MOMENTS_DECIMAL)
            
            self.assertAlmostEqual(M_source_50, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_41, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_32, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_23, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_14, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_05, 0.0, MOMENTS_DECIMAL)
            
            self.assertAlmostEqual(M_source_60, 1.875, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_51, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_42, 0.375, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_33, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_24, 0.375, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_15, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_06, 1.875, MOMENTS_DECIMAL)
            

        
            M_psf_40 = row["ext_shapeHSM_HigherOrderMomentsPSF_40"]
            M_psf_31 = row["ext_shapeHSM_HigherOrderMomentsPSF_31"]
            M_psf_22 = row["ext_shapeHSM_HigherOrderMomentsPSF_22"]
            M_psf_13 = row["ext_shapeHSM_HigherOrderMomentsPSF_13"]
            M_psf_04 = row["ext_shapeHSM_HigherOrderMomentsPSF_04"]

            M_psf_30 = row["ext_shapeHSM_HigherOrderMomentsPSF_30"]
            M_psf_21 = row["ext_shapeHSM_HigherOrderMomentsPSF_21"]
            M_psf_12 = row["ext_shapeHSM_HigherOrderMomentsPSF_12"]
            M_psf_03 = row["ext_shapeHSM_HigherOrderMomentsPSF_03"]
            
            M_psf_50 = row["ext_shapeHSM_HigherOrderMomentsPSF_50"]
            M_psf_41 = row["ext_shapeHSM_HigherOrderMomentsPSF_41"]
            M_psf_32 = row["ext_shapeHSM_HigherOrderMomentsPSF_32"]
            M_psf_23 = row["ext_shapeHSM_HigherOrderMomentsPSF_23"]
            M_psf_14 = row["ext_shapeHSM_HigherOrderMomentsPSF_14"]
            M_psf_05 = row["ext_shapeHSM_HigherOrderMomentsPSF_05"]
            
            M_psf_60 = row["ext_shapeHSM_HigherOrderMomentsPSF_60"]
            M_psf_51 = row["ext_shapeHSM_HigherOrderMomentsPSF_51"]
            M_psf_42 = row["ext_shapeHSM_HigherOrderMomentsPSF_42"]
            M_psf_33 = row["ext_shapeHSM_HigherOrderMomentsPSF_33"]
            M_psf_24 = row["ext_shapeHSM_HigherOrderMomentsPSF_24"]
            M_psf_15 = row["ext_shapeHSM_HigherOrderMomentsPSF_15"]
            M_psf_06 = row["ext_shapeHSM_HigherOrderMomentsPSF_06"]
            

            self.assertAlmostEqual(M_psf_40, 0.75, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_psf_31, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_psf_22, 0.25, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_psf_13, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_psf_04, 0.75, MOMENTS_DECIMAL)
            self.assertFloatsAlmostEqual(M_psf_30, 0.0,atol= 2e-3)
            self.assertFloatsAlmostEqual(M_psf_21, 0.0,atol= 2e-3)
            self.assertFloatsAlmostEqual(M_psf_12, 0.0,atol= 2e-3)
            self.assertFloatsAlmostEqual(M_psf_03, 0.0,atol= 2e-3)
            
            self.assertFloatsAlmostEqual(M_psf_50, 0.0, atol= 8e-3)
            self.assertFloatsAlmostEqual(M_psf_41, 0.0, atol= 8e-3)
            self.assertFloatsAlmostEqual(M_psf_32, 0.0, atol= 8e-3)
            self.assertFloatsAlmostEqual(M_psf_23, 0.0, atol= 8e-3)
            self.assertFloatsAlmostEqual(M_psf_14, 0.0, atol= 8e-3)
            self.assertFloatsAlmostEqual(M_psf_05, 0.0, atol= 8e-3)
            
            self.assertFloatsAlmostEqual(M_psf_60, 1.875, atol= 5e-3)
            self.assertFloatsAlmostEqual(M_psf_51, 0.0, atol= 5e-3)
            self.assertFloatsAlmostEqual(M_psf_42, 0.375, atol= 5e-3)
            self.assertFloatsAlmostEqual(M_psf_33, 0.0, atol= 5e-3)
            self.assertFloatsAlmostEqual(M_psf_24, 0.375, atol= 5e-3)
            self.assertFloatsAlmostEqual(M_psf_15, 0.0, atol= 5e-3)
            self.assertFloatsAlmostEqual(M_psf_06, 1.875, atol= 5e-3)


        self.assertEqual(nFail, 0, "\n"+msg)

