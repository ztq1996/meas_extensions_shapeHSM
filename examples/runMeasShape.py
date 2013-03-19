#!/usr/bin/env python
#
# Original filename: examples/runMeasShape.py
#
# Author: Steve Bickerton
# Email: 
# Date: Mon 2011-07-11 14:20:42
# 
# Summary: 
# 
"""
%prog [options] x y image psf

Calls the HSM shape code through the LSST meas_algorithms code.

algorithms: KSB LINEAR REGAUSS SHAPELET BJ

Sample data in MEAS_EXTENSIONS_HSM_DIR/tests/data/
  - images as:       image.*.fits
  - PSFs as:         psf.*.fits
  - known values in: knownValuesTestHsm.dat (includes x,y centroids)

example to run test data (from build directory):

./examples/runMeasShape.py -a BJ 35.88 19.84 tests/data/image.0.fits tests/data/psf.0.fits

"""

import sys, os, re
import optparse

import lsst.pex.policy        as pexPolicy
import lsst.afw.image         as afwImage
import lsst.meas.algorithms   as algorithms
import lsst.afw.detection     as afwDetection
import lsst.afw.geom          as afwGeom
import lsst.afw.math          as afwMath

import lsst.meas.extensions.shapeHSM.hsmLib as hsmLib


#############################################################
#
# Main body of code
#
#############################################################

def main(x, y, imgFile, psfFile, algorithm, bkgd):

    policyFile = os.path.join(os.getenv("MEAS_EXTENSIONS_SHAPEHSM_DIR"), "policy", "hsmShape.paf")
    policy = pexPolicy.Policy(policyFile).get("shape")
    
    # load the test image
    img = afwImage.ImageF(imgFile)
    img -= bkgd
    nx, ny = img.getWidth(), img.getHeight()
    msk = afwImage.MaskU(afwGeom.Extent2I(nx, ny), 0x0)
    var = afwImage.ImageF(imgFile)
    mimg = afwImage.MaskedImageF(img, msk, var)

    # code won't run if we mask the bad pixels ... worrisome.
    exposure = afwImage.makeExposure(mimg)
    exposure.setWcs(afwImage.makeWcs(afwGeom.Point2D(0.0,0.0), afwGeom.Point2D(1.0,1.0),
                                     1.0/(2.53*3600.0), 0.0, 0.0, 1.0/(2.53*3600.0)))

    # load the corresponding test psf
    psfImg = afwImage.ImageD(psfFile)
    psfImg -= bkgd

    kernel = afwMath.FixedKernel(psfImg)
    kernelPsf = algorithms.KernelPsf(kernel)
    exposure.setPsf(kernelPsf)


    # perform the shape measurement
    shapeFinder = algorithms.makeMeasureShape(None)

    algorithmName = "HSM_" + algorithm
    shapeFinder.addAlgorithm(algorithmName)
    policy.set(algorithmName+".background", bkgd)
    shapeFinder.configure(policy)

    shapeFinder.setImage(exposure)

    x2, y2 = int(x+0.5), int(y+0.5)
    peak = afwDetection.Peak(x2, y2)
    s = shapeFinder.measure(peak).find(algorithmName)

    ##########################################
    # see how we did
    print "e1,e2,shearSig = ", s.getE1(), s.getE2(), 2.0 * s.getE1Err()


    
#############################################################
# end
#############################################################

if __name__ == '__main__':

    ########################################################################
    # command line arguments and options
    ########################################################################
    
    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("-a", "--algorithm", default="REGAUSS",
                      help="specify the algorithm to use (default=%default)")
    parser.add_option("-b", "--bkgd", default=1000.0, type=float,
                      help="specify the background level (default=%default [SDSS atlas image value])")
    opts, args = parser.parse_args()

    if len(args) != 4:
        parser.print_help()
        sys.exit(1)

    x, y, imgFile, psfFile = args
        
    main(float(x), float(y), imgFile, psfFile, opts.algorithm, opts.bkgd)
