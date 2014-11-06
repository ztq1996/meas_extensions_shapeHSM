# Enable HSM measurements
#
# The 'root' should be a SourceMeasurementConfig.
#
# We activate the REGAUSS PSF-corrected shape measurement, and the adaptive moments of the source and PSF.

import lsst.meas.extensions.shapeHSM
root.algorithms.names |= ["shape.hsm.regauss", "shape.hsm.moments", "shape.hsm.psfMoments"]
