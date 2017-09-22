# Enable HSM measurements
#
# The 'config' should be a SourceMeasurementConfig.
#
# We activate the REGAUSS PSF-corrected shape measurement, and the adaptive moments of the source and PSF.

import lsst.meas.extensions.shapeHSM
config.plugins.names |= ["ext_shapeHSM_HsmShapeRegauss", "ext_shapeHSM_HsmSourceMoments",
                         "ext_shapeHSM_HsmPsfMoments"]
config.slots.shape = "ext_shapeHSM_HsmSourceMoments"
config.slots.psfShape = "ext_shapeHSM_HsmPsfMoments"

