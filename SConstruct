# -*- python -*-
#
# Setup our environment
#
import os
import lsst.SConsUtils as scons

# change this to your package name
thisPkg   = "meas_extensions_shapeHSM"

# change this to what you want your python module to be called
# Note: Your swig .i file should be called pythonPkg + ".i"
#  ie. fooLib.i for this example
pythonPkg = "hsmLib"

# this is where you put your fooLib.i file
pythonPath = os.path.join("python", "lsst", "meas", "extensions", "shapeHSM")



############################################################
# You can probably leave everything below here alone.
############################################################
try:
    scons.ConfigureDependentProducts
except AttributeError:
    import lsst.afw.SconsUtils
    scons.ConfigureDependentProducts = lsst.afw.SconsUtils.ConfigureDependentProducts

libs = "meas_algorithms afw daf_base daf_data daf_persistence "
libs += "pex_logging pex_exceptions pex_policy security boost minuit2 utils wcslib"

env = scons.makeEnv(thisPkg, r"$HeadURL:$",
                    scons.ConfigureDependentProducts(thisPkg))
env.libs[thisPkg] +=  env.getlibs(libs)

# stash these values to pass them to the SConscript files
env.thisPkg    = thisPkg
env.pythonPkg  = pythonPkg
env.pythonPath = pythonPath

# farm out all the building to the SConscript files
# lib/SConscript          builds shared object
# python/.../SConscript   does swigging and builds fooLib.py
# tests/SConscript        builds and runs the tests
for d in Split("lib tests examples "+pythonPath):
    SConscript(os.path.join(d, "SConscript"))

scons.CleanTree(r"*~ core *.so *.os *.o *.pyc")

