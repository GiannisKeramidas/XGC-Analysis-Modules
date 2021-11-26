import numpy as np
import core
import angle
#import SRLSF

core.getMeshAndCuts(core.fileDir,core.Rmin,core.Rmax,core.Zmin,core.Zmax)
sh = core.two_d_shear(650)
np.save('shear_mid',sh)
