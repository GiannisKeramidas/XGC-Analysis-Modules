
# coding: utf-8
from datetime import datetime
import numpy as np
import xgc
import math
#import xgcjrm as xgc
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
import sys

startTime = datetime.now()
#limit data to the [Rmin,Rmax,Zmin,Zmax] box, and read only the first two toroidal planes
Rmin= 2.2 #1.9
Rmax= 2.31 #2.2
Zmin= -0.25 #0.45
Zmax= 0.4 #0.7
phi_start=0
phi_end=None

#fileDir='/scratch2/scratchdirs/giannos/particle_pinch_idl'
#fileDir= '/global/cscratch1/sd/u6003/ti253_d3d_Ip_1.5_large'
#fileDir='/global/homes/g/giannos/ti253_d3d_Ip_1.5_large'
#fileDir='/scratch2/scratchdirs/giannos/ti255_d3d_Ip_1.5_med'
fileDir='/scratch1/scratchdirs/giannos/ti255_d3d_Ip_1.5_med'

def getMeshAndCuts(fileDir,nRi=None,nZi=None,Rmin=Rmin,Rmax=Rmax,Zmin=Zmin,Zmax=Zmax):

	global loader,ne,pot,psin,RZ,tri,time,psin1d,psin001d,Te1d,ne1d,pot001d,bfield,pot0,dpot,dpot,pot,Tigc1d,nigc1d,i_T_perp,i_E_para,e_T_perp,e_E_para
	global Nplanes,Ntimes,dt
	global Ri,Zi,RI,ZI,triObj
	global sepInds,psinvalues,jcut,Zcut#,psii
	
	loader=xgc.load(fileDir,Rmin=Rmin,Rmax=Rmax,Zmin=Zmin,Zmax=Zmax,phi_start=0,phi_end=None)
	ne = loader.calcNeTotal()
	pot = loader.calcPotential()
	print('Rmin = ',loader.Rmin)
	print('Rmax = ',loader.Rmax)
	print('Zmin = ',loader.Zmin)
	print('Zmax = ',loader.Zmax)
	print('ne.shape = ',ne.shape)
	
	# (could have imported * from xgc but this way we see explicitly what we are getting)

	# arrays
	psin = loader.psin
	RZ = loader.RZ
	tri = loader.tri
	time = loader.time
	psin1d = loader.psin1d
	psin001d = loader.psin001d
	Te1d = loader.Te1d
	ne1d = loader.ne1d
	pot001d = loader.pot001d
	bfield = loader.bfield
	pot0 = loader.pot0
	dpot = loader.dpot
	pot = loader.pot
        #if hasattr(loader,'i_T_perp'):
        i_T_perp = loader.i_T_perp
        i_E_para = loader.i_E_para
        e_T_perp = loader.e_T_perp
        e_E_para = loader.e_E_para
 
	# scalars
	Zcut = None
	Nplanes = loader.Nplanes
	Ntimes = loader.Ntimes
	dt = loader.dt

	#setup mesh grid
	if nRi == None:
		nRi = 100
	if nZi == None:
		nZi = 100
	Ri = np.linspace(loader.Rmin,loader.Rmax,nRi)
	Zi = np.linspace(loader.Zmin,loader.Zmax,nZi)
	(RI,ZI) = np.meshgrid(Ri,Zi)
	# set up to interpolate using the TriInterpolator class. Should be what tricontourf() uses
	triObj = Triangulation(RZ[:,0],RZ[:,1],tri)

	# get indices of RZ vertices on the separatrix
	sepInds = np.where(np.abs(loader.psin-1.0)<1e-4)[0]

	# flux surface values on the triangular mesh
	psinvalues = np.sort(np.array(list(set(np.round(psin,4))))) # 4 round to 4 decimal places; adjust?
       
	# locate j value where Z = Zcut
	if Zcut == None:
		Zcut = Zi.mean()
	jcut = int(round(np.mean(np.where(np.abs(Zi-Zcut)<1.e-1))))
	#print "Zi[0:5] = ",Zi[0:5]
	#print "Zi[-5:0] = ",Zi[-5:-1]
	#print "Zcut = ",Zcut
	print("jcut = ",jcut)


	# calculate psi(Ri) = psii at the cut
	tci=LinearTriInterpolator(triObj,psin)
	psiout=tci(RI,ZI)
	psii = psiout[jcut,:] # a 1D array of psin values at each Ri along the cut


getMeshAndCuts(fileDir)
print("dt:",datetime.now()-startTime)
