import numpy as np
import xgc
import math
import pylab as py
#import xgcjrm as xgc
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
import sys
import cProfile
#import loader_mod
import core
#import Flux_utilities
#import E_field_utilities
#import Profile_SOL_utilities



def Radial_E_field(iz,it,cut):
    global Zi,Ri,pot
    option=1
    if cut==None:
        cut = core.jcut
    dZ = core.Zi[1]-core.Zi[0]
    dR = core.Ri[1]-core.Ri[0]
    Sep_ip = core.getiparrs_hor(0.95,1.,core.jcut)[-1]

    # dphi/dR
    dphidR_N=[]
    lengthR = len(core.getcutvalue_hor(core.pot[:,iz,it],cut,option))
    iterlistR = list(range(0,lengthR))
    potrm10 = core.getcutvalue_hor(core.pot[:,iz,it],cut,option)[lengthR-1]
    potrp10 = core.getcutvalue_hor(core.pot[:,iz,it],cut,option)[1]
    dphidR_N.append(-(potrp10 - potrm10)/(2.*dR))
    for i in iterlistR[1:-1]:
        potrm1 = core.getcutvalue_hor(core.pot[:,iz,it],cut,option)[i-1]
        potrp1 = core.getcutvalue_hor(core.pot[:,iz,it],cut,option)[i+1]
        dphidR_N.append(-(potrp1 - potrm1)/(2.*dR))
    potrm1L = core.getcutvalue_hor(core.pot[:,iz,it],cut,option)[lengthR-2]
    potrp1L = core.getcutvalue_hor(core.pot[:,iz,it],cut,option)[0]
    dphidR_N.append(-(potrp1L-potrm1L)/2*dR)
    
    #Plotting
    L = np.array([x for x in range(0,len(dphidR_N))])
    unit = (core.Rmax-core.Rmin)/100
    R = np.array([L[ip]*unit+core.Rmin for ip in range(0,len(dphidR_N))])
    plt.title('Radial E-field')
    plt.xlabel('$R(m)$')
    plt.ylabel('$E_r$')
    plt.plot(R[:], dphidR_N[:],'go')
    plt.axvline(x=Sep_ip*unit + core.Rmin, color='b', linestyle='dashed', linewidth=2)
    plt.legend()
    plt.show()



def Local_E_Z(iz,it,cut):
    global Zi,Ri,pot
    option=1
    if cut==None:
        cut = core.jcut
    dZ = core.Zi[1]-core.Zi[0]
    dR = core.Ri[1]-core.Ri[0]
    Sep_ip = core.getiparrs_hor(0.95,1.,core.jcut)[-1]

   
    # dphi/dZ
    potim1 = core.getcutvalue_hor(core.pot[:,iz,it],cut-1,option)
    poti = core.getcutvalue_hor(core.pot[:,iz,it],cut,option)
    potip1 = core.getcutvalue_hor(core.pot[:,iz,it],cut+1,option)
    dphidZ_N = -(potip1-potim1)/(2.*dZ)
    
    #Plotting
    L = np.array([x for x in range(0,len(dphidZ_N))])
    unit = (core.Rmax-core.Rmin)/100
    R = np.array([L[ip]*unit+core.Rmin for ip in range(0,len(dphidZ_N))])
    plt.title('Local $E_y$-field')
    plt.xlabel('$R(m)$')
    plt.ylabel('$E_y$')
    plt.plot(R[:], dphidZ_N[:],'go')
    plt.axvline(x=Sep_ip*unit + core.Rmin, color='b', linestyle='dashed', linewidth=2)
    plt.legend()
    plt.show()


def z_variation_E_z_field(iz,it,cut):
    global Zi,Ri,pot
    option=1
    if cut==None:
        cut = core.jcut
    dZ = core.Zi[1]-core.Zi[0]
    dR = core.Ri[1]-core.Ri[0]
    Sep_ip = core.getiparrs_hor(0.95,1.,core.jcut)[-1]
    
    field_value = []
    for cut_plane in range(cut-20,cut+20):
        # dphi/dZ
        potim1 = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane-1,option)
        poti = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane,option)
        potip1 = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane+1,option)
        dphidZ_N = -(potip1-potim1)/(2.*dZ)
        field_value.append(dphidZ_N[Sep_ip-1])

    #Plotting
    L = np.array([x for x in range(0,len(field_value))])
    unit = (core.Zmax-core.Zmin)/100
    R = np.array([L[ip]*unit+core.Zmin for ip in range(0,len(field_value))])
    plt.title('Local $E_y$-field at a vertical line')
    plt.xlabel('$Z(m)$')
    plt.ylabel('$E_y$')
    plt.plot(R[:], field_value[:],'go')
    plt.legend()
    plt.show()


def z_variation_E_r_field(iz,it,cut):
    global Zi,Ri,pot
    option=1
    if cut==None:
        cut = core.jcut
    dZ = core.Zi[1]-core.Zi[0]
    dR = core.Ri[1]-core.Ri[0]
    Sep_ip = core.getiparrs_hor(0.95,1.,core.jcut)[-1]
    
    field_value = []
    for cut_plane in range(cut-20,cut+20):
        # dphi/dR
        dphidR_N=[]
        lengthR = len(core.getcutvalue_hor(core.pot[:,iz,it],cut_plane,option))
        iterlistR = list(range(0,lengthR))
        potrm10 = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane,option)[lengthR-1]
        potrp10 = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane,option)[1]
        dphidR_N.append(-(potrp10 - potrm10)/(2.*dR))
        for i in iterlistR[1:-1]:
            potrm1 = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane,option)[i-1]
            potrp1 = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane,option)[i+1]
            dphidR_N.append(-(potrp1 - potrm1)/(2.*dR))
        potrm1L = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane,option)[lengthR-2]
        potrp1L = core.getcutvalue_hor(core.pot[:,iz,it],cut_plane,option)[0]
        dphidR_N.append(-(potrp1L-potrm1L)/2*dR)
        field_value.append(dphidR_N[Sep_ip-1])
        
    #Plotting
    L = np.array([x for x in range(0,len(field_value))])
    unit = (core.Zmax-core.Zmin)/100
    R = np.array([L[ip]*unit+core.Zmin for ip in range(0,len(field_value))])
    plt.title('Local $E_r$-field at a vertical line')
    plt.xlabel('$Z(m)$')
    plt.ylabel('$E_r$')
    plt.plot(R[:], field_value[:],'go')
    plt.legend()
    plt.show()
