import numpy as np
import xgc
import math
import pylab as py
import matplotlib.pyplot as plt
import glob
import angle
import numpy.ma as ma
from decimal import Decimal
import core



def phi_tilda(cut):#loads all potential and density matrices that we will take the derivative.
    option = 1
    t_start = 550 
    t_end = 650
    planes = core.Nplanes
    
    #loading density
    n_e=core.loader.calcNeTotal()
    n_tot = np.asarray([[core.getcutvalue_hor(n_e[:,iz,it],(cut),option) for iz in range(0,planes)] for it in range(t_start,t_end)])
    #loading arrays of potential: phi[it,iz,ip]
    phim = np.asarray([[core.getcutvalue_hor(core.pot[:,iz,it],(cut-1)%len(core.Zi),option) for iz in range(0,planes)] for it in range(t_start,t_end)])
    print('array 1 loaded')
    phi = np.asarray([[core.getcutvalue_hor(core.pot[:,iz,it],cut,option) for iz in range(0,planes)] for it in range(t_start,t_end)])
    print('array 2 loaded')
    phip = np.asarray([[core.getcutvalue_hor(core.pot[:,iz,it],(cut+1)%len(core.Zi),option) for iz in range(0,planes)] for it in range(t_start,t_end)])
    print('array 3 loaded')
    
    #toroidal averaging: phi_tor[it,ip]
    phim_tor = phim.mean(axis=1)
    phi_tor = phi.mean(axis=1)
    phip_tor = phip.mean(axis=1)
    n_tor = n_tot.mean(axis=1)

    #time averaging: phi_av[ip]
    phim_av = phim_tor.mean(axis=0)
    phi_av = phi_tor.mean(axis=0)
    phip_av = phip_tor.mean(axis=0)
    n_av = n_tor.mean(axis=0)
    
    #phi_tilda: dphi[it,iz,ip]
    dphim = np.asarray(phim[:,:,:] - phim_av[np.newaxis,np.newaxis,:])
    dphi = np.asarray(phi[:,:,:] - phi_av[np.newaxis,np.newaxis,:])
    dphip = np.asarray(phip[:,:,:] - phip_av[np.newaxis,np.newaxis,:])
    
    return phim,phi,phip,n_tot,n_av,phim_av,phi_av,phip_av


def theta_psi_der(cut):# takes the d/dtheta and d/dpsi derivatives of a phi[it,iz,ip] along a cut. 
    global Zi,Ri,pot,jcut,bfield,B2,unitR
    option=1
    R_len = len(core.Ri)
    if cut==None:
        cut = core.jcut
        
    dZ = core.Zi[1] - core.Zi[0]
    dR = core.Ri[1] - core.Ri[0]
    
    
    phim,phi,phip,n_tot,n_av,phim_av,phi_av,phip_av = phi_tilda(cut)
    
    #d/dZ of total phi
    arrim1 = phim
    arri = phi
    arrip1 = phip    
    darrdZ = (arrip1-arrim1)/(2.*dZ)
    #d/dZ of average phi
    arravm1 = phim_av
    arrav = phi_av
    arravp1 = phip_av    
    darravdZ = (arravp1-arravm1)/(2.*dZ)
    #d/dR of total phi
    arriplus = np.roll(arri,-1,axis=2)
    arriminus = np.roll(arri,1,axis=2)
    darrdR = (arriplus-arriminus)/(2.*dR)
    darrdR[0] = (arri[1]-arri[0])/dR
    darrdR[-1] = (arri[-1]-arri[-2])/dR
    #d/dR of average phi
    arravplus = np.roll(arrav,-1)
    arravminus = np.roll(arrav,1)
    darravdR = (arravplus-arravminus)/(2.*dR)
    darravdR[0] = (arrav[1]-arrav[0])/dR
    darravdR[-1] = (arrav[-1]-arrav[-2])/dR
    
    BRi = core.getcutvalue_hor(core.bfield[:,0],cut,option)
    BZi = core.getcutvalue_hor(core.bfield[:,1],cut,option)
    Bzetai = core.getcutvalue_hor(core.bfield[:,2],cut,option)
    Bp = np.sqrt(np.square(BRi)+np.square(BZi))
    B = np.sqrt(np.square(BRi)+np.square(BZi)+np.square(Bzetai))
    R = np.array([ip*core.unitR+core.Rmin for ip in range(0,R_len)])+core.Rmaj
    
    
    darrdtheta = (1/(Bp))*(BRi*darrdR + BZi*darrdZ)
    darrdpsi = (1/(Bp))*(-BRi*darrdZ + BZi*darrdR)
    darravdtheta = (1/(Bp))*(BRi*darravdR + BZi*darravdZ)
    darravdpsi = (1/(Bp))*(-BRi*darravdZ + BZi*darravdR)
    return -darrdtheta/B, -darrdpsi/B, n_tot, n_av, -darravdtheta/B,-darravdpsi/B

def two_d_theta_psi_der():#creates a 2-D matrix of d/dtheta and d/dpsi derivatives.
    global Zi,Ri,pot,jcut,bfield,B2
    Z_len = len(core.Zi)
    dpotdtheta, dpotdpsi, n_tot,n_av,dpotavdtheta,dpotavdpsi = zip(*[theta_psi_der(cut) for cut in range(0,Z_len)])
    dphidtheta = np.asarray(dpotdtheta)
    dphidpsi = np.asarray(dpotdpsi)#derivative[Z,time,planes,R]
    den = np.asarray(n_tot)
    den_av = np.asarray(n_av)
    dphiavdtheta = np.asarray(dpotavdtheta)
    dphiavdpsi = np.asarray(dpotavdpsi)
    return dphidtheta,dphidpsi,den,den_av,dphiavdtheta,dphiavdpsi

def Reynolds_stress():#calculates the total, equilibrium and turbulent components of the Reynolds stress tensor, -<n E_psi E_theta>/B^2.
    dphidtheta,dphidpsi,den,den_av,dphiavdtheta,dphiavdpsi = two_d_theta_psi_der()
    #Total Reynolds stresses
    Re_tps = den*dphidtheta*dphidpsi
    Re_tt = den*dphidtheta*dphidtheta
    Re_psps = den*dphidpsi*dphidpsi
    Re_tps_tor = Re_tps.mean(axis=2)
    Re_tps_av = Re_tps_tor.mean(axis=1)
    Re_tt_tor = Re_tt.mean(axis=2)
    Re_tt_av = Re_tt_tor.mean(axis=1)
    Re_psps_tor = Re_psps.mean(axis=2)
    Re_psps_av = Re_psps_tor.mean(axis=1)
    #Equilibrium Reynolds stresses
    Re_eq_tps = den_av*dphiavdtheta*dphiavdpsi
    Re_eq_tt = den_av*dphiavdtheta*dphiavdtheta
    Re_eq_psps = den_av*dphiavdpsi*dphiavdpsi
    #Turbulent Reynolds stresses
    Re_tur_tps = Re_tps_av - Re_eq_tps
    Re_tur_tt = Re_tt_av - Re_eq_tt
    Re_tur_psps = Re_psps_av - Re_eq_psps
    
    return Re_tps_av,Re_tt_av,Re_psps_av,Re_eq_tps,Re_eq_tt,Re_eq_psps,Re_tur_tps,Re_tur_tt,Re_tur_psps




   
def R_der(arr):#takes d/dR of array[Z,R]
    global Zi,Ri,pot,jcut,bfield,B2,unitR
    dR = core.Ri[1] - core.Ri[0]
    
    arriplus = np.roll(arr,-1,axis=1)
    arriminus = np.roll(arr,1,axis=1)
    darrdR = (arriplus-arriminus)/(2.*dR)
    darrdR[:,0] = (arr[:,1]-arr[:,0])/dR
    darrdR[:,-1] = (arr[:,-1]-arr[:,-2])/dR

    return darrdR
 
def Z_der(arr,cut):#takes d/dZ of array[Z,R]
    global Zi,Ri,pot,jcut,bfield,B2,unitR
    option=1
    dZ = core.Zi[1] - core.Zi[0]
    Z_dim = arr.shape[0]
    
    arrim1 = arr[(cut-1)%Z_dim,:]
    arri = arr[cut,:]
    arrip1 =  arr[(cut+1)%Z_dim,:]
    darrdZ = (arrip1-arrim1)/(2.*dZ)
    
    return darrdZ

def thita_psi_der(arr,cut):# takes the d/dtheta and d/dpsi derivatives of a [Z,R] array along a cut. Returns simple derivatives and divided by B. 
    global Zi,Ri,pot,jcut,bfield,B2,unitR
    option=1
    
    darrdR = R_der(arr)[cut,:]
    darrdZ = Z_der(arr,cut)
    
    BRi = core.getcutvalue_hor(core.bfield[:,0],cut,option)
    BZi = core.getcutvalue_hor(core.bfield[:,1],cut,option)
    Bzetai = core.getcutvalue_hor(core.bfield[:,2],cut,option)
    Bp = np.sqrt(np.square(BRi)+np.square(BZi))
    B = np.sqrt(np.square(BRi)+np.square(BZi)+np.square(Bzetai))
    
    darrdtheta = (1/(Bp))*(-BRi*darrdZ + BZi*darrdR)
    darrdpsi = (1/(Bp))*(BRi*darrdR + BZi*darrdZ)
    return darrdtheta, darrdpsi


def two_d_thita_psi_der(arr):#creates a 2-D matrix of d/dtheta and d/dpsi derivatives.
    global Zi,Ri,pot,jcut,bfield,B2
    Z_dim = arr.shape[0]
    
    darrdtheta, darrdpsi = zip(*[thita_psi_der(arr,cut) for cut in range(0,Z_dim)])
    ddarrdtheta = np.asarray(darrdtheta)
    ddarrdpsi = np.asarray(darrdpsi)#derivative[Z,R]
    return ddarrdtheta,ddarrdpsi


def Reynolds_force(arr):#takes the d/dtheta and d/dpsi of a Reynolds stress matrix.
    Re_pol, Re_rad = two_d_thita_psi_der(arr)
    return -Re_pol, -Re_rad 



def Re_call_n_save(r):
    Re_tps,Re_tt,Re_psps,Re_eq_tps,Re_eq_tt,Re_eq_psps,Re_tur_tps,Re_tur_tt,Re_tur_psps = Reynolds_stress()
    np.save('Reynolds_tp_%s' %(r),np.nan_to_num(Re_tps))
    np.save('Reynolds_tt_%s' %(r),np.nan_to_num(Re_tt))
    np.save('Reynolds_pp_%s' %(r),np.nan_to_num(Re_psps))
    np.save('Reynolds_eq_tp_%s' %(r),np.nan_to_num(Re_eq_tps))
    np.save('Reynolds_eq_tt_%s' %(r),np.nan_to_num(Re_eq_tt))
    np.save('Reynolds_eq_pp_%s' %(r),np.nan_to_num(Re_eq_psps))
    np.save('Reynolds_tur_tp_%s' %(r),np.nan_to_num(Re_tur_tps))
    np.save('Reynolds_tur_tt_%s' %(r),np.nan_to_num(Re_tur_tt))
    np.save('Reynolds_tur_pp_%s' %(r),np.nan_to_num(Re_tur_psps))
#total force
    rfptp,rfrpt = Reynolds_force(np.nan_to_num(Re_tps))
    rfptt,rfrtt = Reynolds_force(np.nan_to_num(Re_tt))
    rfppp,rfrpp = Reynolds_force(np.nan_to_num(Re_psps))

    np.save('Re_pol_f_%s' %(r), rfrpt+rfptt)
    np.save('Re_rad_f_%s' %(r), rfptp+rfrpp)
#turbulent force    
    rftptp,rftrpt = Reynolds_force(np.nan_to_num(Re_tur_tps))
    rftptt,rftrtt = Reynolds_force(np.nan_to_num(Re_tur_tt))
    rftppp,rftrpp = Reynolds_force(np.nan_to_num(Re_tur_psps))

    np.save('Re_tur_pol_f_%s' %(r), rftrpt+rftptt)
    np.save('Re_tur_rad_f_%s' %(r), rftptp+rftrpp)












