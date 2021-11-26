#!/usr/bin/env python
from mpi4py import MPI
import argparse
import numpy as np
import math
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
import sys
from decimal import Decimal
import xgc
import core
import angle
import Reyn

'''Module that contains the functions that break up the loading of patches to different processes. The different patches cover the entire
   the entire separatrix. It also contains functions that do the calculation of various quantities on separatrix points and store the data
   in text files. Things that need to be checked before different runs:
   1) Correct fileDir for reading in simulation data
   2) Correct separatrix file at sep_file_read()
   3) Coordinates to avoid X-point at loading_checks()
   4) Proper naming of output text files with the correct numbers at all of the functions
   5) It needs to be called with appropriate number of processess, allowed time and correct function name in the parser
'''

#initializing the MPI communicator and the processes rank.
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
r=rank
#r=0
#Creates the parser that reads function name to be called
parser = argparse.ArgumentParser(description='Runs parallely a calculation on patches that cover the entire separatrix.')
parser.add_argument('function', help='function name to be called. Options to use: check, exb_flux, magnetic_flux, diamagnetic_flux, electron_heat_flux, ion_heat_flux, mag, mag_component, dphi_strength, phis, dn_strength, sl, temp, temp_i, den, reynolds, shear_calc, dTe_strength, dTi_strength, dPe_strength, dPi_strength.')
args = parser.parse_args()

phi_start=0
phi_end=None

fileDir = '/global/cscratch1/sd/giannos/ti262_cmod_JTM'
#fileDir = '/global/cscratch1/sd/giannos/ti255_d3d_Ip_1.5_med'


def sep_file_read():
    '''Reads the separatrix node locations. Returns lists with R,Z coordinates of separatrix nodes.'''
    #sep_file = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/Sep-ti255_med.txt","r")
    #sep_file = open("IFSNL.txt","r")
    sep_file = open("/global/homes/g/giannos/xgc_python_dir/ti262_analysis/ti262_Separatrix.txt","r")
    next(sep_file)
    R=[]
    Z=[]
    for columns in (raw.strip().split() for raw in sep_file):
        R.append(columns[0])
        Z.append(columns[1])
    return R,Z

def loading_checks(Rmin,Rmax,Zmin,Zmax):
    '''Checks for reversal of min and max. Returns the R,Z coordinates of the patch to be loaded.'''
    temp1=0
    temp2=0
    if Rmin>=Rmax:
        temp1 = Rmin
        Rmin = Rmax
        Rmax = temp1
    else:
        pass

    if Zmin>=Zmax:
        temp2 = Zmin
        Zmin = Zmax
        Zmax = temp2
    else:
        pass

    if abs(Rmin-Rmax)<0.12: #avoid loading too narrow windows
        Rmin = Rmin-0.06
        Rmax = Rmax+0.06

    if abs(Zmin-Zmax)<0.16:
        Zmin = Zmin-0.08
        Zmax = Zmax+0.08
        #if Zmin-0.08<-1.15: #avoiding the X-point for ti255
            #Zmin = -1.15
        if Zmin-0.08<-0.39: #avoiding the X-point for ti262
            Zmin = -0.39
        else:
            Zmin = Zmin-0.08 
    return Rmin,Rmax,Zmin,Zmax

def check():    
    '''Simple check to see if the values are stored correctly. Useful to see (with verbose option of angle.spline_creation set to true),
    if the angular grid covers faithfully the separatrix.
    '''
    angle.angle_mesh()
    file = open("Check_%s.txt" %(r),'a')
    for i in range(1, len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        file.write(str(Ang)+"\t"+str(R_val)+"\t"+str(Z_val)+"\n")
    file.close()
    
'''Next four functions calculate all the fluxes vs. poloidal angle.'''
    
def exb_flux():
    '''Stores the particle, exb, equilibrium and turbulent fluxes vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti255_particle_fluxes_vs.angle_%s.txt" %(r), 'a')
    file.write("Angle"+"\t"+"Eq.Flux"+"\t"+"Tur.Flux"+"\t"+"Ring"+"\t"+"R"+"\t"+"Z"+"\n")
    
    #calls flux function on diagonal points of the mesh
    for i in range(1, len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Eq, Tur, Ring = angle.exb_fluxes_angle(i,i,False)
        file.write(str(Ang)+"\t"+str(Eq)+"\t"+str(Tur)+"\t"+str(Ring)+"\t"+str(R_val)+"\t"+str(Z_val)+"\n")
    file.close()

def magnetic_flux():
    '''Stores ion magnetic drift fluxes vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_mag_fluxes_vs.angle_%s.txt" %(r),'a')
    file.write("Angle"+"\t"+"Mag.Flux"+"\t"+"R"+"\t"+"Z"+"\n")
    #calls flux function on diagonal points of the mesh
    for i in range(1, len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Mag  = angle.magnetic_flux_angle(i,i)
        file.write(str(Ang)+"\t"+str(Mag)+"\t"+str(R_val)+"\t"+str(Z_val)+"\n")
    file.close()    

def diamagnetic_flux():
    '''Stores ion diamagnetic drift fluxes vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_diamag_fluxes_vs.angle_%s.txt" %(r),'a')
    file.write("Angle"+"\t"+"D.Flux"+"\t"+"R"+"\t"+"Z"+"\n")
    #calls flux function on diagonal points of the mesh
    for i in range(1, len(angle.R_fp)-2):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        dMag  = angle.diamagnetic_drifts_angle(i,i)
        file.write(str(Ang)+"\t"+str(dMag)+"\t"+str(R_val)+"\t"+str(Z_val)+"\n")
    file.close()    
    
def electron_heat_flux():
    '''Stores electron exb heat fluxes vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti255_electron_heat_%s.txt" %(r),'a')
    file.write("Angle"+"\t"+"Eq.Flux"+"\t"+"Tur.Flux"+"\t"+"R"+"\t"+"Z"+"\n")
    for i in range(1, len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Eq, Tur = angle.electron_heat_flux_angle(i,i,True)
        file.write(str(Ang)+"\t"+str(Eq)+"\t"+str(Tur)+"\t"+str(R_val)+"\t"+str(Z_val)+"\n")
        print("Flux point %s done." %(i))
    file.close()

def ion_heat_flux():
    '''Stores ion exb heat fluxes vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti255_ion_heat_%s.txt" %(r),'a')
    file.write("Angle"+"\t"+"Eq.Flux"+"\t"+"Tur.Flux"+"\t"+"R"+"\t"+"Z"+"\n")
    for i in range(1, len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Eq, Tur = angle.ion_heat_flux_angle(i,i,True)
        file.write(str(Ang)+"\t"+str(Eq)+"\t"+str(Tur)+"\t"+str(R_val)+"\t"+str(Z_val)+"\n")
        print("Flux point %s done." %(i))
    file.close()

'''Following functions calculate useful to have quantities vs. poloidal angle.'''
    
def mag(): 
    '''Stores magnitude of total and poloidal magnetic field vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_magnetic_field_vs.angle_%s.txt" %(r), 'a')
    file.write("Angle"+"\t"+"Bp"+"\t"+"B"+"\n")
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Bp, B = angle.magnetic_field(i,i)
        file.write(str(Ang)+"\t"+str(Bp)+"\t"+str(B)+"\n")
    file.close()
    
def mag_component():
    '''Stores components of magnetic field vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_magnetic_field_comp_vs.angle_%s.txt" %(r), 'a')
    file.write("Angle"+"\t"+"Br"+"\t"+"Bz"+"\t"+"Bzeta"+"\n")
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Br, Bz, Bzeta = angle.mag_comp(i,i)
        file.write(str(Ang)+"\t"+str(Br)+"\t"+str(Bz)+"\t"+str(Bzeta)+"\n")
    file.close()    
    
def dphi_strength():
    '''Stores strength of electrostatic fluctuation vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_turbulence_strength_vs.angle_%s.txt" %(r), 'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        _, _, tur = angle.turbulence_strength(i,i)
        file.write(str(Ang)+"\t"+str(tur)+"\n")
    file.close()

def phis():
    '''Stores values potential and electrostatic fluctuation vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_phis_vs.angle_%s.txt" %(r), 'a')
    file.write("Angle"+"\t"+"Phi"+"\t"+"dPhi"+"\n")
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        phi, dphi, _ = angle.turbulence_strength(i,i)
        file.write(str(Ang)+"\t"+str(phi)+"\t"+str(dphi)+"\n")
    file.close()
    
def dn_strength():
    '''Stores strength of density perturbation vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_density_perturbation_strength_vs.angle_%s.txt" %(r), 'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        dn = angle.density_perturbation_strength(i,i)
        file.write(str(Ang)+"\t"+str(dn)+"\n")
    file.close()
    
def sl():
    '''Stores inverse scale lengths vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_Scale_lengths_%s.txt" %(r), 'a')
    #calls flux function on diagonal points of the mesh
    for i in range(2,len(angle.R_fp)-2):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        LTe = angle.LTeinv(i,i)
        Ln = angle.Lninv(i,i)
        LTi = angle.LTiinv(i,i)
        file.write(str(Ang)+"\t"+str(LTe)+"\t"+str(LTi)+"\t"+str(Ln)+"\t"+str(R_val)+"\t"+str(Z_val)+"\n")
    file.close()
    
def temp():
    '''Stores ion and electron total temperatures vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_temperatures_%s.txt" %(r), 'a')
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Ti, Te, _, _ = angle.tmp(i,i) 
        file.write(str(Ang)+"\t"+str(Ti)+"\t"+str(Te)+"\n")
    file.close()

def temp_i():
    '''Stores perpendicular and parallel ion temperatures vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti255_ion_temperatures_%s.txt" %(r), 'a')
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        _, _, Tperp, Ei = angle.tmp(i,i) 
        file.write(str(Ang)+"\t"+str(Tperp)+"\t"+str(2*Ei)+"\n")
    file.close()
    
def den():
    '''Stores density vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti262_densities_%s.txt" %(r), 'a')
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        n = angle.density(i,i) 
        file.write(str(Ang)+"\t"+str(n)+"\n")
    file.close()
    
def dn():
    angle.angle_mesh()
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        dn = angle.density_perturbation(i,i)
        for t in range(2):
            for p in range(core.Nplanes):
                file = open("ti255_dn_%s_t_%s_pl_%s.txt" %(r, t, p),'a')
                file.write(str(Ang)+"\t"+str(dn[t,p])+"\n")
                file.close()

def dTe_strength():
    '''Stores strength of electron temperature perturbation vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti255_Te_perturbation_strength_vs.angle_%s.txt" %(r), 'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        dTe = angle.e_temp_perturbation_strength(i,i)
        file.write(str(Ang)+"\t"+str(dTe)+"\n")
    file.close()
                
def dTi_strength():
    '''Stores strength of ion temperature perturbation vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti255_Ti_perturbation_strength_vs.angle_%s.txt" %(r), 'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        dTi = angle.i_temp_perturbation_strength(i,i)
        file.write(str(Ang)+"\t"+str(dTi)+"\n")
    file.close()

def dPi_strength():
    '''Stores strength of ion pressure perturbation vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti255_Pi_perturbation_strength_vs.angle_%s.txt" %(r), 'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        dPi = angle.i_pressure_perturbation_strength(i,i)
        file.write(str(Ang)+"\t"+str(dPi)+"\n")
    file.close()
    
def dPe_strength():
    '''Stores strength of electron pressure perturbation vs. poloidal angle.'''
    angle.angle_mesh()
    file = open("ti255_Pe_perturbation_strength_vs.angle_%s.txt" %(r), 'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        dPe = angle.e_pressure_perturbation_strength(i,i)
        file.write(str(Ang)+"\t"+str(dPe)+"\n")
    file.close()

'''FUNCTIONS FOR THE CALCULATION OF SHEAR. Important to make sure at which time step shear needs to be calculated.'''

def psi_der(arr,it,cut):#takes the d/dpsi on the potential and divides by RBp
    '''Takes a single derivative wrt psi along a cut.'''
    option=1
    if cut==None:
        cut = core.jcut
        
    dZ = core.Zi[1] - core.Zi[0]
    dR = core.Ri[1] - core.Ri[0]
    
    unitR = (core.loader.Rmax-core.loader.Rmin)/len(core.Ri)
    unitZ = (core.loader.Zmax-core.loader.Zmin)/len(core.Zi)
    
    arrim1 = np.asarray([core.getcutvalue_hor(arr[:,iz,it],(cut-1)%len(core.Zi),option) for iz in range(0,core.Nplanes)]).mean(axis=0)
    arri = np.asarray([core.getcutvalue_hor(arr[:,iz,it],cut,option) for iz in range(0,core.Nplanes)]).mean(axis=0)
    arrip1 = np.asarray([core.getcutvalue_hor(arr[:,iz,it],(cut+1)%len(core.Zi),option) for iz in range(0,core.Nplanes)]).mean(axis=0)
    darrdZ = (arrip1-arrim1)/(2.*dZ)
    
    arriplus = np.roll(arri,-1)
    arriminus = np.roll(arri,1)
    darrdR = (arriplus-arriminus)/(2.*dR)
    darrdR[0] = (arri[1]-arri[0])/dR
    darrdR[-1] = (arri[-1]-arri[-2])/dR
    
    BRi = core.getcutvalue_hor(core.bfield[:,0],cut,option)
    BZi = core.getcutvalue_hor(core.bfield[:,1],cut,option)
    Bzetai = core.getcutvalue_hor(core.bfield[:,2],cut,option)
    Bp = np.sqrt(np.square(BRi)+np.square(BZi))
    L = np.array([x for x in range(0,len(core.Ri))])
    R = np.array([L[ip]*unitR+core.loader.Rmin for ip in range(0,len(core.Ri))])+core.Rmaj

    
    darrdpsi = (1/Bp)*(BRi*darrdZ - BZi*darrdR)
    return -darrdpsi/(R*Bp)

def two_d_psi_der(it):#produces 2D E_psi/RBp matrix
    '''Produces a 2D matrix of psi derivatives.'''
    E_rad = [psi_der(core.pot,it,cut) for cut in range(0,len(core.Zi))]
    E_psi = np.asarray(E_rad)
    return E_psi

def shear(it,cut):#takes the d/dpsi of E_psi/RBp and multiplies with (RBp)^2/B to find the shear.
    '''Calculates the Burrell shear along a cut at a particular time.'''
    option = 1
    c = 3.0e8
    dZ = core.Zi[1] - core.Zi[0]
    dR = core.Ri[1] - core.Ri[0]
    
    unitR = (core.loader.Rmax-core.loader.Rmin)/len(core.Ri)        
    
    dphidpsi = two_d_psi_der(it)
    
    arrim1 = dphidpsi[(cut-1)%len(core.Zi),:]
    arri = dphidpsi[cut,:]
    arrip1 = dphidpsi[(cut+1)%len(core.Zi),:]
    darrdZ = (arrip1-arrim1)/(2.*dZ)
    
    arriplus = np.roll(arri,-1)
    arriminus = np.roll(arri,1)
    darrdR = (arriplus-arriminus)/(2.*dR)
    darrdR[0] = (arri[1]-arri[0])/dR
    darrdR[-1] = (arri[-1]-arri[-2])/dR
    
    BRi = core.getcutvalue_hor(core.bfield[:,0],cut,option)
    BZi = core.getcutvalue_hor(core.bfield[:,1],cut,option)
    Bzetai = core.getcutvalue_hor(core.bfield[:,2],cut,option)
    Bp = np.sqrt(np.square(BRi)+np.square(BZi))
    Bmag = np.sqrt(np.square(BRi)+np.square(BZi)+np.square(Bzetai))
    L = np.array([x for x in range(0,len(core.Ri))])
    R = np.array([L[ip]*unitR+core.loader.Rmin for ip in range(0,len(core.Ri))])+core.Rmaj

    doubledpsi = (1/Bp)*(BRi*darrdZ - BZi*darrdR)
    Shear = ((np.square(R)*np.square(Bp))/Bmag)*doubledpsi
    return Shear

def two_d_shear(it):
    '''Calculates the Burell shear on an RZ plane. Time point needs to be specified.'''
    Shear = [shear(it,cut) for cut in range(0,len(core.Zi))]
    Bur_shear = np.asarray(Shear)
    return Bur_shear

def shear_call_n_save(r):
    '''Function that saves the 2D shear values on an array. Gets called from parallel program.'''
    sh = two_d_shear(0)
    sh2 = np.nan_to_num(sh)
    np.save("shear_%s" %(r),sh2)

def sep_save(r):
    '''Function that saves the separatrix locations on arrays. Gets called from parallel program.''' 
    Rpol = core.RZ[core.sepInds[:],0]
    Zpol = core.RZ[core.sepInds[:],1]
    np.save("Rsep_%s" %(r),Rpol)
    np.save("Zsep_%s" %(r),Zpol)
  
    
'''Those functions are called directly from other modules.'''    
def reynolds():
    '''Calculates all components of Reynolds stress tensor. Reynolds stresses.'''
    Reyn.Re_call_n_save(r)

'''Here, the function shear produces only 2D arrays of the Burell shear rate. Those files need to be postprocesed with postpro_shear to produce shear values on the separatrix. Separatrix arrays are also needed. 
'''  
def shear_calc():
    '''Calculates Burrel Shear rate and stores values on text files.'''
    shear_call_n_save(r)
    sep_save(r)
        

'''Basic loading scheme proceeds in intervals of 20 separatrix node points at a time.'''
R,Z = sep_file_read()
if ((r*20)<(len(R)-1))and (((r+1)*20)<(len(R)-1)):
    Rmin = float(R[r*20])#2.14#
    Rmax = float(R[(r+1)*20])#2.26#
    Zmin = float(Z[r*20])#-0.38#
    Zmax = float(Z[(r+1)*20])#-0.45#
    Rmin,Rmax,Zmin,Zmax = loading_checks(Rmin,Rmax,Zmin,Zmax)
    print(Rmin,Rmax,Zmin,Zmax)
    core.getMeshAndCuts(fileDir,Rmin,Rmax,Zmin,Zmax)
    #runs the function that was given as an input from the command line.
    locals()[args.function]()
    
elif ((r*20)<(len(R)-1))and (((r+1)*20)>(len(R)-1)):
    Rmin = float(R[r*20])
    Rmax = float(R[-1])
    Zmin = float(Z[r*20])
    Zmax = float(Z[-1])

    Rmin,Rmax,Zmin,Zmax = loading_checks(Rmin,Rmax,Zmin,Zmax)
    core.getMeshAndCuts(fileDir,Rmin,Rmax,Zmin,Zmax)
    locals()[args.function]()
    
else:
    print("Process exited without results.") 
 
