'''Module that calculates Normal-to-X-line-fluxes.'''
import numpy as np
import xgc
import math
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
import sys
from decimal import Decimal
import core
import angle

phi_start=0
phi_end=None

fileDir='/global/cscratch1/sd/giannos/ti255_d3d_Ip_1.5_med'

#Parameters for ti255-Case, R-values for LFS leg
#Rmin = 1.4648
#Rmax = 1.8
#Zmin = -1.19
#Zmax = -1.1

#Higher X-point.
#Rmin = 1.4648
#Rmax = 1.99
#Zmin = -1.15
#Zmax = -0.83

#Higher X-point for both legs.
#Rmin = 1.001
#Rmax = 1.99
#Zmin = -1.15
#Zmax = -0.83


#Parameters for ti255-Case, R-values for HFS leg
#Rmax = 1.453
#Rmin = 1.001
#Zmin = -1.19
#Zmax = -1.1

#Parameters for ti-255-Case, R-values for both legs
#Rmin = 1.001
#Rmax = 1.8
#Zmin = -1.19
#Zmax = -1.1

'''Parameter Settings for beggining and ending of turbulence part.'''

#Top point
Rmin = 1.52
Rmax = 2.0
Zmin = 0.89
Zmax = 0.95

#Bottom point
#Rmin = 1.79
#Rmax = 2.0
#Zmin = -0.82
#Zmax = -0.74

#change number according to the accuracy. Double the accuracy when you run both legs.
accuracy = 100
unitR = (Rmax-Rmin)/accuracy
unitZ = (Zmax-Zmin)/accuracy
core.getMeshAndCuts(fileDir,Rmin,Rmax,Zmin,Zmax)

'''Functions that write files with all the normal-to-horizontal plane fluxes, below the X-point. We assume that all these fluxes eventually    end up at the divertors.
''' 

def Particle_N_flux():
    '''Writes on text the Z-component (normal to a horizontal line) of the particle fluxes (exb and parallel components).'''
    option = 1
    start = 0
    time_range = 200
    cut = accuracy/2
    
    #new_file = open("Norm_Particle_Fluxes_higher@%s.txt" %(cut),'a')
    new_file = open("Norm_Particle_Fluxes_Top.txt",'a')
    new_file.write("ExB"+"\t"+"Par"+"\t"+"R"+"\t"+"psi"+"\n")
    br = core.getcutvalue_hor(core.bfield[:,0],cut,1)
    
    #loading arrays
    temp = np.array([[core.getcutvRvZ(iz,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    VZ_all = temp[:,:,1,:]
    print("VZ_all done...")
    
    ne_all = np.array([[core.getcutvalue_hor(core.ne[:,iz,it],cut,option) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("ne_all done...")
    
    Vpar_all = np.array([[core.map_array(core.e_u_para,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Vpar_all done...")
    
    for Flux_ip in range(0,accuracy):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else:
            exb, par = core.particle_fluxes_norm(cut,Flux_ip,ne_all,VZ_all,Vpar_all)
            R = (Rmin + unitR*Flux_ip)
            psi = core.psii[Flux_ip]
            new_file.write(str(exb)+"\t"+str(par)+"\t"+str(R)+"\t"+str(psi)+"\n")
    new_file.close()
    return("Done!")

#Particle_N_flux()

def X_den():
    '''Writes on text the density on a line at the X-point.'''
    
    cut = accuracy/2
    new_file = open("density_LFS@%s.txt" %(cut),'a')
    new_file.write("den"+"\t"+"R"+"\t"+"psi"+"\n")
    br = core.getcutvalue_hor(core.bfield[:,0],cut,1)
    for Flux_ip in range(0,accuracy):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else:
            _, _, _, _,n = core.norm_parallel_flux(cut,Flux_ip)
            R = (Rmin + unitR*Flux_ip)
            psi = core.psii[Flux_ip]
            new_file.write(str(n)+"\t"+str(R)+"\t"+str(psi)+"\n")
    new_file.close()
    return("Done!")

#X_den()

def X_Vpar():
    '''Writes on text the density on a line at the X-point.'''
    cut = accuracy/2
    new_file = open("Vpar_LFS@%s.txt" %(cut),'a')
    new_file.write("Vpar_e"+"\t"+"Vpar_i"+"\t"+"R"+"\t"+"psi"+"\n")
    br = core.getcutvalue_hor(core.bfield[:,0],cut,1)
    for Flux_ip in range(0,accuracy):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else:
            _, _, _, _,_,Ve,Vi = core.norm_parallel_flux(cut,Flux_ip)
            R = (Rmin + unitR*Flux_ip)
            psi = core.psii[Flux_ip]
            new_file.write(str(Ve)+"\t"+str(Vi)+"\t"+str(R)+"\t"+str(psi)+"\n")
    new_file.close()
    return("Done!")

#X_Vpar()

def Mag_N_flux():
    '''Writes on text the Z-component of the magnetic drifts flux.'''
    option = 1
    start = 0
    time_range = 200
    cut = accuracy/2
    
    #new_file = open("Norm_Magnetic_Fluxes_higher@%s.txt" %(cut),'a')
    new_file = open("Norm_Magnetic_Fluxes_Bottom.txt",'a')
    new_file.write("Magnetic"+"\t"+"R"+"\t"+"psi"+"\n")
    br = core.getcutvalue_hor(core.bfield[:,0],cut,1)
    
    ne_all = np.array([[core.getcutvalue_hor(core.ne[:,iz,it],cut,option) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("ne_all done...")
    
    for Flux_ip in range(0,accuracy):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else:
            Mag = core.magnetic_flux_norm(cut,Flux_ip,ne_all)
            R = (Rmin + unitR*Flux_ip)
            psi = core.psii[Flux_ip]
            new_file.write(str(Mag)+"\t"+str(R)+"\t"+str(psi)+"\n")
    new_file.close()
    return("Done!")

#Mag_N_flux()

def e_heat_N_flux():
    '''Writes on text the Z-component of the exb and parallel electron heat flux.'''
    option = 1
    start = 0
    time_range = 200
    cut = accuracy/2
    
    #new_file = open("Norm_e_heat_Fluxes@%s.txt" %(cut),'a')
    new_file = open("Norm_e_heat_Fluxes_Top.txt",'a')
    new_file.write("Eq_exb"+"\t"+"Tur_exb"+"\t"+"Eq_par"+"\t"+"Tur_par"+"\t"+"R"+"\t"+"psi"+"\n")
    br = core.getcutvalue_hor(core.bfield[:,0],cut,1)
    
    temp = np.array([[core.getcutvRvZ(iz,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    VZ_all = temp[:,:,1,:]
    print("VZ_all done...")
    ne_all = np.array([[core.getcutvalue_hor(core.ne[:,iz,it],cut,option) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("ne_all done...")
    Te_perp_all = np.array([[core.map_array(core.e_T_perp,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Te_perp_all done...")
    Ee_para_all = np.array([[core.map_array(core.e_E_para,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Ee_para_all done...")
    Vpar_all = np.array([[core.map_array(core.e_u_para,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Vpar_all done...")
    
    for Flux_ip in range(0,accuracy):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else:
            exb, par = core.e_heat_flux_norm(cut,Flux_ip,VZ_all,ne_all,Te_perp_all,Ee_para_all,Vpar_all)
            R = (Rmin + unitR*Flux_ip)
            psi = core.psii[Flux_ip]
            new_file.write(str(exb)+"\t"+str(par)+"\t"+str(R)+"\t"+str(psi)+"\n")
    new_file.close()
    return("Done!")

e_heat_N_flux()

def i_heat_N_flux():
    '''Writes on text the Z-component of the exb and parallel ion heat flux.'''
    option = 1
    start = 0
    time_range = 200
    cut = accuracy/2
    
    #new_file = open("Norm_i_heat_Fluxes@%s.txt" %(cut),'a')
    new_file = open("Norm_i_heat_Fluxes_Bottom.txt",'a')
    new_file.write("Eq_exb"+"\t"+"Tur_exb"+"\t"+"Eq_par"+"\t"+"Tur_par"+"\t"+"R"+"\t"+"psi"+"\n")
    br = core.getcutvalue_hor(core.bfield[:,0],cut,1)
    
    temp = np.array([[core.getcutvRvZ(iz,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    VZ_all = temp[:,:,1,:]
    print("VZ_all done...")
    ne_all = np.array([[core.getcutvalue_hor(core.ne[:,iz,it],cut,option) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("ne_all done...")
    Ti_perp_all = np.array([[core.map_array(core.i_T_perp,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Ti_perp_all done...")
    Ei_para_all = np.array([[core.map_array(core.i_E_para,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Ei_para_all done...")
    Vpar_all = np.array([[core.map_array(core.e_u_para,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Vpar_all done...")
    
    for Flux_ip in range(0,accuracy):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else:
            exb, par = core.i_heat_flux_norm(cut,Flux_ip,VZ_all,ne_all,Ti_perp_all,Ei_para_all,Vpar_all)
            R = (Rmin + unitR*Flux_ip)
            psi = core.psii[Flux_ip]
            new_file.write(str(exb)+"\t"+str(par)+"\t"+str(R)+"\t"+str(psi)+"\n")
    new_file.close()
    return("Done!")

#i_heat_N_flux()

def e_mag_heat_N_flux():
    '''Writes on text the Z-component of the magnetic drifts electron heat flux.'''
    option = 1
    start = 0
    time_range = 200
    cut = accuracy/2
    
    #new_file = open("Norm_e_mag_heat_Fluxes@%s.txt" %(cut),'a')
    new_file = open("Norm_e_mag_heat_Fluxes_Top.txt",'a')
    new_file.write("Magnetic"+"\t"+"R"+"\t"+"psi"+"\n")
    br = core.getcutvalue_hor(core.bfield[:,0],cut,1)
    
    ne_all = np.array([[core.getcutvalue_hor(core.ne[:,iz,it],cut,option) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("ne_all done...")
    Te_perp_all = np.array([[core.map_array(core.e_T_perp,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Te_perp_all done...")
    Ee_para_all = np.array([[core.map_array(core.e_E_para,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Ee_para_all done...")
    
    for Flux_ip in range(0,accuracy):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else:
            mag = core.mag_e_heat_flux_norm(cut,Flux_ip,ne_all,Te_perp_all,Ee_para_all)
            R = (Rmin + unitR*Flux_ip)
            psi = core.psii[Flux_ip]
            new_file.write(str(mag)+"\t"+str(R)+"\t"+str(psi)+"\n")
    new_file.close()
    return("Done!")

e_mag_heat_N_flux()

def i_mag_heat_N_flux():
    '''Writes on text the Z-component of the magnetic drifts ion heat flux.'''
    option = 1
    start = 0
    time_range = 200
    cut = accuracy/2
    
    #new_file = open("Norm_i_mag_heat_Fluxes@%s.txt" %(cut),'a')
    new_file = open("Norm_i_mag_heat_Fluxes_Bottom.txt",'a')
    new_file.write("Magnetic"+"\t"+"R"+"\t"+"psi"+"\n")
    br = core.getcutvalue_hor(core.bfield[:,0],cut,1)
    
    ne_all = np.array([[core.getcutvalue_hor(core.ne[:,iz,it],cut,option) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("ne_all done...")
    Ti_perp_all = np.array([[core.map_array(core.i_T_perp,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Ti_perp_all done...")
    Ei_para_all = np.array([[core.map_array(core.i_E_para,(iz-1)%core.Nplanes,it,cut) for iz in range(core.Nplanes)] for it in range(start,start+time_range)])
    print("Ei_para_all done...")
    
    for Flux_ip in range(0,accuracy):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else:
            mag = core.mag_i_heat_flux_norm(cut,Flux_ip,ne_all,Ti_perp_all,Ei_para_all)
            R = (Rmin + unitR*Flux_ip)
            psi = core.psii[Flux_ip]
            new_file.write(str(mag)+"\t"+str(R)+"\t"+str(psi)+"\n")
    new_file.close()
    return("Done!")

#i_mag_heat_N_flux()

def Bottom_Integral(Rlist,weight):
    '''Performs the surface integral on a plane at the bottom of the Tokamak.'''
    R_start = 0 #index of starting and ending points. Assume that all four arrays have the same size
    R_end = len(Rlist)-1
    prefactor = (Rlist[R_end]-Rlist[R_start])/(len(Rlist[R_start:R_end])-1)
    #Integration scheme
    weight = np.asarray(weight)
    Rlist = np.asarray(Rlist)
    temp = weight*Rlist
    Summation = (1/2)*weight[R_start]*Rlist[R_start] + (1/2)*weight[R_end]*Rlist[R_end] + np.sum(temp[R_start+1:R_end-2])
    Integral = 2*math.pi*prefactor*Summation
    #print("Integral = ",Integral)
    return Integral

def normalization(Rlist,weight):
    '''#Performs the Integral for the normalizing factor. So far, does only the R-integral and assumes that the integral in the toroidal          dimension extends from 0 to 2pi.
    '''
    weight2 = np.ones(len(weight))
    return Bottom_Integral(Rlist,weight2)

def file_read():
    '''Reads the text files with the results.'''
    norm_file = open("Norm_Particle_Fluxes_higher_LFS@50.txt",'r')
    #norm_file = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/Normal_to_X_line/HFS_divertor/Norm_Particle_Fluxes_HFS@50.txt",'r')
    next(norm_file)
    Eq = []
    Tur = []
    Eqp = []
    Turp = []
    R = []
    for columns in (raw.strip().split() for raw in norm_file):  
        Eq.append(float(columns[0]))
        Tur.append(float(columns[1]))
        Eqp.append(float(columns[2]))
        Turp.append(float(columns[3]))
        R.append(float(columns[4]))
    return Eq, Tur, Eqp, Turp, R

#Eq, Tur, Eqp, Turp, R = file_read()
#Equilibrium = Bottom_Integral(R[:],Eq[:])
#Turbulent = Bottom_Integral(R[:],Tur[:])
#eqp = Bottom_Integral(R[:],Eqp[:])
#turp = Bottom_Integral(R[:],Turp[:])
#print("Normal Eq exb Flux = ", Equilibrium)
#print("Normal Tur exb Flux = ", Turbulent)
#print("eqp = ", eqp)
#print("turp = ", turp)
def lineplot():
    plt.title("Par. Flux at the LFS Bottom")
    plt.xlabel('R(m)')
    plt.ylabel('Flux')
    plt.plot(R,Eqp)
    plt.show()

#lineplot()



