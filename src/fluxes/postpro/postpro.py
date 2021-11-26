import numpy as np
import xgc
import math
import pylab as py
import matplotlib.pyplot as plt
import sys
import numpy.ma as ma
import core
from decimal import Decimal
import angle

def file_read(verbose = False):
    '''Reads all the files and joins the values in lists. Then, it filters duplicates and nans.'''
    tolerance = 0.5
    crs = []
    for i in range(0,77): #ti253-IFS(0,61)/ti255-separatrix (0,60)/ti262 (0,78)
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_particle_fluxes_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti262_magnetic_field_vs.angle_%d.txt" %(i),"r"))
        crs.append(open("/global/homes/g/giannos/xgc_python_dir/Check_%d.txt" %(i),"r"))
    #for i in range(0,77):
        #next(crs[i])
    
    angle = []
    R = []
    Z = []
    #eq_fl = []
    #tur_fl = []
    #ring = []
    #R_value = []
    #Z_value = []
    
    for i in range(0,77):
        for columns in ( raw.strip().split() for raw in crs[i] ):
            angle.append(float(columns[0]))
            R.append(float(columns[1]))
            Z.append(float(columns[2]))
            #eq_fl.append(float(columns[1]))
            #tur_fl.append(float(columns[2]))
            #ring.append(float(columns[3]))
            #R_value.append(float(columns[4]))
            #Z_value.append(float(columns[5]))
    #filter out nan values from fluxes.
    excl_val_2 = []
    for elem in np.where(np.isnan(R))[0][:]:
        excl_val_2.append(elem)
    for elem in sorted(excl_val_2,reverse=True):
        del angle[elem]
        del R[elem]
        del Z[elem]
        #del eq_fl[elem]
        #del tur_fl[elem]
        #del ring[elem]
        #del R_value[elem]
        #del Z_value[elem]
   
    excl_val_3 = []
    for elem in np.where(np.isnan(Z))[0][:]:
        excl_val_3.append(elem)
    for elem in sorted(excl_val_3,reverse=True):
        del angle[elem]
        del R[elem]
        del Z[elem]
        #del eq_fl[elem]
        #del tur_fl[elem]
        #del ring[elem]
        #del R_value[elem]
        #del Z_value[elem]

    #excl_val_4 = []
    #for elem in np.where(np.isnan(ring))[0][:]:
        #excl_val_4.append(elem)
    #for elem in sorted(excl_val_4,reverse=True):
        #del angle[elem]
        #del eq_fl[elem]
        #del tur_fl[elem]
        #del ring[elem]
        #del R_value[elem]
        #del Z_value[elem]
         
    #Keeping only one entry per angle, filtering out all others.
    B = []
    for i in range(0,len(angle)):
        B.append(int(angle[i]))
    
    for i in range(0,len(B)-1):
        excl_val = []
        for j in range(i+1,len(B)):
            diff = abs(B[i]-B[j])
            if (diff<tolerance) or (diff>360-tolerance):
                excl_val.append(j)
        if len(excl_val)>0:
            for elem in sorted(excl_val,reverse=True):
                del B[elem]
                del angle[elem]
                del R[elem]
                del Z[elem]
                #del eq_fl[elem]
                #del tur_fl[elem]
                #del ring[elem]
                #del R_value[elem]
                #del Z_value[elem]
    #Start all values from zero angle    
    nul_pos = angle.index(min(angle))
    angle = np.roll(angle,-nul_pos)
    R = np.roll(R,-nul_pos)
    Z = np.roll(Z,-nul_pos)
    #eq_fl = np.roll(eq_fl,-nul_pos)
    #tur_fl = np.roll(tur_fl,-nul_pos)
    #ring = np.roll(ring,-nul_pos)
    #R_value = np.roll(R_value,-nul_pos)
    #Z_value = np.roll(Z_value,-nul_pos)
    print(angle)
    
    #Returning the R-distance for the integral function. This is the r-distance from the center of the tokamak to the flux point.
    #Rnorm = [] 
    #for r,z in zip(R_value,Z_value):
        #Rnorm.append(r**2 + z**2)
    
    #Returning the distance from the axis of the tokamak for the R dzeta part of the integral.
    #R = R_value + core.Rmaj

    #Plotting option.
    if verbose == True:
        plt.title("Electron Heat Fluxes")
        plt.xlabel("angle")
        plt.ylabel("fluxes")
        plt.xlim([0.0,360.0])
        plt.plot(angle[:], eq_fl[:],'b-', angle[:], tur_fl[:],'r-')
        plt.show()

    #new_file = open("particle_fluxes_only_eq.txt",'a')
    #new_file.write("angle"+"\t"+"Equilibrium"+"\t"+"Turbulent"+"\t"+"Ring"+"\t"+"R"+"\t"+"Z"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(eq_fl[i])+"\t"+str(tur_fl[i])+"\t"+str(ring[i])+"\t"+str(R_value[i])+
                     #"\t"+str(Z_value[i])+"\n")
    #new_file.close()

    
    new_file = open("ti262_R_Z_sep.txt",'a')
    new_file.write("angle"+"\t"+"R"+"\t"+"Z"+"\n")
    for i in range(0,len(angle)):
        new_file.write(str(angle[i])+"\t"+str(R[i])+"\t"+str(Z[i])+"\n")
    new_file.close()
    
    return angle#, eq_fl, tur_fl, Rnorm, R

file_read()

def Line_int_along_curve(rarray,thetalist,weight):#function that returns the line integral of a weight along a curve
    #index of starting and ending angles. Assume that all three arrays have the same size
    theta_start = 0
    theta_end = len(thetalist)-1
    prefactor = (thetalist[theta_end]-thetalist[theta_start])/(len(thetalist[theta_start:theta_end])-1)
    #dr/dtheta
    drdtheta = []
    for i in range(1,theta_end):
        drdtheta.append((rarray[i+1]-rarray[i-1])/(thetalist[i+1]-thetalist[i-1]))
    #linear correction for first and last element of derivative list
    dr2dtheta_start = (rarray[2]-2*rarray[1] + rarray[0])/(thetalist[2]-thetalist[0])**2
    dr2dtheta_end = (rarray[theta_end]-2*rarray[theta_end-1] + rarray[theta_end-2])/(thetalist[theta_end]-thetalist[theta_end-2])**2
    last_elem = drdtheta[-1] + dr2dtheta_end
    drdtheta.append(last_elem)
    first_elem = drdtheta[0] - dr2dtheta_start
    a = first_elem
    drdtheta_new = [a]+drdtheta
    der = np.asarray(drdtheta_new)
    #Integration scheme of the full square root 
    temp = weight*rarray*np.sqrt(1 + (1/np.power(rarray,2))*np.power(der,2))
    Summation = (1/2)*weight[theta_start]*rarray[theta_start]*np.sqrt(1+(1/rarray[theta_start]**2)*der[theta_start]**2) +\
    (1/2)*weight[theta_end]*rarray[theta_end]*np.sqrt(1+(1/rarray[theta_end]**2)*der[theta_end]**2) +\
    np.sum(temp[theta_start+1:theta_end-2])
    Integral = prefactor*Summation
    #print("Integral = ",Integral)
    return Integral

#while the previous function gives the line integral along a curve and is useful, this one does the surface integral of a tokamak, i.e., contains
#the R dzeta integration as well. The rarray still refers to the distance from the center of the plane, while the Rarray is the distance from the 
#axis of the tokamak. Also contains the 2pi factor from the dzeta integration.
def Tokamak_Integral(Rarray,rarray,thetalist,weight):
    #index of starting and ending angles. Assume that all four arrays have the same size
    theta_start = 0
    theta_end = len(thetalist)-1
    prefactor = (thetalist[theta_end]-thetalist[theta_start])/(len(thetalist[theta_start:theta_end])-1)
    #dr/dtheta
    drdtheta = []
    for i in range(1,theta_end):
        drdtheta.append((rarray[i+1]-rarray[i-1])/(thetalist[i+1]-thetalist[i-1]))
    #linear correction for first and last element of derivative list
    dr2dtheta_start = (rarray[2]-2*rarray[1] + rarray[0])/(thetalist[2]-thetalist[0])**2
    dr2dtheta_end = (rarray[theta_end]-2*rarray[theta_end-1] + rarray[theta_end-2])/(thetalist[theta_end]-thetalist[theta_end-2])**2
    last_elem = drdtheta[-1] + dr2dtheta_end
    drdtheta.append(last_elem)
    first_elem = drdtheta[0] - dr2dtheta_start
    a = first_elem
    drdtheta_new = [a]+drdtheta
    der = np.asarray(drdtheta_new)
    #Integration scheme of the full square root 
    temp = weight*rarray*Rarray*np.sqrt(1 + (1/np.power(rarray,2))*np.power(der,2))
    Summation = (1/2)*weight[theta_start]*rarray[theta_start]*Rarray[theta_start]*np.sqrt(1+(1/rarray[theta_start]**2)*der[theta_start]**2) +\
    (1/2)*weight[theta_end]*rarray[theta_end]*Rarray[theta_end]*np.sqrt(1+(1/rarray[theta_end]**2)*der[theta_end]**2) +\
    np.sum(temp[theta_start+1:theta_end-2])
    Integral = 2*math.pi*prefactor*Summation
    #print("Integral = ",Integral)
    return Integral

#Performs the surface integral on a plane at the bottom of the Tokamak.
def Bottom_Integral(Rlist,weight):
    #index of starting and ending points. Assume that all four arrays have the same size
    R_start = 0
    R_end = len(Rlist)-1
    prefactor = (Rlist[R_end]-Rlist[R_start])/(len(Rlist[R_start:R_end])-1)
    #Integration scheme 
    temp = weight*Rlist
    Summation = (1/2)*weight[R_start]*Rlist[R_start] + (1/2)*weight[R_end]*Rlist[R_end] + np.sum(temp[R_start+1:R_end-2])
    Integral = 2*math.pi*prefactor*Summation
    #print("Integral = ",Integral)
    return Integral





def running_integral(angle,eq_fl,tur_fl,Rnorm):
    Equilibrium=[]
    Turbulent=[]
    angle_stop=[]
    new_angle =[]
    new_eq_fl =[]
    new_tur_fl =[]
    new_Rnorm =[]
 
    for i in range(4,len(angle)):
        angle_stop.append(angle[i])
        new_angle = angle[0:i]
        new_eq_fl = eq_fl[0:i]
        new_tur_fl = tur_fl[0:i]
        new_Rnorm = Rnorm[0:i]
        Equilibrium.append(Line_int_along_curve(new_Rnorm,new_angle,new_eq_fl))
        Turbulent.append(Line_int_along_curve(new_Rnorm,new_angle,new_tur_fl))
        new_angle =[]
        new_eq_fl =[]
        new_tur_fl =[]
        new_Rnorm =[]
    plt.title("Running total of flux integral")
    plt.plot(angle_stop,Equilibrium,'b--',angle_stop,Turbulent,'r--')
    plt.xlabel("Angle")
    plt.ylabel("Integrated Flux")
    plt.show() 

def stagnation_point_lists():
    #stagnation point integration starts from 259 degrees and ends at 6 degrees. Therefore, we roll the matrices and return only the relevant values.
    angle, eq_fl, tur_fl, r, R = file_read(False)
    nul_pos = 102
    angle = np.roll(angle,nul_pos)
    eq_fl = np.roll(eq_fl,nul_pos)
    tur_fl = np.roll(tur_fl,nul_pos)
    r = np.roll(r,nul_pos)
    R = np.roll(R,nul_pos)
    return angle[0:108], eq_fl[0:108],tur_fl[0:108],r[0:108],R[0:108]

def fluctuation():
    angle, eq_fl, tur_fl, r, R = file_read(False)
    eq = np.asarray(eq_fl)
    tur = np.asarray(tur_fl)
    eq_m = np.mean(eq)
    tur_m = np.mean(tur)
    eq_std = np.std(eq)
    tur_std = np.std(tur)
    eq_fluct = eq_std/eq_m
    tur_fluct = tur_std/tur_m
    #eq_m = sum(eq_fl)/float(len(eq_fl))
    #tur_m = sum(tur_fl)/float(len(tur_fl))
    print('equilibrium fluctuation = ',eq_fluct)
    print('turbulent fluctuation = ',tur_fluct)
    

#fluctuation()
#angle, eq_fl, tur_fl, r, R = file_read(True)
#angle, eq_fl, tur_fl, r, R = stagnation_point_lists()
#running_integral(angle,eq_fl,tur_fl,R)
#new_weight = np.ones(len(angle))
#max_weight_eq = np.max(eq_fl)
#max_weight_tur = np.max(tur_fl)
#Equilibrium = Tokamak_Integral(R,r,angle,eq_fl)
#Turbulent = Tokamak_Integral(R,r,angle,tur_fl)
#Normalization = Tokamak_Integral(R,r,angle,new_weight)
#print("Equilibrium Integrated Flux = ",Equilibrium)
#print("Turbulent Integrated Flux = ",Turbulent)
#print("FSA Equilibrium = ",Equilibrium/Normalization)
#print("FSA Turbulent = ",Turbulent/Normalization)
#print("Equilibrium cancelation = ",Equilibrium/(Normalization*max_weight_eq))
#print("Turbulent cancelation= ",Turbulent/(Normalization*max_weight_tur))

