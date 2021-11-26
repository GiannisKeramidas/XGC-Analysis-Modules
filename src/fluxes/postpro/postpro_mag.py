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
    for i in range(0,77): #ti253-IFS(0,61)/ti255-separatrix (0,60)
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_mag_fluxes_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_diamag_fluxes_vs.angle_%d.txt" %(i),"r"))
        crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti262_mag_fluxes_vs.angle_%d.txt" %(i),"r"))
    for i in range(0,77):
        next(crs[i])
      
    angle = []
    mag_fl = []
    R_value = []
    Z_value = []
    for i in range(0,77):
        for columns in ( raw.strip().split() for raw in crs[i] ):  
            angle.append(float(columns[0]))
            mag_fl.append(float(columns[1]))
            R_value.append(float(columns[2]))
            Z_value.append(float(columns[3]))
   
    #filter out nan values from fluxes
    excl_val_2 = []
    for elem in np.where(np.isnan(mag_fl))[0][:]:
        excl_val_2.append(elem)
    for elem in sorted(excl_val_2,reverse=True):
        del angle[elem]
        del mag_fl[elem]
        del R_value[elem]
        del Z_value[elem]
   
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
                del mag_fl[elem]
                del R_value[elem]
                del Z_value[elem]
    
    #start all values from zero angle.
    nul_pos = angle.index(min(angle))
    angle = np.roll(angle,-nul_pos)
    mag_fl = np.roll(mag_fl,-nul_pos)
    R_value = np.roll(R_value,-nul_pos)
    Z_value = np.roll(Z_value,-nul_pos)
    print(angle)
    
    #returning the R-distance for the integral function.
    Rnorm = []
    for r,z in zip(R_value,Z_value):
        Rnorm.append(r**2 + z**2)

    R = R_value + core.Rmaj
    
    #Plotting option
    if verbose == True:
        plt.xlabel("angle")
        plt.ylabel("fluxes")
        plt.plot(angle[:], mag_fl[:],'b-')
        plt.show()

    #new_file = open("DiaMag_flux_20_smoothes.txt",'a')
    #new_file.write("angle"+"\t"+"DiaMagnetic"+"\t"+"R"+"\t"+"Z"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(mag_fl[i])+"\t"+str(R_value[i])+"\t"+str(Z_value[i])+"\n")
    #new_file.close()
    
    new_file = open("ti262_Mag_flux.txt",'a')
    new_file.write("angle"+"\t"+"Magnetic"+"\t"+"R"+"\t"+"Z"+"\n")
    for i in range(0,len(angle)):
        new_file.write(str(angle[i])+"\t"+str(mag_fl[i])+"\t"+str(R_value[i])+"\t"+str(Z_value[i])+"\n")
    new_file.close()
    
    return angle, mag_fl, Rnorm, R

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

def stagnation_point_lists():
    #stagnation point integration starts from 259 degrees and ends at 6 degrees. Therefore, we roll the matrices and return only the relevant values.
    angle, mag_fl, r, R = file_read(False)
    nul_pos = 102
    angle = np.roll(angle,nul_pos)
    mag_fl = np.roll(mag_fl,nul_pos)
    r = np.roll(r,nul_pos)
    R = np.roll(R,nul_pos)
    return angle[0:108], mag_fl[0:108],r[0:108],R[0:108]



#angle, mag_fl,r, R = file_read(False)
#angle, mag_fl, r, R = stagnation_point_lists()
#new_weight = np.ones(len(angle))
#max_weight = np.max(mag_fl)
#Magnetic = Tokamak_Integral(R,r,angle,mag_fl)
#Normalization = Tokamak_Integral(R,r,angle,new_weight)
#print("Magnetic Integrated Flux = ",Magnetic)
#print("Flux Surface Averaged Flux = ",Magnetic/Normalization)
#print("FSA Flux/Maximum Integrand = ", Magnetic/(Normalization*max_weight))
