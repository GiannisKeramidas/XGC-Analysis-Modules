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

#reads all the files and joins the values in lists. Then, it filters duplicates and nans.
def file_read(verbose = True):
    tolerance = 0.5
    crs = []
    for i in range(0,60):#ti253-separatrix(0,60)/ti253-IFS(0,61)/ti255-separatrix (0,60)
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti253_analysis/ti253_sep_flux/flux_data_new/ti253_sep_The_fluxes_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti253_analysis/ti253_IFS_flux/ti253_IFS_The_fluxes_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/sep_flux/long_run/ti255_The_fluxes_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/sep_flux/short_run/ti255_sep_short_The_fluxes_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/sep_flux/shorter_run/ti255_sep_short_The_fluxes_vs.angle_%d.txt" %(i),"r"))
        crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_sep_short_All_fluxes_vs.angle_%d.txt" %(i),"r"))
    for i in range(0,60):
        next(crs[i])
    angle = []
    eq_fl = []
    mean_fl = []
    tur_fl = []
    ttf = []
    total = []
    check = []
    R_value = []
    Z_value = []
    for i in range(0,60):
        for columns in ( raw.strip().split() for raw in crs[i] ):  
            angle.append(float(columns[0]))
            eq_fl.append(float(columns[1]))
            mean_fl.append(float(columns[2]))
            tur_fl.append(float(columns[3]))
            ttf.append(float(columns[4]))
            total.append(float(columns[5]))
            check.append(float(columns[6]))
            R_value.append(float(columns[7]))
            Z_value.append(float(columns[8]))
   
    excl_val_2 = []
    for elem in np.where(np.isnan(eq_fl))[0][:]:#filter out nan values from fluxes
        excl_val_2.append(elem)
    for elem in sorted(excl_val_2,reverse=True):
        del angle[elem]
        del eq_fl[elem]
        del mean_fl[elem]
        del tur_fl[elem]
        del ttf[elem]
        del total[elem]
        del check[elem]
        del R_value[elem]
        del Z_value[elem]
   

    excl_val_3 = []
    for elem in np.where(np.isnan(tur_fl))[0][:]:
        excl_val_3.append(elem)
    for elem in sorted(excl_val_3,reverse=True):
        del angle[elem]
        del eq_fl[elem]
        del mean_fl[elem]
        del tur_fl[elem]
        del ttf[elem]
        del total[elem]
        del check[elem]
        del R_value[elem]
        del Z_value[elem]
    
    B = []#Keeping only one entry per angle, filtering out all others.
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
                del angle[elem]
                del eq_fl[elem]
                del mean_fl[elem]
                del tur_fl[elem]
                del ttf[elem]
                del total[elem]
                del check[elem]
                del R_value[elem]
                del Z_value[elem]
                del B[elem]
        
    nul_pos = angle.index(min(angle))#start all values from zero angle
    angle = np.roll(angle,-nul_pos)
    eq_fl = np.roll(eq_fl,-nul_pos)
    mean_fl = np.roll(eq_fl,-nul_pos)
    tur_fl = np.roll(tur_fl,-nul_pos)
    ttf_fl = np.roll(eq_fl,-nul_pos)
    total_fl = np.roll(eq_fl,-nul_pos)
    check_fl = np.roll(eq_fl,-nul_pos)
    R_value = np.roll(R_value,-nul_pos)
    Z_value = np.roll(Z_value,-nul_pos)
    #print(angle)
    
    Rnorm = []#returning the R-distance for the integral function. This is the r-distance from the center of the tokamak to the flux point. 
    for r,z in zip(R_value,Z_value):
        Rnorm.append(r**2 + z**2)
    #Returning the distance from the axis of the tokamak for the R dzeta part of the integral.
    R = R_value + core.Rmaj
    check2 = eq_fl - mean_fl - ttf
    if verbose == True:# Plotting option
        fig, ax =plt.subplots()
        plt.xlabel("angle")
        plt.ylabel("fluxes")
        ax.plot(angle[:], eq_fl[:],'b-',label = 'Equilibrium')
        ax.plot(angle[:], tur_fl[:],'r-',label = 'Turbulent')
        ax.plot(angle[:],mean_fl[:],'g-',label = 'Mean Field')
        ax.plot(angle[:], total[:],'k-',label = 'Total')
        legend = ax.legend()
        #plt.plot(angle[:], check[:],'b-',angle[:], tur_fl[:],'r-')
        plt.show()
        fig,ax = plt.subplots()
        ax.plot(angle[:],check[:],'b-',label='check')
        ax.plot(angle[:],tur_fl[:],'r-',label='turbulent')
        legend=ax.legend()
        plt.show()

        fig,ax = plt.subplots()
        ax.plot(angle[:],eq_fl[:],'b-',label='Equilibrium')
        ax.plot(angle[:],mean_fl[:],'r-',label='Mean Field')
        ax.plot(angle[:],mean_fl[:]-ttf[:],'g',label = 'check')
        legend=ax.legend()
        plt.show()

        new_file = open("All_new_fluxes_total.txt",'a')
        new_file.write("angle"+"\t"+"Equilibrium"+"\t"+"Mean Field"+"\t"+"Turbulent"+"\t"+"ttf"+"\t"+"Total"+"\t"+"check"+"\t"+"R"+"\t"+"Z"+"\n")
        for i in range(0,len(angle)):
            new_file.write(str(angle[i])+"\t"+str(eq_fl[i])+"\t"+str(mean_fl[i])+"\t"+str(tur_fl[i])+"\t"+str(ttf[i])+"\t"+str(total[i])+"\t"+str(check[i])+"\t"+str(R_value[i])+"\t"+str(Z_value[i])+"\n")
        new_file.close()


    #return angle, eq_fl, tur_fl, Rnorm, R



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





#angle, eq_fl, tur_fl, r, R = file_read(True)
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

