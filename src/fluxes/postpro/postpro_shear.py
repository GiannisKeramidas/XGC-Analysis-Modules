import numpy as np
import xgc
import math
import pylab as py
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys
import numpy.ma as ma
import core
from decimal import Decimal
import angle

'''Module that contains functions that will take as input shear arrays (2D arrays of shear values on an RZ plane-toroidally averaged values    have been used in the shear calculation-they have been calculated using a parallel program), Rsep and Zsep (also calculated parallely)      arrays and will produce a text file with the shear values on the separatrix. List of things that need to be checked before running          different cases:
   1) Correct separatrix file at sep_file_read
   2) Correct cut-off to avoid X-point at loading_checks
   3) Make sure that shear arrays with proper numbers exist as well as Rsep, Zsep
   4) Write correct run number on filename at sh_on_sep, file_read and plotting
'''

def sep_file_read():
    '''Function that reads the separatrix node locations'''
    #sep_file = open("/global/homes/g/giannos/xgc_python_dir/ti253_analysis/ti253_sep_flux/Separatrix nodes locations","r")
    #sep_file = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/Sep-ti255_med.txt","r")
    sep_file = open("/global/homes/g/giannos/xgc_python_dir/ti262_analysis/ti262_Separatrix.txt","r")
    #sep_file = open("IFSNL.txt","r")
    next(sep_file)
    R=[]
    Z=[]
    for columns in (raw.strip().split() for raw in sep_file):
        R.append(columns[0])
        Z.append(columns[1])
    return R,Z

def loading_checks(Rmin,Rmax,Zmin,Zmax):
    '''Originally checks for reversal of min and max. Here is used to extract the min's and max's from the parallel loading method.'''
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

    if abs(Rmin-Rmax)<0.12:#avoid loading too narrow windows
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

Rsep,Zsep = sep_file_read()

def sh_on_sep(r):
    '''For a particular process number, reads in the matrices of R,Z separatrix points, the respective shear matrix and writes out a file
    that contains the R,Z,shear information.'''
    #reading the matrices
    Rpol = np.load("Rsep_%s.npy"%(r))
    Zpol = np.load("Zsep_%s.npy"%(r))
    arr = np.load("shear_%s.npy"%(r))
    
    #finding min/max of R,Z.
    if ((r*20)<(len(Rsep)-1))and (((r+1)*20)<(len(Rsep)-1)):
        Rmin = float(Rsep[r*20])#2.14#
        Rmax = float(Rsep[(r+1)*20])#2.26#
        Zmin = float(Zsep[r*20])#-0.38#
        Zmax = float(Zsep[(r+1)*20])#0.45#
        Rmin,Rmax,Zmin,Zmax = loading_checks(Rmin,Rmax,Zmin,Zmax)
    
    elif ((r*20)<(len(Rsep)-1))and (((r+1)*20)>(len(Rsep)-1)):
        Rmin = float(Rsep[r*20])
        Rmax = float(Rsep[-1])
        Zmin = float(Zsep[r*20])
        Zmax = float(Zsep[-1])
    else:
        pass
    unitR = (Rmax-Rmin)/100
    unitZ = (Zmax-Zmin)/100
    
    Zs = np.array([x for x in range(0,arr.shape[0])])
    Rs = np.array([x for x in range(0,arr.shape[1])])
    file = open("shear_files_%s.txt" %(r),'a')
    #file = open("Reynolds_turpf_mid.txt",'a')
    #file = open("Reynolds_force_midplane.txt",'a')
    #file = open("pol_flow_mid.txt",'a')
    #file = open("Re_rad_f_mid.txt",'a')
    #file = open("Re_rad_f_mid.txt",'a')
    Septrix = list(zip(Rpol,Zpol))
    Z_coord = []
    Z_loc = []
    R_loc = []
    R_coord = []
    for Z in Zpol:
        diff=[]
        for zeta in Zs:
            diff.append(abs(Z-(zeta*unitZ+Zmin)))#find units and min/max from parallel.
        tmp = np.amin(diff)
        Sep_loc_z = diff.index(tmp)
        Z_loc.append(Sep_loc_z)
        Z_coord.append(Zs[Sep_loc_z]*unitZ+Zmin)
    for R in Rpol:
        diff=[]
        for rho in Rs:
            diff.append(abs(R-(rho*unitR+Rmin)))
        Sep_loc_r = np.argmin(diff)
        R_loc.append(Sep_loc_r)
        R_coord.append(Rs[Sep_loc_r]*unitR+Rmin-core.Rmaj)
    #plt.plot(R_coord,Z_coord,'r--',Rpol,Zpol,'b*')
    #plt.show()
    #file.write("R"+"\t"+"Z"+"\t"+"Shear"+"\n")
    #file.write("R"+"\t"+"Z"+"\t"+"Re_S"+"\n")
    file.write("angle"+"\t"+"R"+"\t"+"Z"+"\t"+"shear"+"\n")
    coordinates = list(zip(R_coord,Z_coord))
    ANG = []
    for i in range(0,len(coordinates)):
        angie = angle.norm_atan(Z_coord[i],R_coord[i])
        ANG.append(angie)
    for i in range(0,len(coordinates)):
        if arr[Z_loc[i],R_loc[i]] > 0:
            file.write(str(ANG[i])+"\t"+str(coordinates[i][0])+"\t"+str(coordinates[i][1])+"\t"+str(arr[Z_loc[i],R_loc[i]])+"\n")
        else:
            pass
    file.close()
    
'''Smoothing algorithms'''

def smoothing_alg(L):
    LF = np.roll(L,1)
    LB = np.roll(L,-1)
    LFF = np.roll(L,2)
    LBB = np.roll(L,-1)
    newL = (2.0*(LF + LB)+LFF+LBB)/6.0
    return newL

def smooth(L):
    sL = L
    for i in range(0,100):
        sL = smoothing_alg(sL)
    return sL

def file_read(verbose = True):
    '''Reads all the files and joins the values in lists. Then, it filters duplicates and nans.'''
    tolerance = 0.2
    crs = []
    for i in range(0,77):#ti253-separatrix(0,60)/ti253-IFS(0,61)/ti255-separatrix (0,60)
        crs.append(open("/global/homes/g/giannos/xgc_python_dir/shear_files_%d.txt" %(i),"r"))
    for i in range(0,77):
        next(crs[i])
    angle = []
    R_value = []
    Z_value = []
    shear = []
    for i in range(0,77):
        try:
            for columns in (raw.strip().split() for raw in crs[i]):  
                angle.append(float(columns[0]))
                R_value.append(float(columns[1]))
                Z_value.append(float(columns[2]))
                shear.append(float(columns[3]))
        except ValueError:
            print("error on file",i)
   
    excl_val_2 = []
    for elem in np.where(np.isnan(shear))[0][:]:#filter out nan values from shear
    #for elem in np.where(shear < 1e2)[0][:]:
        excl_val_2.append(elem)
        #print(np.where(shear == 0.0))
    for elem in sorted(excl_val_2,reverse=True):
        del angle[elem]
        del R_value[elem]
        del Z_value[elem]
        del shear[elem]
    
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
                del B[elem]
                del angle[elem]
                del shear[elem]
                del R_value[elem]
                del Z_value[elem]
        
    nul_pos = angle.index(min(angle))#start all values from zero angle
    angle = np.roll(angle,-nul_pos)
    shear = np.roll(shear,-nul_pos)
    R_value = np.roll(R_value,-nul_pos)
    Z_value = np.roll(Z_value,-nul_pos)
    #print(angle)
    
    Rnorm = []#returning the R-distance for the integral function.
    for r,z in zip(R_value,Z_value):
        Rnorm.append(r**2 + z**2)

    R = R_value + core.Rmaj

    if verbose == True:# Plotting option
        #new_angle = smooth(angle)
        fig,ax = plt.subplots()
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
        ax.set_xlim([0,360])
        ax.set_xlabel("angle")
        ax.set_ylabel("shear rate")
        ax.plot(angle[:], shear[:],'b-')
        ax.plot(np.unique(angle), np.poly1d(np.polyfit(angle, shear, 18))(np.unique(angle)),'r')
        plt.show()

    new_file = open("Shear_total.txt",'a')
    new_file.write("angle"+"\t"+"Shear"+"\t"+"R"+"\t"+"Z"+"\n")
    for i in range(0,len(angle)):
        new_file.write(str(angle[i])+"\t"+str(shear[i])+"\t"+str(R_value[i])+"\t"+str(Z_value[i])+"\n")
    new_file.close()

    return angle, shear, Rnorm, R

def plotting():
    '''Plotting function'''
    crs = open("Shear_total.txt",'r')
    angle = []
    shear = []
    for columns in (raw.strip().split() for raw in crs):  
        angle.append(float(columns[0]))
        shear.append(float(columns[1]))
    fig,ax = plt.subplots()
    ax.set_xlim([0,360])
    ax.xlabel("angle")
    ax.ylabel("shear rate")
    ax.plot(angle[:], shear[:],'b-')
    plt.show()
            
'''First run the for loop to create the text files, then run the file_read to join them in one. Plotting is optional.'''
#for r in range(0,77):
   # sh_on_sep(r)
file_read()
#plotting()

