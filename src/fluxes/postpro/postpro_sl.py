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

#set up parameters for plotting:
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12


#reads all the files and joins the values in lists. Then, it filters duplicates and nans.
def file_read(verbose = False):
    tolerance = 0.5
    crs = []
    for i in range(0,60): #ti253-separatrix(0,60)/ti253-IFS(0,61)/ti255-separatrix (0,60)/ti262 (0,78)
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_Scale_lengths_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_temperatures_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_densities_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_phis_vs.angle_%d.txt" %(i),"r"))
        crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_magnetic_field_vs.angle_%d.txt" %(i),"r"))
    for i in range(0,60):
        next(crs[i])
    angle = []
    Bp = []
    B_mag = []
    #phi = []
    #dphi = []
    #den = []
    #lte = []
    #lti = []
    #ln = []
    #R_value = []
    #Z_value = []
    for i in range(0,60):
        for columns in ( raw.strip().split() for raw in crs[i] ):  
            angle.append(float(columns[0]))
            Bp.append(float(columns[1]))
            B_mag.append(float(columns[2]))
            #phi.append(float(columns[1]))
            #dphi.append(float(columns[2]))
            #den.append(float(columns[1]))
            #lte.append(float(columns[1]))
            #lti.append(float(columns[2]))
            #ln.append(float(columns[3]))
            #R_value.append(float(columns[4]))
            #Z_value.append(float(columns[5]))
   
    excl_val_2 = []
    for elem in np.where(np.isnan(Bp))[0][:]:#filter out nan values from fluxes
        excl_val_2.append(elem)
    for elem in sorted(excl_val_2,reverse=True):
        del angle[elem]
        del Bp[elem]
        del B_mag[elem]
        #del phi[elem]
        #del dphi[elem]
        #del den[elem]
        #del lte[elem]
        #del lti[elem]
        #del ln[elem]
        #del R_value[elem]
        #del Z_value[elem]
   

    excl_val_3 = []
    for elem in np.where(np.isnan(B_mag))[0][:]:
        excl_val_3.append(elem)
    for elem in sorted(excl_val_3,reverse=True):
        del angle[elem]
        del Bp[elem]
        del B_mag[elem]
        #del phi[elem]
        #del dphi[elem]
        #del lte[elem]
        #del lti[elem]
        #del ln[elem]
        #del R_value[elem]
        #del Z_value[elem]
        
    #excl_val_4 = []
    #for elem in np.where(np.isnan(ln))[0][:]:
        #excl_val_4.append(elem)
    #for elem in sorted(excl_val_4,reverse=True):
        #del angle[elem]
        #del lte[elem]
        #del lti[elem]
        #del ln[elem]
        #del R_value[elem]
        #del Z_value[elem]    
    
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
                del Bp[elem]
                del B_mag[elem]
                #del phi[elem]
                #del dphi[elem]
                #del den[elem]
                #del lte[elem]
                #del lti[elem]
                #del ln[elem]
                #del R_value[elem]
                #del Z_value[elem]
        
    nul_pos = angle.index(min(angle))#start all values from zero angle
    angle = np.roll(angle,-nul_pos)
    Bp = np.roll(Bp,-nul_pos)
    B_mag = np.roll(B_mag,-nul_pos)
    #phi = np.roll(phi, -nul_pos)
    #dphi = np.roll(dphi, -nul_pos)
    #den = np.roll(den, -nul_pos)
    #lte = np.roll(lte,-nul_pos)
    #lti = np.roll(lti,-nul_pos)
    #ln = np.roll(ln,-nul_pos)
    #R_value = np.roll(R_value,-nul_pos)
    #Z_value = np.roll(Z_value,-nul_pos)
    print(angle)
    
    #Rnorm = []#returning the R-distance for the integral function. This is the r-distance from the center of the tokamak to the flux point. 
    #for r,z in zip(R_value,Z_value):
        #Rnorm.append(r**2 + z**2)
    #Returning the distance from the axis of the tokamak for the R dzeta part of the integral.
    #R = R_value + core.Rmaj
    #etai = lti/ln
    #etae = lte/ln
    #new_file = open("etas.txt",'a')
    #new_file.write("angle"+"\t"+"eta_e"+"\t"+"eta_i"+"\t"+"R"+"\t"+"Z"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(etae[i])+"\t"+str(etai[i])+"\t"+str(R_value[i])+"\t"+str(Z_value[i])+"\n")
    #new_file.close()

    #new_file = open("density.txt",'a')
    #new_file.write("angle"+"\t"+"density"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(den[i])+"\n")
    #new_file.close()
    
    #new_file = open("phis_eq.txt",'a')
    #new_file.write("angle"+"\t"+"phi"+"\t"+"dphi"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(phi[i])+"\t"+str(dphi[i])+"\n")
    #new_file.close()
    
    new_file = open("Sep_Mag_field.txt",'a')
    new_file.write("angle"+"\t"+"Bp"+"\t"+"B"+"\n")
    for i in range(0,len(angle)):
        new_file.write(str(angle[i])+"\t"+str(Bp[i])+"\t"+str(B_mag[i])+"\n")
    new_file.close()
    
    
    #plt.style.use('ggplot')
    #plt.title("Stability Parameters")
    #plt.xlabel("angle")
    #plt.ylabel(r"$\eta$")
    #plt.xlim([0.0,360.0])
    #plt.ylim([0.0,3.0])
    #plt.plot(angle[:], etai[:],'b',marker='o',label=r'$\eta_i$')
    #plt.plot(angle[:], etae[:],'r',marker='o',label=r'$\eta_e$')
    #plt.legend()
    #plt.show()

    
    
    if verbose == True:# Plotting option
        plt.style.use('ggplot')
        plt.title("Inverse Scale Lenghts")
        plt.xlabel("angle")
        plt.ylabel(r"$L^{-1}(m)^{-1}$")
        plt.xlim([0.0,360.0])
        plt.plot(angle[:], lti[:],'b',marker='o',label=r'$L^{-1}_{T_i}$')
        plt.plot(angle[:], lte[:],'r',marker='o',label=r'$L^{-1}_{T_e}$')
        plt.plot(angle[:], ln[:],'g',marker='o',label=r'$L^{-1}_{n}$')
        plt.legend()
        plt.show()

        new_file = open("Scale_lengths_2.txt",'a')
        new_file.write("angle"+"\t"+"LTe"+"\t"+"LTi"+"\t"+"Ln"+"\t"+"R"+"\t"+"Z"+"\n")
        for i in range(0,len(angle)):
            new_file.write(str(angle[i])+"\t"+str(lte[i])+"\t"+str(lti[i])+"\t"+str(ln[i])+"\t"+str(R_value[i])+"\t"+str(Z_value[i])+"\n")
        new_file.close()


    return angle#, lte, lti, ln, Rnorm, R
file_read()