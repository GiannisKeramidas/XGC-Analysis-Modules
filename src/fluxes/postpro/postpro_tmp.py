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
def file_read(verbose = True):
    tolerance = 0.5
    crs = []
    for i in range(0,60):#ti253-separatrix(0,60)/ti253-IFS(0,61)/ti255-separatrix (0,60)/ti262 (0,78)
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_temperatures_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_densities_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_turbulence_strength_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_magnetic_field_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_geom_magnetic_drift_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_ion_temperatures_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_magnetic_field_comp_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_magnetic_drifts_vs.angle_%d.txt" %(i),"r"))
        #crs.append(open("/global/homes/g/giannos/xgc_python_dir/ti255_density_perturbation_strength_vs.angle_%d.txt" %(i),"r"))
        crs.append(open("/global/homes/g/giannos/xgc_python_dir/Check_%d.txt" %(i),"r"))
    #for i in range(0,60):
        #next(crs[i])
    
        
    angle = []
    #gradBR = []
    #gradBZ = []
    #curvR = [] 
    #curvZ = []
    #Br = []
    #Bz = []
    #T_perp = []
    #E_par = []
    #geom = []
    #Bp = []
    #Bm = []
    #den=[]
    #te = []
    #ti = []
    #den_str = []
    R = []
    Z = []
    for i in range(0,60):
        for columns in ( raw.strip().split() for raw in crs[i] ):  
            angle.append(float(columns[0]))
            R.append(float(columns[1]))
            Z.append(float(columns[2]))
            #den_str.append(float(columns[1]))
            #gradBR.append(float(columns[1]))
            #gradBZ.append(float(columns[2]))
            #curvR.append(float(columns[3]))
            #curvZ.append(float(columns[4]))
            #Br.append(float(columns[1]))
            #Bz.append(float(columns[2]))
            #T_perp.append(float(columns[1]))
            #E_par.append(float(columns[2]))
            #geom.append(float(columns[1]))
            #Bp.append(float(columns[1]))
            #Bm.append(float(columns[2]))
            #den.append(float(columns[1]))
            #ti.append(float(columns[1]))
            #te.append(float(columns[2]))
            
   
    excl_val_2 = []
    for elem in np.where(np.isnan(R))[0][:]:#filter out nan values from fluxes
        excl_val_2.append(elem)
    for elem in sorted(excl_val_2,reverse=True):
        del angle[elem]
        del R[elem]
        del Z[elem]
        #del den_str[elem]
        #del gradBR[elem]
        #del gradBZ[elem]
        #del curvR[elem]
        #del curvZ[elem]
        #del Br[elem]
        #del Bz[elem]
        #del Bzeta[elem]
        #del T_perp[elem]
        #del E_par[elem]
        #del geom[elem]
        #del Bp[elem]
        #del Bm[elem]
        #del den[elem]
        #del ti[elem]
        #del te[elem]
        

    excl_val_3 = []
    for elem in np.where(np.isnan(Z))[0][:]:
        excl_val_3.append(elem)
    for elem in sorted(excl_val_3,reverse=True):
        del angle[elem]
        del R[elem]
        del Z[elem]
        #del gradBR[elem]
        #del gradBZ[elem]
        #del curvR[elem]
        #del curvZ[elem]
        #del Br[elem]
        #del Bz[elem]
        #del Bzeta[elem]
        #del T_perp[elem]
        #del E_par[elem]
        #del Bp[elem]
        #del Bm[elem]
        #del ti[elem]
        #del te[elem]
        
    #excl_val_4 = []
    #for elem in np.where(np.isnan())[0][:]:
        #excl_val_4.append(elem)
    #for elem in sorted(excl_val_4,reverse=True):
        #del angle[elem]
        #del Br[elem]
        #del Bz[elem]
        #del Bzeta[elem]     
            
    
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
                del R[elem]
                del Z[elem]
                #del den_str[elem]
                #del gradBR[elem]
                #del gradBZ[elem]
                #del curvR[elem]
                #del curvZ[elem]
                #del Br[elem]
                #del Bz[elem]
                #del Bzeta[elem]
                #del T_perp[elem]
                #del E_par[elem]
                #del geom[elem]
                #del Bp[elem]
                #del Bm[elem]
                #del den[elem]
                #del ti[elem]
                #del te[elem]
                
        
    nul_pos = angle.index(min(angle))#start all values from zero angle
    angle = np.roll(angle,-nul_pos)
    R = np.roll(R,-nul_pos)
    Z = np.roll(Z,-nul_pos)
    #den_str = np.roll(den_str,-nul_pos)
    #gradBR = np.roll(gradBR,-nul_pos)
    #gradBZ = np.roll(gradBZ,-nul_pos)
    #curvR = np.roll(curvR,-nul_pos)
    #curvZ = np.roll(curvZ,-nul_pos)
    #Br = np.roll(Br,-nul_pos)
    #Bz = np.roll(Bz,-nul_pos)
    #Bzeta = np.roll(Bzeta,-nul_pos)
    #T_perp = np.roll(T_perp,-nul_pos)
    #E_par = np.roll(E_par,-nul_pos)
    #geom = np.roll(geom,-nul_pos)
    #Bp = np.roll(Bp,-nul_pos)
    #Bm = np.roll(Bm,-nul_pos)
    #den = np.roll(den,-nul_pos)
    #ti = np.roll(ti,-nul_pos)
    #te = np.roll(te,-nul_pos)
    print(angle)
   
    #new_file = open("Sep_temperatures.txt",'a')
    #new_file.write("angle"+"\t"+"t_i"+"\t"+"t_e"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(ti[i])+"\t"+str(te[i])+"\n")
    #new_file.close()

    #new_file = open("Sep_densities.txt",'a')
    #new_file.write("angle"+"\t"+"den"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(den[i])+"\n")
    #new_file.close()
    
    #new_file = open("Sep_turbulence_strength.txt",'a')
    #new_file.write("angle"+"\t"+"rms(dphi)/<phi>"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(den[i])+"\n")
    #new_file.close()
    
    #new_file = open("Sep_Mag_field.txt",'a')
    #new_file.write("angle"+"\t"+"B_p"+"\t"+"B"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(Bp[i])+"\t"+str(Bm[i])+"\n")
    #new_file.close()
    
    #new_file = open("Sep_Mag_geom.txt",'a')
    #new_file.write("angle"+"\t"+"geom"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(geom[i])+"\n")
    #new_file.close()
    
    #new_file = open("Ion_temperatures.txt",'a')
    #new_file.write("angle"+"\t"+"T_perp"+"\t"+"E_para"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(T_perp[i])+"\t"+str(E_par[i])+"\n")
    #new_file.close()
    
    #new_file = open("Mag_components.txt",'a')
    #new_file.write("angle"+"\t"+"Br"+"\t"+"Bz"+"\t"+"Bzeta"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(Br[i])+"\t"+str(Bz[i])+"\n")
    #new_file.close()
    
    #new_file = open("Mag_curv_drifts.txt",'a')
    #new_file.write("angle"+"\t"+"curvR"+"\t"+"curvZ"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(curvR[i])+"\t"+str(curvZ[i])+"\n")
    #new_file.close()
    
    #new_file = open("Den_pert_str.txt",'a')
    #new_file.write("angle"+"\t"+"dn"+"\n")
    #for i in range(0,len(angle)):
        #new_file.write(str(angle[i])+"\t"+str(den_str[i])+"\n")
    #new_file.close()

    new_file = open("Check.txt",'a')
    new_file.write("angle"+"\t"+"R"+"\t"+"Z"+"\n")
    for i in range(0,len(angle)):
        new_file.write(str(angle[i])+"\t"+str(R[i])+"\t"+str(Z[i])+"\n")
    new_file.close()

    
    #plt.style.use('ggplot')
    #plt.xlabel(r"angle$^{\circ}$")
    #plt.ylabel(r"T (eV)")
    #plt.xlim([0.0,360.0])
    #plt.plot(angle[:], ti[:],'b',marker='o',label=r'$T_i$')
    #plt.plot(angle[:], te[:],'r',marker='o',label=r'$T_e$')
    #plt.legend()
    #plt.show()

    #plt.style.use('ggplot')
    #plt.xlabel(r"angle$^{\circ}$")
    #plt.ylabel(r"n ($m^{-3}$)")
    #plt.xlim([0.0,360.0])
    #plt.plot(angle[:], den[:],'b',marker='o')
    #plt.legend()
    #plt.show()

    return angle

file_read()
