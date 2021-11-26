import numpy as np
import xgc
import math
import pylab as py
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import postpro

def total_read():
    #total_file_sep = open("Sep_fluxes_total_new.txt",'r')
    #total_file_sep = open("ti255_shorter_run_total.txt",'r')
    #total_file_sep = open("ti255_Sep_fluxes_total.txt",'r')
    #total_file_ifs = open("IFS_fluxes_total.txt",'r')
    #total_file_sep = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/sep_flux/long_run/ti255_Sep_fluxes_total.txt","r")
    total_file_sep = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/sep_flux/short_run/ti255_Short_fluxes_total.txt","r")

    next(total_file_sep)
    angle_s=[]
    Eq_s=[]
    Tur_s=[]
    R_s=[]
    Z_s=[]

    for columns in (raw.strip().split() for raw in total_file_sep):  
        angle_s.append(float(columns[0]))
        Eq_s.append(float(columns[1]))
        Tur_s.append(float(columns[2]))
        R_s.append(float(columns[3]))
        Z_s.append(float(columns[4]))
    Rnorm = []
    for r,z in zip(R_s,Z_s):
        Rnorm.append(r**2+z**2)

    return angle_s,Eq_s,Tur_s,R_s,Z_s,Rnorm

#angle,Eq,Tur,R,Z,Rnorm = total_read()

def smoothed_file():
    smooth = open("Smoothed_Eq_flux.txt","a")
    new_Eq = savgol_filter(Eq, 5, 2)
    for val in new_Eq:
        smooth.write(str(val)+"\n")
    smooth.close()



def eq_plot():# creates arrow plots of Equilibrium flux. Has option for even/uneven aspect ratio.
    new_an = np.arctan2(Z,R)
    Eqx,Eqy = Eq*np.cos(new_an), Eq*np.sin(new_an)
    plt.title('Equilibrium Flux around the Separatrix')
    #plt.title('Equilibrium Flux around the Inner Flux Surface')
    plt.quiver(R, Z, Eqx, Eqy,Eq, alpha=.5)
    plt.quiver(R, Z, Eqx, Eqy,edgecolor='k', facecolor='None', linewidth=.5)
    plt.xlabel('R(m)')
    plt.ylabel('Z(m)')
    plt.xlim(-0.6,0.9)
    plt.ylim(-1.57,1.1)
    plt.axes().set_aspect('equal')
    plt.show()


def tur_plot():# creates arrow plots of Turbulent flux
    Turx,Tury = Tur*np.cos(new_an), Tur*np.sin(new_an)
    plt.title('Turbulent Flux around the Separatrix')
    #plt.title('Turbulent Flux around the Inner Flux Surface')
    plt.quiver(R, Z, Turx, Tury,Tur, alpha=.5)
    plt.quiver(R, Z, Turx, Tury,edgecolor='k', facecolor='None', linewidth=.5)
    plt.xlabel('R(m)')
    plt.ylabel('Z(m)')
    plt.xlim(-0.6,0.79)
    plt.ylim(-1.3,1.3)
    #plt.axes().set_aspect('equal')
    plt.show()

def smooth_plot():#Smoothes Equilibrium data using a Savitzky-Golay filter
    new_Eq = savgol_filter(Eq, 5, 2)

    #creates arrow plots of Smoothed Equlibrium flux
    filt_Eqx,filt_Eqy = new_Eq*np.cos(new_an), new_Eq*np.sin(new_an)
    plt.title('Smoothed Equilibrium Flux around the Separatrix')
    #plt.title('Smoothed Equilibrium Flux around the Inner flux surface')
    plt.quiver(R, Z, filt_Eqx, filt_Eqy,new_Eq, alpha=.5)
    plt.quiver(R, Z, filt_Eqx, filt_Eqy,edgecolor='k', facecolor='None', linewidth=.5)
    plt.xlabel('R(m)')
    plt.ylabel('Z(m)')
    plt.xlim(-0.6,0.87)
    plt.ylim(-1.57,1.1)
    plt.axes().set_aspect('equal')
    plt.show()

def int_check():
    angle, eq_fl, tur_fl,_,_, Rnorm = total_read()
    R = np.asarray(Rnorm)
    print(type(R))
    print(type(angle))
    new_Eq = savgol_filter(eq_fl,5,2)
    old_eq = postpro.Line_int_along_curve(R,angle,eq_fl)
    new_eq = postpro.Line_int_along_curve(R,angle,new_Eq)
    print("Old value:",old_eq)
    print("New value",new_eq)
     
