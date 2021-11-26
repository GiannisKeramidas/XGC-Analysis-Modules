import numpy as np
import xgc
import math
import pylab as py
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev, splprep, interp1d
from scipy.optimize import curve_fit
import sys
import numpy.ma as ma
from decimal import Decimal
import core
import angle

phi_start=0
phi_end=None


def main_Ad(Rmin):
    angle.angle_mesh()
    file = open("ti255_Adiabatic_fluxes_vs.angle.txt",'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Eq_Gamma_time_avg, Tur_Gamma_time_avg = angle.two_fluxes_angle(i,i,Ang,0,Rmin,1)
        file.write(str(Ang)+"\t"+str(Eq_Gamma_time_avg)+"\t"+str(Tur_Gamma_time_avg)+"\t"+str(R_val) + "\t"+ str(Z_val)+ "\n")
    file.close()

def main_NAd(Rmin):
    angle.angle_mesh()
    file = open("ti255_NAdiabatic_fluxes_vs.angle.txt",'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Eq_Gamma_time_avg, Tur_Gamma_time_avg = angle.two_fluxes_angle(i,i,Ang,0,Rmin,1)
        file.write(str(Ang)+"\t"+str(Eq_Gamma_time_avg)+"\t"+str(Tur_Gamma_time_avg)+"\t"+str(R_val) + "\t"+ str(Z_val)+ "\n")
    file.close()

def main_Tot(Rmin):
    angle.angle_mesh()
    file = open("ti255_Total_fluxes_vs.angle.txt",'a')
    #calls flux function on diagonal points of the mesh
    for i in range(1,len(angle.R_fp)-1):
        R_val = angle.R_fp[i]
        Z_val = angle.Z_fp[i]
        Ang = angle.norm_atan(Z_val,R_val)
        Eq_Gamma_time_avg, Tur_Gamma_time_avg = angle.two_fluxes_angle(i,i,Ang,0,Rmin,3)
        file.write(str(Ang)+"\t"+str(Eq_Gamma_time_avg)+"\t"+str(Tur_Gamma_time_avg)+"\t"+str(R_val) + "\t"+ str(Z_val)+ "\n")
    file.close()

core.getMeshAndCuts(core.fileDir,core.Rmin,core.Rmax,core.Zmin,core.Zmax)
main_Ad(core.Rmin)
main_NAd(core.Rmin)
main_Tot(core.Rmin)
