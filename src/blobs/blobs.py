import numpy as np
import xgc
import math
from scipy.interpolate import interp2d
from scipy.interpolate import RegularGridInterpolator
import scipy
import scipy.integrate
from scipy.integrate import odeint
import sqlite3
import sys
import numpy.ma as ma
import os
import pandas as pd
from skimage import measure
from numpy.linalg import eig, inv
from matplotlib import animation, rc
import matplotlib.animation as animation
from matplotlib.ticker import FormatStrFormatter
import types
import random
from matplotlib import colors as mcolors
from IPython.display import HTML
import warnings
import seaborn as sns
import matplotlib
import pylab as py
import matplotlib.pyplot as plt
from scipy.stats import kurtosis
from scipy.stats import skew
from scipy.stats import tstd
plt.style.use('ggplot')


"""Discharge, frame and database parameters."""
"""CMOD"""
dt = 3.113e-8 #9.454e-8 #(s)
Rmin = 0.872 #0.88 # (m)
Rmax = 0.90 #0.908 # (m)
Zmin = -0.1 # (m)
Zmax = 0.1 # (m)
big_Rmax = 0.909968
big_Rmin = 0.4402
big_Zmax = 0.4328
big_Zmin = -0.5959

Nplanes = 16

#R_sep_list = contours.R_sep[305:-1139]
#Z_sep_list = contours.Z_sep[305:-1139]

#magnetic field of the small frame
#BR = np.load('/global/cscratch1/sd/giannos/ti262_arrays/B_R_mid.npy')
#BZ = np.load('/global/cscratch1/sd/giannos/ti262_arrays/B_Z_mid.npy')
#Bzeta = np.load('/global/cscratch1/sd/giannos/ti262_arrays/B_zeta_mid.npy')

BR = np.load('/global/cscratch1/sd/giannos/ti344_arrays/B_R_mid.npy')
BZ = np.load('/global/cscratch1/sd/giannos/ti344_arrays/B_Z_mid.npy')
Bzeta = np.load('/global/cscratch1/sd/giannos/ti344_arrays/B_zeta_mid.npy')


#field line coordinates of next and previous frame
#Rp = np.load('/global/cscratch1/sd/giannos/ti262_arrays/Rp1.npy')
#Zp = np.load('/global/cscratch1/sd/giannos/ti262_arrays/Zp1.npy')
#Rm = np.load('/global/cscratch1/sd/giannos/ti262_arrays/Rm1.npy')
#Zm = np.load('/global/cscratch1/sd/giannos/ti262_arrays/Zm1.npy')

#psi array (it's usually masked and we have to reconstruct it)
psi = np.load('/global/cscratch1/sd/giannos/ti344_arrays/psi.npy')
psi_mask = np.load('/global/cscratch1/sd/giannos/ti344_arrays/mask_psi.npy')
psi = ma.array(psi,mask=psi_mask)
psi = np.swapaxes(psi,0,1)

#rho_s = 0.09 # In (cm)
#rho_i = 0.04#0.07 # In (cm)
#delta_star = 0.44 # In (cm)
#alpha = 0.5 # In (cm)
#k_perp = 314 # In (m^{-1})

#db_name = "/global/cscratch1/sd/giannos/ti262_arrays/blobs_rc"
db_name = "/global/cscratch1/sd/giannos/ti344_arrays/blobs3"


"""DIII-D"""
"""dt = 2.33357e-7 #(s)
Rmin = 2.19 #(m)
Rmax = 2.27 #(m)
Zmin = -0.25 #(m)
Zmax = 0.25 #(m)
Nplanes = 32
big_Rmax = 2.377
big_Rmin = 1.001
big_Zmax = 1.348
big_Zmin = -1.363

#R_sep_list = contours.R_sep
#Z_sep_list = contours.Z_sep

# Magnetic field of the small frame
BR = np.load('/global/cscratch1/sd/giannos/ti255_arrays/B_R_mid.npy')
BZ = np.load('/global/cscratch1/sd/giannos/ti255_arrays/B_Z_mid.npy')
Bzeta = np.load('/global/cscratch1/sd/giannos/ti255_arrays/B_zeta_mid.npy')

# Field line coordinates of next and previous frame
Rp = np.load('/global/cscratch1/sd/giannos/ti255_arrays/Rp.npy')
Zp = np.load('/global/cscratch1/sd/giannos/ti255_arrays/Zp.npy')
Rm = np.load('/global/cscratch1/sd/giannos/ti255_arrays/Rm.npy')
Zm = np.load('/global/cscratch1/sd/giannos/ti255_arrays/Zm.npy')

#psi array
psi = np.load('/global/cscratch1/sd/giannos/ti255_arrays/psi_mid.npy')
psi = np.swapaxes(psi,0,1)

rho_s = 0.19 # In (cm)
rho_i = 0.2 # In (cm)
delta_star = 1.61 # In (cm)
alpha = 1.27 # In (cm)
k_perp = 123 # In (m^{-1})

db_name = '/global/cscratch1/sd/giannos/ti255_arrays/blobs_fls'
"""
"""PPinch"""
'''#dt = 1.5656919e-7 #(s)
#big_Rmin = 1.016 #(m)
#big_Rmax = 2.3767 #(m)
#big_Zmin = -1.36 #(m)
#big_Zmax = 1.35 #(m)
#Rmin = 2.23
#Rmax = 2.31
#Zmin = -0.25
#Zmax = 0.4
#Nplanes = 16
#maxR_abs = 2.3767
#minR_abs = 1.016
#maxZ_abs = 1.3515
#minZ_abs = -1.3675

#magnetic field of the small frame
#BR = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\B_R_mid_pinch.npy')
#BZ = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\B_Z_mid_pinch.npy')
#Bzeta = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\B_zeta_mid_pinch.npy')

#field line coordinates of next and previous frame
#Rp = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\Rp.npy')
#Zp = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\Zp.npy')
#Rm = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\Rm.npy')
#Zm = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\Zm.npy')

#database
#db_name = 'blobs_pinch'
#db_name = 'blobs_pinch2'
'''
"""Common parameters and units calculations"""
Bp = np.sqrt(np.power(BR,2)+np.power(BZ,2))

# The frames have a resolution of 100x100 but the magnetic field arrays for the full tokamak have a resolution of 1000x1000.
Lr = np.array([x for x in range(0,100)])
unitR = (Rmax-Rmin)/len(Lr)
R_points = np.array([Lr[ip]*unitR+Rmin for ip in range(0,len(Lr))])
Lz = np.array([x for x in range(0,100)])
unitZ = (Zmax-Zmin)/len(Lz)
Z = np.array([y for y in range(0,100)])
Z_points = np.array([Z[ip]*unitZ+Zmin for ip in range(0,len(Lz))])
abs_unitR = (big_Rmax-big_Rmin)/1000
abs_unitZ = (big_Zmax-big_Zmin)/1000

'''Functions that returns real space coordinates from grid coordinates x, y.'''
def R_space(x):
    
    return Rmin + unitR*x

def Z_space(y):
    
    return Zmin + unitZ*y


"""Plotting helping functions."""
'''Inverse transformation of separatrix coordinates to mesh coordinates.'''
def separatrix_transform(R_sep_list,Z_sep_list):
    
    R_sep = [(x-Rmin)/unitR for x in R_sep_list]
    Z_sep = [(x-Zmin)/unitZ for x in Z_sep_list]
    
    return R_sep, Z_sep

'''Returns x,y axes plot labels.'''
def plot_labels(Rmin,Rmax,Zmin,Zmax,xtick_num,ytick_num):
    
    x_labels = np.linspace(Rmin,Rmax,xtick_num+1)
    x_lbl = [str("%.3f" % round(x,3))  for x in x_labels]
    y_labels = np.linspace(Zmin,Zmax,ytick_num+1)
    y_lbl = [str("%.2f" % round(y,3))  for y in y_labels]
    
    return(x_lbl,y_lbl)

"""Field Line Smoothing Routines"""
"""Preparing interpolation functions for full BR,BZ,Bzeta."""
# Loading global fields
"""PPinch"""
#BR_g = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\B_R_pinch.npy')
#BZ_g = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\B_Z_pinch.npy')
#Bzeta_g = np.load('C:\\Users\\giannis\\Desktop\\Python_scrpts\\pinch\\B_zeta_pinch.npy')
"""DIII-D"""
#BR_g = np.load('/global/cscratch1/sd/giannos/ti255_arrays/B_R_glob.npy')
#BZ_g = np.load('/global/cscratch1/sd/giannos/ti255_arrays/B_Z_glob.npy')
#Bzeta_g = np.load('/global/cscratch1/sd/giannos/ti255_arrays/B_zeta_glob.npy')
"""CMOD"""
#BR_g = np.load('/global/cscratch1/sd/giannos/ti262_arrays/B_R_glob.npy')
#BZ_g = np.load('/global/cscratch1/sd/giannos/ti262_arrays/B_Z_glob.npy')
#Bzeta_g = np.load('/global/cscratch1/sd/giannos/ti262_arrays/B_zeta_glob.npy')
# Swapaxes to make them into [R,Z] arrays
#BR_g = np.swapaxes(BR_g,0,1)
#BZ_g = np.swapaxes(BZ_g,0,1)
#Bzeta_g = np.swapaxes(Bzeta_g,0,1)
# Prepare the interpolation points
#r = [big_Rmin + abs_unitR*x for x in range(1000)] 
#z = [big_Zmin + abs_unitZ*x for x in range(1000)]
# Prepare the interpolation functions/we only use the interpolated global fields denoted by an underscore.
#B_R = RegularGridInterpolator((r,z), BR_g[:,:], method='linear', bounds_error=False, fill_value = 0)
#B_Z = RegularGridInterpolator((r,z), BZ_g[:,:], method='linear', bounds_error=False, fill_value = 0)
#B_zeta = RegularGridInterpolator((r,z), Bzeta_g[:,:], method='linear', bounds_error=False, fill_value = 0)

"""Defining the field line equations."""
def field_line(zeta, state_vec):
    
    R = state_vec[0]
    Z = state_vec[1]
    dRdzeta = R*(B_R([R,Z])[0]/B_zeta([R,Z])[0])
    dZdzeta = R*(B_Z([R,Z])[0]/B_zeta([R,Z])[0])
    
    return [dRdzeta, dZdzeta]

"""Making a change of variable s=-zeta so we can integrate backwards in zeta."""
def field_line_inverse(s, state_vec):
    
    R = state_vec[0]
    Z = state_vec[1]
    dRds = -R*(B_R([R,Z])[0]/B_zeta([R,Z])[0])
    dZds = -R*(B_Z([R,Z])[0]/B_zeta([R,Z])[0])
    
    return [dRds, dZds]        

"""Function that gives the coordinates of the previous plane that the grid points of the central plane map to, following the field lines.
"""
def backward_interpolation():
    
    solver = scipy.integrate.ode(field_line_inverse).set_integrator('dopri5')
    zeta_stop = (2.0*math.pi)/Nplanes
    dzeta = zeta_stop/10
    Rm = []
    Zm = []
    
    for i in range(100):
        print("Row",i)
        
        for j in range(100):
            R_init = Rmin + unitR*i
            Z_init = Zmin + unitZ*j
            
            state_vec0,zeta0 = [R_init,Z_init],0.0
            solver.set_initial_value(state_vec0,zeta0)
            y,t = [], []
            warnings.filterwarnings("ignore")
            
            while solver.successful() and solver.t<zeta_stop:
                solver.set_initial_value([solver.y[0],solver.y[1]],solver.t)
                solver.integrate(solver.t+dzeta)
                y.append(solver.y)
                t.append(solver.t)
            
            warnings.resetwarnings()
            y = np.array(y)
            t = np.array(t)
            Rm.append(y[-1,0])
            Zm.append(y[-1,1])
            
    return np.reshape(Rm,(100,100)), np.reshape(Zm,(100,100))

"""Function that gives the coordinates of the forward plane that the grid points of the central plane map to, following the field lines.
"""
def forward_interpolation():
    
    solver = scipy.integrate.ode(field_line).set_integrator('dopri5')
    zeta_stop = (2.0*math.pi)/Nplanes
    dzeta = zeta_stop/10
    Rp = []
    Zp = []
    
    for i in range(100):
        print("Row",i)
        
        for j in range(100):
            R_init = Rmin + unitR*i
            Z_init = Zmin + unitZ*j
            
            state_vec0,zeta0 = [R_init,Z_init],0.0
            solver.set_initial_value(state_vec0,zeta0)
            y,t = [], []
            warnings.filterwarnings("ignore")
            
            while solver.successful() and solver.t<zeta_stop:
                solver.set_initial_value([solver.y[0],solver.y[1]],solver.t)
                solver.integrate(solver.t+dzeta)
                y.append(solver.y)
                t.append(solver.t)
            
            warnings.resetwarnings()
            y = np.array(y)
            t = np.array(t)
            Rp.append(y[-1,0])
            Zp.append(y[-1,1])
            
    return np.reshape(Rp,(100,100)), np.reshape(Zp,(100,100))

def interp_values():
   
    print("Starting backward...")
    Rm,Zm = backward_interpolation()
    print("Starting forward...")
    Rp,Zp = forward_interpolation()
    np.save('/global/cscratch1/sd/giannos/ti262_arrays/Rm1', Rm)
    np.save('/global/cscratch1/sd/giannos/ti262_arrays/Zm1', Zm)
    np.save('/global/cscratch1/sd/giannos/ti262_arrays/Rp1', Rp)
    np.save('/global/cscratch1/sd/giannos/ti262_arrays/Zp1', Zp)

    
"""Takes three smoothed frames from adjacent planes and does smoothing along the field line."""
def field_line_smooth(arr_c,arr_p1,arr_m1):
    
    # Creating grid and interpolation functions
    big_unitR = (big_Rmax-big_Rmin)/arr_p1.shape[0]
    big_unitZ = (big_Zmax-big_Zmin)/arr_p1.shape[1]
    
    r = [big_Rmin + big_unitR*x for x in range(arr_p1.shape[0])]
    z = [big_Zmin + big_unitZ*x for x in range(arr_p1.shape[1])]

    interp_m1 = RegularGridInterpolator((r,z),arr_m1,bounds_error=False,fill_value = 0)
    interp_p1 = RegularGridInterpolator((r,z),arr_p1,bounds_error=False,fill_value = 0)
    
    # Smoothing coefficients and field line coordinates
    w1 = 0.5
    w2 = 0.25
    
    smooth_arr = []
    
    warnings.filterwarnings("ignore")
    
    for i in range(arr_c.shape[0]):
        for j in range(arr_c.shape[1]):
            smooth_arr.append(w1*arr_c[i,j] + w2*interp_p1([Rp[i][j], Zp[i][j]]) + w2*interp_m1([Rm[i][j], Zm[i][j]]))            
    return np.reshape(smooth_arr,(arr_c.shape[0],arr_c.shape[1]))  

'''Function to be used in the decision whether to count a value or not in the smoothing sum.'''
def sg(x):
    
    return abs(np.sign(x))

'''Nine-point smoothing.'''
def frame_smooth(arr):
    
    # User-defined weights:
    w1 = 2.
    w2 = 1.
    w3 = 0.75
    tally1 = w1 + w2*4 + w3*4 #for center
    tally2 = w1 + w2*3 + w3*2 #for edges
    tally3 = w1 + w2*2 +w3 #for corners
    smooth_arr = []
    
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if i-1>=0 and j-1>=0 and i+1<arr.shape[0] and j+1<arr.shape[1]:#center
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i-1,j]+arr[i+1,j]+arr[i,j-1]+arr[i,j+1])+w3*(arr[i-1,j-1]+arr[i+1,j-1]+arr[i-1,j+1]+arr[i+1,j+1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i-1,j])+sg(arr[i+1,j])+sg(arr[i,j-1])+sg(arr[i,j+1]))+w3*(sg(arr[i-1,j-1])+sg(arr[i+1,j-1])+sg(arr[i-1,j+1])+sg(arr[i+1,j+1]))))
            if i-1<0 and j-1>=0 and i+1<arr.shape[0] and j+1<arr.shape[1]:#leftmost column
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i+1,j]+arr[i,j-1]+arr[i,j+1])+w3*(arr[i+1,j-1]+arr[i+1,j+1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i+1,j])+sg(arr[i,j-1])+sg(arr[i,j+1]))+w3*(sg(arr[i+1,j-1])+sg(arr[i+1,j+1]))))
            if i-1>=0 and j-1<0 and i+1<arr.shape[0] and j+1<arr.shape[1]:#uppermost row
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i+1,j]+arr[i-1,j]+arr[i,j+1])+w3*(arr[i-1,j+1]+arr[i+1,j+1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i+1,j])+sg(arr[i-1,j])+sg(arr[i,j+1]))+w3*(sg(arr[i-1,j+1])+sg(arr[i+1,j+1]))))
            if i-1>=0 and j-1>=0 and i+1>=arr.shape[0] and j+1<arr.shape[1]:#rightmost column
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i-1,j]+arr[i,j-1]+arr[i,j+1])+w3*(arr[i-1,j-1]+arr[i-1,j+1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i-1,j])+sg(arr[i,j-1])+sg(arr[i,j+1]))+w3*(sg(arr[i-1,j-1])+sg(arr[i-1,j+1]))))
            if i-1>=0 and j-1>=0 and i+1<arr.shape[0] and j+1>=arr.shape[1]:#bottom row
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i-1,j]+arr[i,j-1]+arr[i+1,j])+w3*(arr[i-1,j-1]+arr[i+1,j-1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i-1,j])+sg(arr[i,j-1])+sg(arr[i+1,j]))+w3*(sg(arr[i-1,j-1])+sg(arr[i+1,j-1]))))
            if i-1<0 and j-1<0 and i+1<arr.shape[0] and j+1<arr.shape[1]:#up left corner
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i+1,j]+arr[i,j+1])+w3*(arr[i+1,j+1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i+1,j])+sg(arr[i,j+1]))+w3*(sg(arr[i+1,j+1]))))
            if i-1>=0 and j-1<0 and i+1>=arr.shape[0] and j+1<arr.shape[1]:#up right corner
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i-1,j]+arr[i,j+1])+w3*(arr[i-1,j+1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i-1,j])+sg(arr[i,j+1]))+w3*(sg(arr[i-1,j+1]))))
            if i-1<0 and j-1>=0 and i+1<arr.shape[0] and j+1>=arr.shape[1]:#down left corner
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i,j-1]+arr[i+1,j])+w3*(arr[i+1,j-1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i,j-1])+sg(arr[i+1,j]))+w3*(sg(arr[i+1,j-1]))))
            if i-1>=0 and j-1>=0 and i+1>=arr.shape[0] and j+1>=arr.shape[1]:#down right corner
                smooth_arr.append((w1*arr[i,j]+w2*(arr[i-1,j]+arr[i,j-1])+w3*(arr[i-1,j-1]))/(w1*sg(arr[i,j])+w2*(sg(arr[i-1,j])+sg(arr[i,j-1]))+w3*(sg(arr[i-1,j-1]))))
            
    return np.reshape(smooth_arr,(arr.shape[0],arr.shape[1])) 

"""Apply the smoothing function multiple times."""
def smooth(arr,n):
    
    sarr = arr
    
    for i in range(n):
        sarr = frame_smooth(sarr)
    
    return sarr
    
"""Checks if a point is inside a given contour."""
def check_in(i,j, contour):
    
    inside = False # Initialize the state
    inter_count = 0
    double_count = 0
    single_count = 0
    x0 = i # Define coordinates in terms of x,y to improve readability
    y0 = j
    
    for k in range(contour.shape[0]-1):
                # Define the end points of each line segment of the polygon
        x1 = contour[k,1]
        y1 = contour[k,0]
        x2 = contour[k+1,1]
        y2 = contour[k+1,0]
        
        #print(x1,y1,'\n',x2,y2)
        if y0<max(y1,y2) and y0>min(y1,y2):
            x_i = (y0*x2-y0*x1-y1*x2+y2*x1)/(y2-y1)
            if x_i>=x0:
                single_count = single_count+1
        if (y0==max(y1,y2) or y0==min(y1,y2)) and x0<=min(x1,x2):
            double_count = double_count+1
        
        inter_count = single_count+(double_count/2)
    
    if inter_count % 2 != 0: #If number of intersections is odd then point lies inside contour
        inside = True
    
    return inside

"""Takes a contour matrix and returns the boundaries of a box that enscribes it so that we can narrow the search for the peak."""
def box_contour(arr,C):
    
    x_max = np.amax(C[:,1])
    
    if ~np.isnan(x_max):
        if math.ceil(x_max)<=arr.shape[0]:
            x_max = math.ceil(x_max)
        else:
            x_max = math.floor(x_max)
    else:
        return "Nan encountered"
    
    x_min = np.amin(C[:,1])
    if ~np.isnan(x_min):
        if math.floor(x_min)>=0:
            x_min = math.floor(x_min)
        else:
            x_min = 0
    else:
        return "Nan encountered"
    
    y_max = np.amax(C[:,0])
    if ~np.isnan(y_max):
        if math.ceil(y_max)<=arr.shape[1]:
            y_max = math.ceil(y_max)
        else:
            y_max = math.floor(y_max)
    else:
        return "Nan encountered"
    
    y_min = np.amin(C[:,0])
    if ~np.isnan(y_min):
        if math.floor(y_min)>=0:
            y_min = math.floor(y_min)
        else:
            y_min = 0
    else:
        return "Nan encountered"
    
    if ~np.isnan(x_min) and ~np.isnan(x_max) and ~np.isnan(y_min) and ~np.isnan(y_max):
        return x_min,x_max,y_min,y_max
    else:
        return "Nan encountered"
    
"""Takes a contour matrix and finds its peak."""
#### Find also the distance to the FWHM
def find_peak(arr,C):
    
    temp = []
    
    if type(box_contour(arr,C)) == str:
        return (0,0,0,0)
    else:
        x_min,x_max,y_min,y_max = box_contour(arr,C)
        level = arr[int(C[0,1]),int(C[0,0])]
        for i in range(int(x_min),int(x_max)+1):
            for j in range(int(y_min),int(y_max)+1):
                if check_in(i,j,C) == True:
                    temp.append((arr[i,j],i,j))
        
        stemp = sorted(temp, key = lambda x:x[0], reverse = True)
        
        return stemp[0]

"""Takes an array and finds its contours between a minimum and a maximum value, in increments of 0.01."""
def find_contours(arr, c_min, c_max):
    
    arr_T = arr.T
    contour_list = []
    value_list = np.arange(c_min,c_max,0.01)
    
    for value in value_list:
        contours = measure.find_contours(arr.T,value)
        contour_list.append(contours)
    
    return contour_list    

"""Finds all the peaks and the contours within a range of a given array. Removes duplicate peaks."""
def peaks(arr,c_min,c_max):
    
    peak_list = []
    C_list = find_contours(arr,c_min,c_max)
    
    for contours in C_list:
        for contour in contours:
            peak_list.append(find_peak(arr,contour))
    
    return list(set(peak_list)), C_list

"""Rejects contours that enclose multiple peaks."""
def double_peak_rejection(peak_list, C_list):
    
    contour_list = []
    
    for contours in C_list:
        for contour in contours:
            counter = 0
            for peak in peak_list:
                if check_in(peak[1],peak[2],contour) == True:
                    counter = counter +1
            if counter<=1:
                contour_list.append(contour)
            else:
                pass
    
    return contour_list

"""Implements Shoelace formula to find the area of a contour. Taken from stackexchange."""
def PolyArea(x,y):
    
    return 0.5*np.abs(np.dot((x*unitR)+Rmin,np.roll((y*unitZ)+Zmin,1))-np.dot((y*unitZ)+Zmin,np.roll((x*unitR)+Rmin,1)))*10000.0

"""Takes each peak, finds all contours that share it and keeps only the one with the largest area."""
def contour_reduction(peak_list, C_list):
    
    contour_list = []
    pc_list=[]
    
    for peak in peak_list:
        for i,contour in enumerate(C_list):
            if check_in(peak[1],peak[2],contour) == True:
                area = PolyArea(contour[:,1],contour[:,0])
                pc_list.append((i,area))
        spc = sorted(pc_list, key = lambda x:x[1], reverse = True) 
        if len(spc)>0:
            contour_list.append(spc[0])
        pc_list=[]
    
    return contour_list

"""Functions for fitting ellipse on contours. Taken from http://nicky.vanforeest.com. Fitzgibbon's algorithm."""
def fitEllipse(x,y):
    
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    determinant = np.linalg.det(S)
    
    if determinant == 0.0:
        
        return "Singular Matrix"
    
    else:
        C = np.zeros([6,6])
        C[5,5] = 0; C[0,2] = C[2,0] = -2 ; C[1,1] = 1
        E, V =  eig(np.dot(inv(S), C))
        n = np.argmax(np.abs(E))
        a = V[:,n]
    
    return a

def ellipse_center(a):
    
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    
    return np.real(np.array([x0,y0]))

def ellipse_axis_length(a):
    
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    
    return np.real(np.array([res1, res2]))


def ellipse_angle_of_rotation(a):
    
    major, minor = ellipse_axis_length(a)
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    
    if major > minor:
        if b == 0:
            if a > c:
                return 0
            else:
                return np.pi/2
        else:
            if a > c:
                return np.arctan(2*b/(a-c))/2
            else:
                return np.pi + np.arctan(2*b/(a-c))/2
    else:
        if b == 0:
            if a > c:
                return np.pi/2
            else:
                return 0
        else:
            if a > c:
                return np.pi/2 + np.arctan(2*b/(a-c))/2
            else:
                return np.pi/2 + np.arctan(2*b/(a-c))/2
        

def draw_ellipse(x,y,a,b,t_rot):
    
    t = np.linspace(0, 2*math.pi, 100)
    Ell = np.array([a*np.cos(t) , b*np.sin(t)])  
     # x,y removed to keep the same center location
    R_rot = np.array([[math.cos(t_rot) , -math.sin(t_rot)],[math.sin(t_rot) , math.cos(t_rot)]])  
     # 2-D rotation matrix

    Ell_rot = np.zeros((2,Ell.shape[1]))
    
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
    
    return Ell_rot

"""Distance RMS error that shows the goodness of fit between a contour and its ellipse."""
def goodness_of_fit(x_cont,y_cont,x_el,y_el):
    
    temp = []
    error = []
    
    for i in range(len(x_cont)):
        for j in range(len(x_el)):
            temp.append(np.power(x_cont[i]-x_el[j],2)+np.power(y_cont[i]-y_el[j],2))
        error.append(min(temp))
        temp = []
    
    return np.sqrt(sum(error))

"""Takes the x,y coordinates of a contour and returns the parameters and coordinates of a fitted ellipse."""
def ellipse(x,y):
    
    el = fitEllipse(x,y)
    
    if el == "Singular Matrix":
        return "No ellipse"
    else:
        x_0, y_0 = ellipse_center(el)
        t_rot_0 = ellipse_angle_of_rotation(el)
        a_0, b_0 = ellipse_axis_length(el)
        Ell_rot = draw_ellipse(x_0,y_0,a_0,b_0,t_rot_0)
        return Ell_rot, x_0, y_0, t_rot_0, a_0, b_0

"""First find the arr_avg_glob to give as input to the db_creation function"""    
def large_array_avg(arr_glob):
    
    global arr_avg_glob
    arr_tor_glob = arr_glob.mean(1)
    arr_avg_glob = arr_tor_glob.mean(0)
    
    return arr_avg_glob

def prepare_data(arr_mid, arr_glob, arr_avg_glob, time):
    
    arr_tor = arr_mid.mean(1)
    arr_avg = arr_tor.mean(0)
    new_arr = arr_mid[time,0,:,:]
    for_arr = arr_glob[time,1,:,:]
    bac_arr = arr_glob[time,-1,:,:]
    work = (new_arr-arr_avg)/arr_avg # Definition of perturbation
    for_work = (for_arr-arr_avg_glob)/arr_avg_glob
    bac_work = (bac_arr-arr_avg_glob)/arr_avg_glob
    work = np.nan_to_num(work) # Avoid missing data/Nans
    for_work = np.nan_to_num(for_work)
    bac_work = np.nan_to_num(bac_work)
    sw = smooth(work,3) # Smoothing 
    for_sw = smooth(for_work,3)
    bac_sw = smooth(bac_work,3)
    sw = np.swapaxes(sw,0,1) # Watch out initial form
    for_sw = np.swapaxes(for_sw,0,1)
    bac_sw = np.swapaxes(bac_sw,0,1)
    sw = np.nan_to_num(sw)
    for_sw = np.nan_to_num(for_sw)
    bac_sw = np.nan_to_num(bac_sw)
    sw = field_line_smooth(sw,for_sw,bac_sw) # Field-line smoothing
    
    return sw

"""Takes a full array in the form [time, plane, Z, R] and a time and gives the n-times smoothed perturbed array in the form [R,Z] without doing the field-line smoothing. Faster, in order to be used for creating animations.
"""
def prepare_data2(arr_mid, plane, time):
    
    arr_tor = arr_mid.mean(1)
    arr_avg = arr_tor.mean(0)
    new_arr = arr_mid[time,plane,:,:]
    work = (new_arr-arr_avg)/arr_avg # Definition of perturbation
    sw = smooth(work,3) # Smoothing 
    sw = np.swapaxes(sw,0,1) # Watch out initial form
    
    return sw


"""Takes an [R,Z] array and returns the interpolation function of it."""
def array_interp(sw):
    
    x =  np.arange(0,sw.shape[0],1)
    y = np.arange(0,sw.shape[1],1)
    X,Y = np.meshgrid(x,y)
    nsw = np.nan_to_num(sw)
    f = interp2d(x, y, nsw[X,Y], kind='linear')
    
    return f
    
"""Takes a frame and a range of contour values and applies all the analysis. Returns four lists of tuples, whose length is equal to the number of blobs in that frame. So the first element of each list corresponds to the blob number. The contents of those lists are:

final_contour_list = [(contour)[y,x]]
final_peak_list = [(peak_value,peak_x,peak_y)]
level_list = [contour_level]
area_list = [area]
ellipse_list = [(el_rot[x,y], x, y, theta, a, b, goodness_of_fit)]

In this notation, to find the rise of blob i, we do: rise[i] = final_peak_list[i][0] - level_list[i]
"""
def analyze_frame(arr_mid, arr_glob, arr_avg_glob, time, bottom, top):
    
    sw = prepare_data(arr_mid, arr_glob, arr_avg_glob, time)
    f = array_interp(sw)
    peak_list, C_list = peaks(sw, bottom, top)
    contour_list = double_peak_rejection(peak_list,C_list)
    ncl = contour_reduction(peak_list, contour_list)
    ellipse_list = []
    level_list = []
    area_list = []
    final_contour_list = []
    final_peak_list = []
    
    for i in range(len(ncl)):
        contour_ind = ncl[i][0]
        contour = contour_list[contour_ind]
        final_contour_list.append(contour)
        level = f(contour[0,1],(contour[0,0]))
        level_list.append(level)
        area = PolyArea(contour[:,1],contour[:,0])
        area_list.append(area)    
        if np.isnan(contour).any() == False:
            if ellipse(contour[:,1]*unitR + Rmin,contour[:,0]*unitZ + Zmin) == "No ellipse":
                el_rot, x, y, t, a, b, good = (np.nan,np.nan,np.nan,np.nan,np.nan,np.nan, np.nan)
            else:
                el_rot, x, y, t, a, b = ellipse(contour[:,1]*unitR + Rmin,contour[:,0]*unitZ + Zmin)
                good = goodness_of_fit(contour[:,1],contour[:,0],x+el_rot[0,:], y+el_rot[1,:])
                ellipse_list.append((el_rot, x, y, t, a, b, good))
        else: 
            el_rot, x, y, t, a, b, good = (np.nan,np.nan,np.nan,np.nan,np.nan,np.nan, np.nan)
            ellipse_list.append((el_rot, x, y, t, a, b, good))
        for peak in peak_list:
            if check_in(peak[1],peak[2],contour) == True:
                final_peak_list.append(peak)       

    return final_contour_list, final_peak_list, level_list, ellipse_list, area_list

def analyze_frame2(arr_mid, plane, time, bottom, top):
    
    sw = prepare_data2(arr_mid, plane, time)
    f = array_interp(sw)
    peak_list, C_list = peaks(sw, bottom, top)
    contour_list = double_peak_rejection(peak_list,C_list)
    ncl = contour_reduction(peak_list, contour_list)
    ellipse_list = []
    level_list = []
    area_list = []
    final_contour_list = []
    final_peak_list = []
    
    for i in range(len(ncl)):
        contour_ind = ncl[i][0]
        contour = contour_list[contour_ind]
        final_contour_list.append(contour)
        level = f(contour[0,1],(contour[0,0]))
        level_list.append(level)
        area = PolyArea(contour[:,1],contour[:,0])
        area_list.append(area)    
        if np.isnan(contour).any() == False:
            if ellipse(contour[:,1],contour[:,0]) == "No ellipse":
                el_rot, x, y, t, a, b, good = (np.nan,np.nan,np.nan,np.nan,np.nan,np.nan, np.nan)
            else:
                el_rot, x, y, t, a, b = ellipse(contour[:,1],contour[:,0])
                good = goodness_of_fit(contour[:,1],contour[:,0],x+el_rot[0,:], y+el_rot[1,:])
                ellipse_list.append((el_rot, x, y, t, a, b, good))
        else: 
            el_rot, x, y, t, a, b, good = (np.nan,np.nan,np.nan,np.nan,np.nan,np.nan, np.nan)
            ellipse_list.append((el_rot, x, y, t, a, b, good))
        for peak in peak_list:
            if check_in(peak[1],peak[2],contour) == True:
                final_peak_list.append(peak)       

    return final_contour_list, final_peak_list, level_list, ellipse_list, area_list


"""Function that draws the frame with all the contours, peaks and ellipses."""
def draw_frame(sw, C_l, P_l, E_l):
    
    x =  np.arange(0,sw.shape[0],1)
    y = np.arange(0,sw.shape[1],1)
    X,Y = np.meshgrid(x,y)
    fig, ax = plt.subplots(figsize=(20,10))
    ax.pcolor(X,Y,sw[X,Y])
    
    for contour in C_l:
        ax.plot(contour[:, 1],contour[:, 0], linewidth=2)
    for j in range(len(E_l)):
        if ~np.isnan(E_l[j][0]).any():
            ax.plot(E_l[j][1]+E_l[j][0][0], E_l[j][2]+E_l[j][0][1],linewidth = 2, color = 'k',linestyle = "--")
        else:
            pass
    for peak in P_l:
        ax.plot(peak[1],peak[2],"b", marker='x',markersize = 6)    
    
    # Setting the axes labels
    x_lbl,y_lbl = plot_labels(Rmin,Rmax,Zmin,Zmax,5,5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.set_xticklabels(x_lbl)
    ax.set_yticklabels(y_lbl)
    
    # Drawing the separatrix (Need the right units)
    #R_sep, Z_sep = separatrix_transform(R_sep_list,Z_sep_list)
    #ax.plot(R_sep,Z_sep,'k--')#ti262 midplane
    ax.set_xlabel('R(m)')
    ax.set_ylabel("Z(m)")
        
"""Function that draws adjacent in time frames, side by side."""        
def draw_frames(arr_all, time1, time2, bottom, top):
    
    C_l_1, P_l_1, L_l_1, E_l_1, A_l_1 = analyze_frame(arr_all, time1, bottom, top)
    C_l_2, P_l_2, L_l_2, E_l_2, A_l_2 = analyze_frame(arr_all, time2, bottom, top)
    sw_1 = prepare_data(arr_all,time1)
    sw_2 = prepare_data(arr_all,time2)
    labels1 = ['blob{0}'.format(i) for i in range(len(E_l_1))]
    labels2 = ['blob{0}'.format(i) for i in range(len(E_l_2))]
    x =  np.arange(0,sw_1.shape[0],1)
    y = np.arange(0,sw_1.shape[1],1)
    X,Y = np.meshgrid(x,y)
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(20,10))
    ax1.pcolor(X,Y,sw_1[X,Y])
    
    for contour in C_l_1:
        ax1.plot(contour[:, 1],contour[:, 0], linewidth=2)
    for j in range(len(E_l_1)):
        if ~np.isnan(E_l_1[j][0]).any():
            ax1.plot(E_l_1[j][1]+E_l_1[j][0][0], E_l_1[j][2]+E_l_1[j][0][1],linewidth = 2, color = 'k',linestyle = "--")
            ax1.annotate(labels1[j], xy=(E_l_1[j][1],E_l_1[j][2]),xytext=(-20,20),textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        else:
            pass
    
    for peak in P_l_1:
        ax1.plot(peak[1],peak[2],"b", marker='x',markersize = 6)
        
    ax2.pcolor(X,Y,sw_2[X,Y])
    
    for contour in C_l_2:
        ax2.plot(contour[:, 1],contour[:, 0], linewidth=2)
    
    for j in range(len(E_l_2)):
        if ~np.isnan(E_l_2[j][0]).any():
            ax2.plot(E_l_2[j][1]+E_l_2[j][0][0], E_l_2[j][2]+E_l_2[j][0][1],linewidth = 2, color = 'k',linestyle = "--")
            ax2.annotate(labels2[j], xy=(E_l_2[j][1],E_l_2[j][2]),xytext=(-20,20),textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        else:
            pass
    
    for peak in P_l_2:
        ax2.plot(peak[1],peak[2],"b", marker='x',markersize = 6)
    
    # Setting the axes labels
    x_lbl,y_lbl = plot_labels(Rmin,Rmax,Zmin,Zmax,5,5)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax1.set_xticklabels(x_lbl)
    ax1.set_yticklabels(y_lbl)
    ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.set_xticklabels(x_lbl)
    ax2.set_yticklabels(y_lbl)
    
    # Drawing the separatrix (Need the right units)
    #R_sep, Z_sep = separatrix_transform(R_sep_list,Z_sep_list)
    #ax1.plot(R_sep,Z_sep,'k--')#ti262 midplane
    ax1.set_xlabel('R(m)')
    ax1.set_ylabel("Z(m)")
    #ax2.plot(R_sep,Z_sep,'k--')
    ax2.set_xlabel('R(m)')
    ax2.set_ylabel("Z(m)")
        
    
"""Database functions"""
"""Function that creates a table of blob values within the database at a particular time."""        
def table_creation(arr_mid, arr_glob, arr_avg_glob, time, bottom, top):
    
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("SELECT name FROM sqlite_master WHERE type = 'table' AND name = 't_"+str(time)+"'")
    a = c.fetchall()
    
    if not a:
        for plane in range(Nplanes):
        #C_l, P_l, L_l, E_l, A_l = analyze_frame(arr_mid, arr_glob, arr_avg_glob, time, bottom, top)
            Cl, P_l, L_l, E_l, A_l = analyze_frame2(arr_mid, plane, time, bottom, top)
            c.execute('''CREATE TABLE IF NOT EXISTS t_'''+str(time)+'''(peak_value float, peak_x float, peak_y float, area float, level float, el_x float, el_y float, tilt float, maj_ax float, min_ax float, goodness float, plane)''')
            for i in range(len(E_l)):
                (c.execute('''INSERT INTO t_'''+str(time)+''' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',(float(P_l[i][0]), float(P_l[i][1]), float(P_l[i][2]), float(A_l[i]), float(L_l[i][0]), float(E_l[i][1]), float(E_l[i][2]), float(E_l[i][3]), float(E_l[i][4]), float(E_l[i][5]), float(E_l[i][6]), int(plane),)))
    else:
        pass
    conn.commit()
    conn.close()

"""Function that creates tables of blob values within a range of times."""    
def db_creation(arr_mid, arr_glob, arr_avg_glob, t_start, t_end, bottom, top):
    
    for time in range(t_start, t_end):
        table_creation(arr_mid, arr_glob, arr_avg_glob, time, bottom, top)

"""Function that adds new columns to the db that contain the rescaled major/minor axes."""    
#def rescaling(t_start,t_end):
#    
#    conn = sqlite3.connect('/global/cscratch1/sd/giannos/ti255_arrays/db2')
#    c = conn.cursor()
#
#    for time in range(t_start,t_end):
#        c.execute('ALTER TABLE t_'+str(time)+ ' ADD COLUMN maj_p')
#        c.execute('ALTER TABLE t_'+str(time)+ ' ADD COLUMN min_p')
#        c.execute('ALTER TABLE t_'+str(time)+ ' ADD COLUMN tilt_p')
#        blobs = []
#        for row in c.execute('SELECT ROWID, * FROM t_'+str(time)):
#            blobs.append(row)
#        for blob in blobs:
#            if ~pd.isnull([blob[8],blob[9],blob[10]]).any():
#                a = blob[9]
#                b = blob[10]
#                t = blob[8]
#            
#                kx = (1./np.power(unitR,2))*(np.power(np.cos(t)/a,2)+np.power(np.sin(t)/b,2))
#                ky = (1./np.power(unitZ,2))*(np.power(np.sin(t)/a,2)+np.power(np.cos(t)/b,2))
#                kxy = ((np.cos(t)*np.sin(t))/(unitR*unitZ))*(1./np.power(a,2) - 1./np.power(b,2))
#           
#                tp = (1./2.)*np.arctan(2*kxy/(kx-ky))
#                tmp_a = kx + kx + 2*kxy/np.sin(2*tp)
#                tmp_b = kx + kx - 2*kxy/np.cos(2*tp)
#                ap = np.sqrt(2./tmp_a)
#                bp = np.sqrt(2./tmp_b)
#            
#                app = max(ap,bp)
#                bpp = min(ap,bp)
#            
#                if ap == app:
#                    tpp = tp
#                else:
#                    tpp = tp + np.pi/2.
#                
#                c.execute('''UPDATE t_'''+str(time)+''' SET maj_p = ?, min_p = ?, tilt_p = ? WHERE ROWID = #'''+str(blob[0]),(app, bpp, tpp))
#            else:
#                c.execute('''UPDATE t_'''+str(time)+''' SET maj_p = ?, min_p = ?, tilt_p = ? WHERE ROWID = #'''+str(blob[0]),(None, None, None))
#    
#    conn.commit()
#    conn.close()

"""Function that draws a frame using the right database table."""
def db_draw_frame(arr_all,plane,time):
    sw = prepare_data2(arr_all,plane,time)
    x =  np.arange(0,sw.shape[0],1)
    y = np.arange(0,sw.shape[1],1)
    X,Y = np.meshgrid(x,y)

    conn = sqlite3.connect(db_name)
    c = conn.cursor()

    rows = []
    for row in c.execute("SELECT ROWID,* FROM t_"+str(time)):
        rows.append(row)

    fig, ax = plt.subplots(figsize=(20,10))
    ax.pcolor(X,Y,sw[X,Y])
    for row in rows:
        if row[4]>=1.0:
            if ~pd.isnull([row[6],row[7],row[9],row[10],row[8]]).any():
                Ell_rot = draw_ellipse(row[6],row[7],row[9],row[10],row[8])
                ax.plot(row[6]+Ell_rot[0,:],row[7]+Ell_rot[1,:],linewidth = 2, color = 'k',linestyle = "--")
                ax.plot(row[2],row[3],"b", marker='x',markersize = 6)
                ax.annotate('blob_'+str(row[0]), xy=(row[6],row[7]),xytext=(-20,20),textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

    # Setting the axes labels
    x_lbl,y_lbl = plot_labels(Rmin,Rmax,Zmin,Zmax,5,5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.set_xticklabels(x_lbl)
    ax.set_yticklabels(y_lbl)

    #R_sep, Z_sep = separatrix_transform(R_sep_list,Z_sep_list)
    #ax.plot(R_sep,Z_sep,'k--')#ti262 midplane

    ax.set_xlabel('R(m)')
    ax.set_ylabel("Z(m)")

    conn.close()

"""Function that draws two dataframes from the database."""    
def db_draw_frames(arr_all,plane,time1,time2):
    
    sw1 = prepare_data2(arr_all,plane,time1)
    sw2 = prepare_data2(arr_all,plane,time2)
    x =  np.arange(0,sw1.shape[0],1)
    y = np.arange(0,sw1.shape[1],1)
    X,Y = np.meshgrid(x,y)
    
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    rows1 = []
    for row in c.execute("SELECT ROWID,* FROM t_"+str(time1)):
        rows1.append(row)
    
    rows2 = []
    for row in c.execute("SELECT ROWID,* FROM t_"+str(time2)):
        rows2.append(row)
    
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(20,10))
    ax1.pcolor(X,Y,sw1[X,Y])
    
    for row in rows1:
        if row[4]>=0.2:
            if ~pd.isnull([row[6],row[7],row[9],row[10],row[8]]).any():
                Ell_rot = draw_ellipse(row[6],row[7],row[9],row[10],row[8])
                ax1.plot(row[6]+Ell_rot[0,:],row[7]+Ell_rot[1,:],linewidth = 2, color = 'k',linestyle = "-")
                ax1.plot(row[2],row[3],"b", marker='x',markersize = 6)
                ax1.annotate('b'+str(row[0])+'_'+"%.2f"%(row[1]), xy=(row[6],row[7]),xytext=(-20,20),textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
                
    ax2.pcolor(X,Y,sw2[X,Y])
    for row in rows2:
        if row[4]>=0.2:
            if ~pd.isnull([row[6],row[7],row[9],row[10],row[8]]).any():
                Ell_rot = draw_ellipse(row[6],row[7],row[9],row[10],row[8])
                ax2.plot(row[6]+Ell_rot[0,:],row[7]+Ell_rot[1,:],linewidth = 2, color = 'k',linestyle = "-")
                ax2.plot(row[2],row[3],"b", marker='x',markersize = 6)
                ax2.annotate('b'+str(row[0])+'_'+"%.2f"%(row[1]), xy=(row[6],row[7]),xytext=(-20,20),textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    
    # Setting the axes labels
    x_lbl,y_lbl = plot_labels(Rmin,Rmax,Zmin,Zmax,5,5)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax1.set_xticklabels(x_lbl)
    ax1.set_yticklabels(y_lbl)
    ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.set_xticklabels(x_lbl)
    ax2.set_yticklabels(y_lbl)
    
    # Drawing the separatrix (Need the right units)
    #R_sep, Z_sep = separatrix_transform(R_sep_list,Z_sep_list)
    #ax1.plot(R_sep,Z_sep,'k--')#ti262 midplane
    
    ax1.set_xlabel('R(m)')
    ax1.set_ylabel("Z(m)")
    #ax2.plot(R_sep,Z_sep,'k--')
    ax2.set_xlabel('R(m)')
    ax2.set_ylabel("Z(m)")
                
    conn.close()
    
"""Returns only the exponential part of the normal distribution without the normalization."""
def Normal(x,mu,sigma):
    
    var = sigma*sigma
    
    return np.exp(-np.power(x-mu,2)/(2*var))

"""Returns the pdf of radial blob velocity."""    
def pdf_Rad(x):
    
    global mu_rad
    # In units of m/s
    mu_rad = 300
    sigma = 1000 
    # Choose normal distribution scipy.stats.norm.pdf(x,mu,sigma)
    
    return Normal(x,mu_rad,sigma)

pdf_Rad(0)

"""Returns the pdf of poloidal blob velocity."""    
def pdf_Pol(x):
    
    global mu_pol
    # In units of m/s
    mu_pol = 4000 
    sigma = 1000 
    # Choose normal distribution
    
    return Normal(x,mu_pol,sigma)

pdf_Pol(0)

"""Returns the pdf of the area difference."""
def pdf_A(A1,A2):
    
    dA = np.abs(A1-A2)
    sigma = np.sqrt(A1*A2)
    
    return Normal(dA,0.,sigma)

"""Returns the pdf of the height difference."""
def pdf_h(h1,h2):
    
    dh = np.abs(h1-h2)
    sigma = np.sqrt(h1*h2)/3.
    
    return Normal(dh,0,sigma)

"""Returns the normalization factor coming from a perfect match."""
def score_calibration():
    
    global mu_rad, mu_pol
    calibration = pdf_Rad(mu_rad)*pdf_Pol(mu_pol)*pdf_A(10,10)
    
    return calibration

C = score_calibration()

"""Score function: Takes two blobs (lists) from adjacent frames and calculates the probability that they are the same blob."""
def score(blobA,blobB):

    vR = ((blobB[6]-blobA[6])*unitR)/dt
    vZ = ((blobB[7]-blobA[7])*unitZ)/dt
    Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2.
    Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2.
    Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2.
    
    V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
    V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
    score = pdf_Rad(np.abs(V_psi))*pdf_Pol(np.abs(V_pol))*pdf_A(blobA[4],blobB[4])*(1/C)
    #print('blob',blobA[0],'-','blob',blobB[0],'V_R:',V_psi,'score:',score, 'rad:',pdf_Rad(np.abs(vR)),'pol:',pdf_Pol(np.abs(vZ)),'area:',pdf_A(blobA[4],blobB[4]),'norm:',1/C)
    
    return score

"""Takes two adjacent frames and computes the score matrix between their individual blobs. Also returns the area dictionary of the first frame.
"""
def score_matrix(time1,time2,plane):

    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    blobs_1 = []
    for row in c.execute("SELECT ROWID,* FROM t_"+str(time1)):
        if row[12] == plane:
            blobs_1.append(row)
    blobs_2 = []
    for row in c.execute("SELECT ROWID,* FROM t_"+str(time2)):
        if row[12] == plane:
            blobs_2.append(row)
    
    SM = [0.0]
    for blob in blobs_2:
        if type(blob[4]) == type(blob[6]) == type(blob[7]) == type(blob[8]) == type(blob[9]) == type(blob[10]) == float and blob[4]>=5e-5:
            SM.append(blob[0])
    
    SM = np.reshape(np.array(SM),(-1,1))
    np.set_printoptions(formatter={'float': '{: 0.2f}'.format})
    col = []
    areas = {}
    for blob1 in blobs_1:
        if type(blob1[4]) == type(blob1[6]) == type(blob1[7]) == type(blob1[8]) == type(blob1[9]) == type(blob1[10]) == float and blob1[4]>=5e-5:
            
            col.append(blob1[0])
            areas[blob1[0]] = blob1[4]
            for blob2 in blobs_2:
                if type(blob2[4]) == type(blob2[6]) == type(blob2[7]) == type(blob2[8]) == type(blob2[9]) == type(blob2[10]) == float and blob2[4]>=5e-5:
                    col.append(score(blob1,blob2))
            SM = np.insert(SM,SM.shape[1],col,axis=1)
        col = []
        
    conn.close()
    #for row in SM:
        #print(row)
    #for line in SM:
        #print(np.array_str(line,precision=9))
    
    return SM,areas

"""Function that searches the score array and creates forward and backward dictionaries for blob relations."""
def search_array(S):
    
    blob_dict = {}
    
    for i in range(1,S.shape[1]):
        key = int(S[0][i])
        blob_dict.setdefault(key,[])
        for j in range(1,S.shape[0]):
            if S[j][i]>= 0.01:
                blob_dict[key].append(int(S[j][0])) 
    
    reverse_dict = {}
    
    for i in range(1,S.shape[0]):
        key = int(S[i][0])
        reverse_dict.setdefault(key,[])
        for j in range(1,S.shape[1]):
            if S[i][j]>= 0.01:
                reverse_dict[key].append(int(S[0][j]))
        
    return blob_dict, reverse_dict

"""Function that checks for mergers. If it detects one, discontinues all merged blobs except the one with the largest area."""
def mergers(blob_dict,reverse_dict,area_dict):
    
    for merged_blob in reverse_dict:
        if len(reverse_dict[merged_blob])>1:
            area_list = []
            for init_blob in reverse_dict[merged_blob]:
                area_list.append(area_dict[init_blob])
            blob_with_max_area = reverse_dict[merged_blob][area_list.index(max(area_list))]
            for init_blob in reverse_dict[merged_blob]:
                if init_blob != blob_with_max_area:
                    blob_dict[init_blob] = [0]
    
    return blob_dict    

"""Function that creates time histories of blobs from the database. Each blob appears as a row in the matrix, the index is the time point and the value is the number of the plot in the database table of that time point. Zero denotes the absence of the blob.
"""
def blob_history(t_start,t_end):

    time_window = (t_end-t_start)+1
    history = np.asarray(np.zeros(time_window)) # Array initialization
    history = history[np.newaxis,:] # Making it a [1,time_window] array
    
    for time in range(t_start,t_end):
        for plane in range(Nplanes):
            S, area_dict = score_matrix(time,time+1,plane)
            blob_dict, reverse_dict = search_array(S)
            fixed_blob_dict = mergers(blob_dict, reverse_dict, area_dict)
        
            for blob in fixed_blob_dict:
                if len(fixed_blob_dict[blob]) == 0: # Case when blob ceizes to exist
                    for row in history:
                        if row[time-t_start] == blob:
                            row[time-t_start+1] = 0
                
                if len(fixed_blob_dict[blob]) == 1: # Case when blob is mapped to exactly one blob
                    exists = 0
                    for row in history: # Case when blob already exists
                        if row[time-t_start] == blob:
                            row[time-t_start+1] = fixed_blob_dict[blob][0]
                            exists = exists+1
                    if exists == 0: # Case when blob just appeared
                        init_zeros = list(np.zeros(time-t_start))
                        blob_values = [blob,fixed_blob_dict[blob][0]]
                        end_zeros = list(np.zeros(time_window-((time-t_start)+2)))
                        new_blob = np.asarray(init_zeros+blob_values+end_zeros)
                        new_blob = new_blob[np.newaxis,:]
                        history = np.vstack([history,new_blob])
                
                if len(fixed_blob_dict[blob]) > 1: # Case when blob splits in many blobs
                    for split_blob in fixed_blob_dict[blob]:
                        init_zeros = list(np.zeros(time-t_start))
                        blob_values = [blob,split_blob]
                        end_zeros = list(np.zeros(time_window-((time-t_start)+2)))
                        new_blob = np.asarray(init_zeros+blob_values+end_zeros)
                        new_blob = new_blob[np.newaxis,:]
                        history = np.vstack([history,new_blob])
    
    history = np.delete(history,(0),axis=0) # Deleting first row of zeros
    int_history = history.astype(int) # Turning all elements to integers
    
    return int_history 

"""Function that detects splits of blobs."""
def split_detection(H):
    
    splits = []
    
    for i in range(0,H.shape[0]):
        split_list = [i]
        ind = next((index for index,value in enumerate(H[i]) if value != 0), None)
        for j in range(0,H.shape[0]):
            if i != j and H[j][ind] == H[i][ind]:
                split_list.append(j)
        splits.append(split_list)
    
    return splits

"""Function that detects the longest persisting blobs."""
def longest_blob(H):
    
    nz = []
    for i in range(0,H.shape[0]):
        nz.append(np.count_nonzero(H[i]))
    first = nz.index(max(nz))
    first_len = max(nz)
    nz[first] = 0
    second = nz.index(max(nz))
    second_len = max(nz)
    nz[second] = 0
    third = nz.index(max(nz))
    third_len = max(nz)
    nz[third] = 0
    fourth = nz.index(max(nz))
    fourth_len = max(nz)
    
    return first, first_len, second, second_len, third, third_len, fourth, fourth_len

"""Animates the time history of a single blob."""
def blob_animation(arr_mid, blob_history_array, blob_number, t_offset, t_start, t_end):
    """Here, t_offset refers to the time point that the database starts. For example, if the database contains time steps 600-880, then t_offset=600. But if we only want a movie of 20 time steps, e.g., 720-740, we need t_start = 720, t_end = 740
    """
    conn = sqlite3.connect(db_name) # Connection to database
    c = conn.cursor()
    
    # Preparation of grid for frame
    sw = prepare_data2(arr_mid,0) # Use array without field-line smoothing to save time.
    x =  np.arange(0,sw.shape[0],1)
    y = np.arange(0,sw.shape[1],1)
    X,Y = np.meshgrid(x,y)
    
    # Picking the time history row of the blob
    blob_history = blob_history_array[blob_number]
    time_range = t_end-t_start
    offset2 = t_start-t_offset
    #time_range = len(blob_history_array[blob_number])
    
    fig,ax = plt.subplots(figsize=(8,8))
    ims=[]
    
    #for time in range(t_offset,t_offset+time_range):
    for time in range(t_start, t_end):
        sw = prepare_data2(arr_mid,time)
        im = ax.pcolor(X,Y,sw[X,Y],cmap='seismic')
        blob_id = blob_history[time-t_offset]
        if blob_id > 0: #if blob exists
            blob = c.execute('SELECT ROWID, * FROM t_'+str(time)+' WHERE ROWID='+str(blob_id)) # Pick blob from database
            for row in blob:
                Ell_rot = draw_ellipse(row[6],row[7],row[9],row[10],row[8])
                bl, = ax.plot(row[6]+Ell_rot[0,:],row[7]+Ell_rot[1,:],linewidth = 1, color = 'k',linestyle = "-")
        
        else: # Do nothing if blob doesnt exist
            bl = ax.annotate("", xy = (0.0,0.0), xycoords='axes fraction',fontsize=14) 
        
        text = 't = '+str(time)
        an = ax.annotate(text, xy=(0.8, 0.94), xycoords='axes fraction',fontsize=14)
        
        '''#################################################################
    ## Bug fix for Quad Contour set not having attribute 'set_visible', 'set_animated'
        def setvisible(self,vis):
            for c in self.collections: c.set_visible(vis)'/global/cscratch1/sd/giannos/ti255_arrays/db2'
        def setanimated(self,ani):
            for c in self.collections: c.set_animated(ani)
        im.set_visible = types.MethodType(setvisible,im)
        im.set_animated = types.MethodType(setanimated,im)
        im.axes = plt.gca()
        im2 = im.collections
        #im.figure=fig
          ####################################################################'''
            
        ims.append([im]+[an]+[bl])
    
    conn.close()
    # Setting the axes labels
    x_lbl,y_lbl = plot_labels(Rmin,Rmax,Zmin,Zmax,5,5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.set_xticklabels(x_lbl)
    ax.set_yticklabels(y_lbl)

    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)    
    art_ani = animation.ArtistAnimation(fig, ims, interval=1000, repeat_delay=1000,blit=False)
    
    # Drawing the separatrix (Need the right units)
    #R_sep, Z_sep = separatrix_transform(R_sep_list,Z_sep_list)
    #plt.plot(R_sep,Z_sep,'k--')#ti262 midplane
    plt.xlabel('R(m)')
    plt.ylabel("Z(m)")
    plt.show()
    
    return art_ani

"""Animates a list of blobs."""
def list_animation(arr_all,blob_history_array, blob_list, t_offset):
    
    conn = sqlite3.connect(db_name) # Connection to database
    c = conn.cursor()
    
    # Preparation of grid for frame
    sw = prepare_data(arr_all,0)
    x =  np.arange(0,sw.shape[0],1)
    y = np.arange(0,sw.shape[1],1)
    X,Y = np.meshgrid(x,y)
    
    time_range = len(blob_history_array[0])
    num_of_blobs = len(blob_list)    
    
    colours = list(mcolors.BASE_COLORS.keys())
    while len(colours)<= num_of_blobs:
        colours = colours + colours
    
    fig,ax = plt.subplots(figsize=(8,8))
    ims=[]
    
    for time in range(t_offset,t_offset+time_range):
        sw = prepare_data(arr_all,time)
        im = ax.pcolor(X,Y,sw[X,Y])
        
        imb = []
        for i,blob_position in enumerate(blob_list):
            blob_history = blob_history_array[blob_position]
            blob_id = blob_history[time-t_offset]
            if blob_id > 0: #if blob exists
                blob = c.execute('SELECT ROWID, * FROM t_'+str(time)+' WHERE ROWID='+str(blob_id)) # Pick blob from database
                for row in blob:
                    if  type(row[6]) == type(row[7]) == type(row[8]) == type(row[9]) == type(row[10]) == float:
                        Ell_rot = draw_ellipse(row[6],row[7],row[9],row[10],row[8])
                        bl, = ax.plot(row[6]+Ell_rot[0,:],row[7]+Ell_rot[1,:],linewidth = 1.5, color = colours[i],linestyle = "-")
                    else:
                        bl = ax.annotate("", xy = (0.0,0.0), xycoords='axes fraction',fontsize=14) 
            else: # Do nothing if blob doesnt exist
                bl = ax.annotate("", xy = (0.0,0.0), xycoords='axes fraction',fontsize=14) 
            imb = imb + [bl] 
        
        text = 't = '+str(time)
        an = ax.annotate(text, xy=(0.8, 0.94), xycoords='axes fraction',fontsize=14)
        
        '''#################################################################
    ## Bug fix for Quad Contour set not having attribute 'set_visible', 'set_animated'
        def setvisible(self,vis):
            for c in self.collections: c.set_visible(vis)
        def setanimated(self,ani):
            for c in self.collections: c.set_animated(ani)
        im.set_visible = types.MethodType(setvisible,im)
        im.set_animated = types.MethodType(setanimated,im)
        im.axes = plt.gca()
        im2 = im.collections
        #im.figure=fig
          ####################################################################'''
            
        ims.append([im]+[an]+imb)

    conn.close()
    
    # Setting the axes labels
    x_lbl,y_lbl = plot_labels(Rmin,Rmax,Zmin,Zmax,5,5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.set_xticklabels(x_lbl)
    ax.set_yticklabels(y_lbl)      
    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)    
    art_ani = animation.ArtistAnimation(fig, ims, interval=1000, repeat_delay=1000,blit=False)
    
    # Drawing the separatrix (Need the right units)
    R_sep, Z_sep = separatrix_transform(R_sep_list,Z_sep_list)
    plt.plot(R_sep,Z_sep,'k--') # ti262 midplane
    plt.xlabel('R(m)')
    plt.ylabel("Z(m)")
    plt.show()
    
    return art_ani

"""Functions that plot statistical distributions from the blobs database."""
"""Probability distribution of blob sizes."""
def area_dist(time_start,time_end):
    
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    areas = []
    count = 0
    for time in range(time_start,time_end):
        for row in c.execute("SELECT ROWID,* FROM t_"+str(time)):
            if ~pd.isnull(row[4]) and row[4] is not None and row[4]>5e-5:
                areas.append(row[4])
                count = count+1
            else:
                pass
    print("Blob count:", count)
    print("mean", np.mean(np.sqrt(areas)))
    plt.title(r"$Area^{\frac{1}{2}}$ Distribution")
    sns.distplot(np.sqrt(areas),hist=True,kde=True,color='blue')
    #plt.axvline(x=delta_star,linestyle='--',label=r'$\delta_{\star}$')
    #plt.axvline(x=rho_s,label=r'$\rho_s$')
    plt.xlabel('cm')
    plt.legend()
    plt.show()

def hwhm(d,level,peak):
    """Takes the diameter of the blob in meters and returns the HWHM in cm."""
    Ang1 = np.arcsin((level+1.)/(peak+1.))
    
    hwhm = (np.pi*d)/(3*(np.pi-2*Ang1))
    
    return hwhm*100

    
def proj_length(maj_r, tilt, Br, Bz, Bzz, Bpol, Bmag):
    """Takes blob information and returns the projected length to be used by the hwhm function. Returns in meters."""
    # Projection on R and Z axes in m
    d_r = np.cos(tilt)*maj_r
    d_z = np.sin(tilt)*maj_r
    
    # Projection on binormal direction
    pref = 1./(Bmag*Bpol)
    d_bn = pref*(-Bzz*Br*d_r - Bzz*Bz*d_z)
    
    return d_bn   

def size_dist_nrs(time_start,time_end):
    
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    maj = []
    proj_sizes1 = []
    proj_sizes2 = []
    proj_sizes3 = []
    sqrt_area = []
    counter = 0
    
    for time in range(time_start,time_end):
        for row in c.execute("SELECT ROWID,* FROM t_"+str(time)):
            try:
                if ~pd.isnull(row[4]) and row[4] is not None and row[4]>5e-5 and ~pd.isnull(row[1]) and row[1] is not None and ~pd.isnull(row[10]) and row[10] is not None and ~pd.isnull(row[5]) and row[5] is not None and ~pd.isnull(row[8]) and row[8] is not None and ~pd.isnull(row[9]) and row[9] is not None:
                    # Define blob characteristics.
                    x = int(row[2]) # Blobs peak coordinates
                    y = int(row[3])
                    peak = float(row[1])
                    level = float(row[5])
                    area = float(row[4])
                    tilt = float(row[8])
                    # Find major/minor axes.
                    MajR = float(row[9])
                    MinR = float(row[10])
                    if MinR > MajR:
                        temp = MajR
                        MajR = MinR
                        MinR = temp
                    # Evaluate the magnetic field components at the blob peak.
                    Br = BR[y][x] # Careful with the peculiar definition of the B-field components. 
                    Bz = BZ[y][x]
                    Bzz = Bzeta[y][x]
                    Bpol = Bp[y][x]
                    Bmag = np.sqrt(Br**2 + Bz**2 + Bzz**2)
                    
                    # Quantities to plot.
                    sqrt_area.append(np.sqrt(area))
                    # Reject unnaturally large blobs.
                    if MajR*100 < 0.5 and peak < 2:
                        counter = counter + 1 
                        maj.append(MajR*100)
                        d1 = np.sqrt(np.pi*MajR*MinR)
                        d2 = np.sqrt(area)/100
                        d3 = proj_length(MajR, tilt, Br, Bz, Bzz, Bpol, Bmag)
                        d_z = np.sin(tilt)*MajR
                        proj_sizes1.append(hwhm(d3,level,peak))
                        proj_sizes2.append(hwhm(d_z,level,peak))
                        proj_sizes3.append(hwhm(MajR,level,peak))
                    else:
                        pass
                else:
                    pass
            except IndexError:
                pass
    
    WL = (2*np.pi/k_perp)*100 # Wavelength of linear mode in cm
    HWHM = WL/6. #  half width at half maximum size of linear mode in cm
    maj = [x for x in maj if str(x)!='nan' and float(x)<10]
    proj_sizes1 = [x for x in proj_sizes1 if str(x)!='nan' and float(x)<3. and float(x)>0]
    
    #file = open("/global/cscratch1/sd/giannos/cmod_sizes(fig.10.b).txt",'a')
    #file.write("List of blob sizes"+"\n")
    #for i in range(len(proj_sizes1)):
        #file.write(str(proj_sizes1[i])+"\n")
    #file.close()    
    
    proj_sizes2 = [x for x in proj_sizes2 if str(x)!='nan' and float(x)<3. and float(x)>0]
    proj_sizes3 = [x for x in proj_sizes3 if str(x)!='nan' and float(x)<3. and float(x)>0]
    sqrt_area = [x for x in sqrt_area if str(x)!='nan']
    
    print("Number of blobs:", counter)
    plt.title(r"HWHM of Blobs Distribution") #D=binormal
    sns.distplot(proj_sizes1, hist=True, kde=True, color='blue')
    plt.axvline(x = delta_star, linestyle='--', label=r'$\delta_{\star}$')
    plt.axvline(x = rho_s, label=r'$\rho_s$')
    plt.axvline(x = HWHM, linestyle='-.', color='black', label=r'HWHM_LM')
    plt.xlabel('cm') 
    plt.legend()
    plt.show()
    plt.title(r"HWHM of Blobs Distribution, D=d_z")
    sns.distplot(proj_sizes2, hist=True, kde=True, color='blue')
    plt.axvline(x = delta_star, linestyle='--', label=r'$\delta_{\star}$')
    plt.axvline(x = rho_s, label=r'$\rho_s$')
    plt.axvline(x = HWHM, linestyle='-.', color='black', label=r'HWHM_LM')
    plt.xlabel('cm') 
    plt.legend()
    plt.show()
    plt.title(r"HWHM of Blobs Distribution, D=R_M")
    sns.distplot(proj_sizes3, hist=True, kde=True, color='blue')
    plt.axvline(x = delta_star, linestyle='--', label=r'$\delta_{\star}$')
    plt.axvline(x = rho_s, label=r'$\rho_s$')
    plt.axvline(x = HWHM, linestyle='-.', color='black', label=r'HWHM_LM')
    plt.xlabel('cm') 
    plt.legend()
    plt.show()
    plt.title(r"Major Radius Distribution")
    sns.distplot(maj, hist=True, kde=True, color='blue')
    plt.axvline(x = delta_star, linestyle='--', label=r'$\delta_{\star}$')
    plt.axvline(x = rho_s, label=r'$\rho_s$')
    plt.axvline(x = HWHM, linestyle='-.', color='black', label=r'HWHM_LM')
    plt.xlabel('cm')
    plt.legend()
    plt.show()
    plt.title(r"Sqrt Area Distribution")
    sns.distplot(sqrt_area, hist=True, kde=True, color='blue')
    plt.axvline(x = delta_star, linestyle='--', label=r'$\delta_{\star}$')
    plt.axvline(x = rho_s, label=r'$\rho_s$')
    plt.axvline(x = HWHM, linestyle='-.', color='black', label=r'HWHM_LM')
    plt.xlabel('cm')
    plt.legend()
    plt.show()

    
"""Probability distributions of blob amplitudes."""
def amplitude_dist(time_start,time_end):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    amplitudes = []
    for time in range(time_start,time_end):
        for row in c.execute("SELECT ROWID,* FROM t_"+str(time)):
            if ~pd.isnull(row[1]) and row[4] is not None and float(row[1])<1.5:
                amplitudes.append(row[1])
            else:
                pass
    #file = open("/global/cscratch1/sd/giannos/cmod_amplitudes(fig.9.5).txt",'a')
    #file.write("List of blob amplitudes"+"\n")
    #for i in range(len(amplitudes)):
    #    file.write(str(amplitudes[i])+"\n")
    #file.close()    
    plt.title("Amplitude distribution")
    sns.distplot(amplitudes,color='blue',hist=True,kde=True)
    plt.xlabel(r'$\frac{\delta n}{n}$')
    plt.show()
    
"""R,Z,radial and poloidal velocity distributions."""    
def total_velocity_distributions(H):
    off = 0
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    vR_l = []
    vZ_l = []
    vPSI_l = []
    vPOL_l = []
    row_counter=0
    for row in H:
        row_counter=row_counter+1
        if np.count_nonzero(row)>=2:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+off)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1+off)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
                if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                    pass
                else:
                    V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                    V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                    
                    if V_psi>10000:
                        print('time1:',fnz+ind+off-1,'rowid:',int(row[fnz+ind-1]),blobA[0],blobA[1],blobA[2],blobA[3],blobA[4],blobA[5],blobA[6],blobA[7])
                        print('time2:',fnz+ind+off+1-1,'rowid:',int(row[fnz+ind+1-1]),blobB[0],blobB[1],blobB[2],blobB[3],blobB[4],blobB[5],blobB[6],blobB[7])
                        print('row',row_counter-1)
                
                    vR_l.append(vR)
                    vZ_l.append(vZ)
                    vPSI_l.append(V_psi)
                    vPOL_l.append(V_pol)
    
    vrmean = np.mean(vR_l)
    vrmed = np.median(vR_l)
    vrkurt = kurtosis(vR_l)
    vrskew = skew(vR_l)
    vrstd = tstd(vR_l)
    
    vzmean = np.mean(vZ_l)
    vzmed = np.median(vZ_l)
    vzkurt = kurtosis(vZ_l)
    vzskew = skew(vZ_l)
    vzstd = tstd(vZ_l)
    
    vpsimean = np.mean(vPSI_l)
    vpsimed = np.median(vPSI_l)
    vpsikurt = kurtosis(vPSI_l)
    vpsiskew = skew(vPSI_l)
    vpsistd = tstd(vPSI_l)
    
    vpolmean = np.mean(vPOL_l)
    vpolmed = np.median(vPOL_l)
    vpolkurt = kurtosis(vPOL_l)
    vpolskew = skew(vPOL_l)
    vpolstd = tstd(vPOL_l)
    
    print('VRm',vrmean)
    print('VRmed',vrmed)
    print('VRkurt',vrkurt)
    print('VRsk',vrskew)
    print('VRstd',vrstd)
    
    print('VZm',vzmean)
    print('VZmed',vzmed)
    print('VZkurt',vzkurt)
    print('VZsk',vzskew)
    print('VZstd',vzstd)
    
    print('VPSIm',vpsimean)
    print('VPSImed',vpsimed)
    print('VPSIkurt',vpsikurt)
    print('VPSIsk',vpsiskew)
    print('VPSIstd',vpsistd)
    
    print('VPOLm',vpolmean)
    print('VPOLmed',vpolmed)
    print('VPOLkurt',vpolkurt)
    print('VPOLsk',vpolskew)
    print('VPOLstd',vpolstd)
                    
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(20,20))
    ax1.set_title(r"$V_R$ distribution")
    ax1.hist(vR_l,bins='auto',color='blue')
    ax1.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax2.set_title(r"$V_Z$ distribution")
    ax2.hist(vZ_l,bins='auto',color='blue')
    ax2.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax3.set_title("Radial velocity distribution")
    ax3.hist(vPSI_l,bins='auto',color='blue')
    ax3.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax4.set_title("Poloidal velocity distribution")
    ax4.hist(vPOL_l,bins='auto',color='blue')
    ax4.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    plt.show()

"""Area distribution vs lifetime."""    
def area_dist_persistance(H):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    areas = []
    areas2 = []
    areas3 = []
    areas4 = []
    areas5 = []
    areas6 = []
    for time in range(0,688):
        for row in c.execute("SELECT ROWID,* FROM t_"+str(time)):
            if ~pd.isnull(row[4]) and row[4] is not None and row[4]>5e-5:
                areas.append(row[4])
    for row in H:
        if np.count_nonzero(row)==2:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(2):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[4]) and blob[4] is not None and blob[4]>5e-5:
                    areas2.append(blob[4])
        if np.count_nonzero(row)==3:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(3):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[4]) and blob[4] is not None and blob[4]>5e-5:
                    areas3.append(blob[4])
        if np.count_nonzero(row)==4:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(4):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[4]) and blob[4] is not None and blob[4]>5e-5:
                    areas4.append(blob[4])
        if np.count_nonzero(row)==5:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(5):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[4]) and blob[4] is not None and blob[4]>5e-5:
                    areas5.append(blob[4])
        if np.count_nonzero(row)==6:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(6):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[4]) and blob[4] is not None and blob[4]>5e-5:
                    areas6.append(blob[4])
    plt.figure(figsize = (10,10))
    plt.title(r"$Area^{\frac{1}{2}}$ distribution vs. persistance")    
    sns.distplot(np.sqrt(areas),hist=True,kde=True,label='all')
    sns.distplot(np.sqrt(areas2),hist=True,kde=True,label='t=2')
    sns.distplot(np.sqrt(areas3),hist=True,kde=True,label='t=3')
    sns.distplot(np.sqrt(areas4),hist=True,kde=True,label='t=4')
    sns.distplot(np.sqrt(areas5),hist=True,kde=True,label='t=5')
    sns.distplot(np.sqrt(areas6),hist=True,kde=True,label='t=6')
    plt.xlabel('cm')
    plt.legend()
    plt.show()

"""Amplitude distribution vs. lifetime."""    
def amplitude_dist_persistance(H):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    ampl = []
    ampl2 = []
    ampl3 = []
    ampl4 = []
    ampl5 = []
    ampl6 = []
    for time in range(0,688):
        for row in c.execute("SELECT ROWID,* FROM t_"+str(time)):
            if ~pd.isnull(row[1]) and row[4] is not None and row[4]>5e-5:
                ampl.append(row[1])
    for row in H:
        if np.count_nonzero(row)==2:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(2):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[1]) and blob[1] is not None and blob[4]>5e-5:
                    ampl2.append(blob[1])
        if np.count_nonzero(row)==3:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(3):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[1]) and blob[1] is not None and blob[4]>5e-5:
                    ampl3.append(blob[1])
        if np.count_nonzero(row)==4:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(4):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[1]) and blob[1] is not None and blob[4]>5e-5:
                    ampl4.append(blob[1])
        if np.count_nonzero(row)==5:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(5):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[1]) and blob[1] is not None and blob[4]>5e-5:
                    ampl5.append(blob[1])
        if np.count_nonzero(row)==6:
            fnz = next((i for i,x in enumerate(row) if x),None)
            for i in range(6):
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+i)+" WHERE ROWID="+str(int(row[fnz+i])))
                blob = c.fetchall()[0]
                if ~pd.isnull(blob[1]) and blob[1] is not None and blob[4]>5e-5:
                    ampl6.append(blob[1])
    plt.figure(figsize = (10,10))
    plt.title("Amplitude distribution vs. persistance")    
    sns.distplot(ampl,hist=True,kde=True,label='all')
    sns.distplot(ampl2,hist=True,kde=True,label='t=2')
    sns.distplot(ampl3,hist=True,kde=True,label='t=3')
    sns.distplot(ampl4,hist=True,kde=True,label='t=4')
    sns.distplot(ampl5,hist=True,kde=True,label='t=5')
    sns.distplot(ampl6,hist=True,kde=True,label='t=6')
    plt.xlabel(r'$\frac{\delta n}{n}$')
    plt.legend()
    plt.show()

"""Scatter plots of area/amplitude vs. velocity."""    
def amplitude_velocity(H):
    off=0
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    vPSI_l = []
    vPOL_l = []
    ampl = []
    area = []
    for row in H:
        if np.count_nonzero(row)>=3:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+off)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1+off)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
    
                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                if blobA[1]>=0.2:
                    if type(V_psi) == type(V_pol) == float:
                        vPSI_l.append(V_psi)
                        vPOL_l.append(V_pol)
                    ampl.append(blobA[1])
                    area.append(blobA[4])
                else:
                    pass
                
    #fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2, figsize=(20,20))
    fig, ax1 = plt.subplots(figsize=(10,10))
    #ax1.scatter(ampl,vPSI_l)
    #ax1.set_ylabel(r'$V_{\psi}$',fontsize=20)
    #ax1.set_xlabel('Amplitude',fontsize=20)
    #ax2.scatter(ampl,vPOL_l)
    #ax2.set_ylabel(r'$V_{\theta}$',fontsize=20)
    #ax2.set_xlabel(r'Amplitude',fontsize=20)
    #ax3.scatter(np.sqrt(area),vPSI_l)
    #ax3.set_ylabel(r'$V_{\psi}$',fontsize=20)
    #ax3.set_xlabel(r'$Area^{\frac{1}{2}}\,(cm)$',fontsize=20)
    #ax4.scatter(np.sqrt(area),vPOL_l)
    #ax4.set_ylabel(r'$V_{\theta}$',fontsize=20)
    #ax4.set_xlabel(r'$Area^{\frac{1}{2}}\,(cm)$',fontsize=20)
    ax1.scatter(np.sqrt(area),ampl)
    ax1.set_xlabel(r'$Area^{\frac{1}{2}}\,(cm)$',fontsize=20)
    ax1.set_ylabel('Amplitude',fontsize=20)
    #ax2.scatter(np.sqrt(area),ampl)
    plt.show()

"""Velocity distributions vs. lifetime."""    
def velocity_persistance_distributions(H):
    warnings.filterwarnings("ignore")
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    vPSI_2 = []
    vPSI_3 = []
    vPSI_4 = []
    vPSI_5 = []
    vPSI_6 = []
    for row in H:
        if np.count_nonzero(row)==2:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
    
                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                vPSI_2.append(V_psi)
        if np.count_nonzero(row)==3:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
    
                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                vPSI_3.append(V_psi)
        if np.count_nonzero(row)==4:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
    
                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                vPSI_4.append(V_psi)
        if np.count_nonzero(row)==5:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
    
                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                vPSI_5.append(V_psi)
        if np.count_nonzero(row)==6:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
    
                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                vPSI_6.append(V_psi)
                
    plt.figure(figsize = (10,5))
    plt.title("Radial velocity distribution vs. persistance")    
    sns.distplot(vPSI_2,hist=True,kde=True,label='t=2')
    sns.distplot(vPSI_3,hist=True,kde=True,label='t=3')
    sns.distplot(vPSI_4,hist=True,kde=True,label='t=4')
    sns.distplot(vPSI_5,hist=True,kde=True,label='t=5')
    sns.distplot(vPSI_6,hist=True,kde=True,label='t=6')
    plt.xlabel(r'$\frac{m}{s}$')
    plt.legend()
    plt.show()
    
"""DIII-D: dR/dpsi = 3-4 mm. Smallest blob size is almost half that. You have blobs up to 3-4 times this size. Play with bin sizes for psi between dpsi = 0.01-0.04"""

def position_distribution(time_start,time_end):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    positions = []
    for time in range(time_start,time_end):
        for row in c.execute("SELECT ROWID,* FROM t_"+str(time)):
            if ~pd.isnull(row[1]) and ~np.isnan(row[2]) and ~np.isnan(row[3]) and row[4] is not None:
                rad_pos = psi[int(row[2]),int(row[3])]
                if np.ma.is_masked(rad_pos) == False:
                    positions.append(rad_pos)
                else:
                    pass
            else:
                pass
    plt.title("Positions distribution")
    sns.distplot(positions,color='blue',hist=True,kde=True)
    plt.xlabel(r'$\Psi_N$')
    plt.show()
    
def birth_distribution(H):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    off = 0
    
    positions = []
    for row in H:
        fnz = next((i for i,x in enumerate(row) if x),None)
        c.execute("SELECT ROWID,* FROM t_"+str(fnz+off)+" WHERE ROWID="+str(int(row[fnz])))
        blob = c.fetchall()[0]
        if ~pd.isnull(blob[1]) and ~np.isnan(blob[2]) and ~np.isnan(blob[3]) and blob[4] is not None:
                rad_pos = psi[int(blob[2]),int(blob[3])]
                if np.ma.is_masked(rad_pos) == False:
                    positions.append(rad_pos)
                else:
                    pass
        else:
            pass
    plt.title("Birth position distribution")
    sns.distplot(positions,color='blue',hist=True,kde=True)
    plt.xlabel(r'$\Psi_N$')
    plt.show()        
            
def vel_vs_pos(H):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    off = 0
    
    #break the results in 5 blocks: psi<0.995, 0.995-1.000, 1.000-1.005, 1.005-1.010, psi>1.010
    vR_1 = []
    vZ_1 = []
    vPSI_1 = []
    vPOL_1 = []
    vR_2 = []
    vZ_2 = []
    vPSI_2 = []
    vPOL_2 = []
    vR_3 = []
    vZ_3 = []
    vPSI_3 = []
    vPOL_3 = []
    vR_4 = []
    vZ_4 = []
    vPSI_4 = []
    vPOL_4 = []
    vR_5 = []
    vZ_5 = []
    vPSI_5 = []
    vPOL_5 = []
    for row in H:
        if np.count_nonzero(row)>=2:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+off)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1+off)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos<=0.995:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                
                                vR_1.append(vR)
                                vZ_1.append(vZ)
                                vPSI_1.append(V_psi)
                                vPOL_1.append(V_pol)
                        if rad_pos>0.995 and rad_pos<=1.000:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                
                                vR_2.append(vR)
                                vZ_2.append(vZ)
                                vPSI_2.append(V_psi)
                                vPOL_2.append(V_pol)
                        if rad_pos>1.000 and rad_pos<=1.005:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                
                                vR_3.append(vR)
                                vZ_3.append(vZ)
                                vPSI_3.append(V_psi)
                                vPOL_3.append(V_pol)
                        if rad_pos>1.005 and rad_pos<=1.010:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                                
                                vR_4.append(vR)
                                vZ_4.append(vZ)
                                vPSI_4.append(V_psi)
                                vPOL_4.append(V_pol)
                                
                        if rad_pos>1.010:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                
                                vR_5.append(vR)
                                vZ_5.append(vZ)
                                vPSI_5.append(V_psi)
                                vPOL_5.append(V_pol)

                        
                    else:
                        pass
                else:
                    pass
    
    pol_m1 = np.mean(vPOL_1)
    pol_m2 = np.mean(vPOL_2)
    pol_m3 = np.mean(vPOL_3)
    pol_m4 = np.mean(vPOL_4)
    pol_m5 = np.mean(vPOL_5)
    rad_m1 = np.mean(vPSI_1)
    rad_m2 = np.mean(vPSI_2)
    rad_m3 = np.mean(vPSI_3)
    rad_m4 = np.mean(vPSI_4)
    rad_m5 = np.mean(vPSI_5)
    pol_med1 = np.median(vPOL_1)
    pol_med2 = np.median(vPOL_2)
    pol_med3 = np.median(vPOL_3)
    pol_med4 = np.median(vPOL_4)
    pol_med5 = np.median(vPOL_5)
    rad_med1 = np.median(vPSI_1)
    rad_med2 = np.median(vPSI_2)
    rad_med3 = np.median(vPSI_3)
    rad_med4 = np.median(vPSI_4)
    rad_med5 = np.median(vPSI_5)
    pol_sk1 = skew(vPOL_1)
    pol_sk2 = skew(vPOL_2)
    pol_sk3 = skew(vPOL_3)
    pol_sk4 = skew(vPOL_4)
    pol_sk5 = skew(vPOL_5)
    rad_sk1 = skew(vPSI_1)
    rad_sk2 = skew(vPSI_2)
    rad_sk3 = skew(vPSI_3)
    rad_sk4 = skew(vPSI_4)
    rad_sk5 = skew(vPSI_5)
    pol_kt1 = kurtosis(vPOL_1)
    pol_kt2 = kurtosis(vPOL_2)
    pol_kt3 = kurtosis(vPOL_3)
    pol_kt4 = kurtosis(vPOL_4)
    pol_kt5 = kurtosis(vPOL_5)
    rad_kt1 = kurtosis(vPSI_1)
    rad_kt2 = kurtosis(vPSI_2)
    rad_kt3 = kurtosis(vPSI_3)
    rad_kt4 = kurtosis(vPSI_4)
    rad_kt5 = kurtosis(vPSI_5)

                    
    fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8),(ax9,ax10)) = plt.subplots(5,2, figsize=(20,40))
    ax1.set_title(r"$V_{\Psi}$ distribution, $\Psi_N<0.995$")
    ax1.hist(vPSI_1,bins='auto',color='blue')
    ax1.axvline(x=rad_m1,color='r',label='mean')
    ax1.axvline(x=rad_med1,color='k',label='median')
    ax1.text(0.8,1,r'skewness=%.5f'%(rad_sk1),transform=ax1.transAxes)
    ax1.text(0.8,0.96,r'kurtosis=%.5f'%(rad_kt1),transform=ax1.transAxes)
    ax1.text(0.8,0.92,r'mean=%.0f'%(rad_m1),transform=ax1.transAxes)
    ax1.text(0.8,0.88,r'median=%.0f'%(rad_med1),transform=ax1.transAxes)
    ax1.legend(loc=7)
    ax1.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax2.set_title(r"$V_{\theta}$ distribution, $\Psi_N<0.995$")
    ax2.hist(vPOL_1,bins='auto',color='blue')
    ax2.axvline(x=pol_m1,color='r',label='mean')
    ax2.axvline(x=pol_med1,color='k',label='median')
    ax2.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax2.text(0.8,1,r'skewness=%.5f'%(pol_sk1),transform=ax2.transAxes)
    ax2.text(0.8,0.96,r'kurtosis=%.5f'%(pol_kt1),transform=ax2.transAxes)
    ax2.text(0.8,0.92,r'mean=%.0f'%(pol_m1),transform=ax2.transAxes)
    ax2.text(0.8,0.88,r'median=%.0f'%(pol_med1),transform=ax2.transAxes)
    ax2.legend(loc=6)
    ax3.set_title(r"$V_{\Psi}$ distribution, $0.995<\Psi_N<1.000$")
    ax3.hist(vPSI_2,bins='auto',color='blue')
    ax3.axvline(x=rad_m2,color='r',label='mean')
    ax3.axvline(x=rad_med2,color='k',label='median')
    ax3.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax3.text(0.8,1,r'skewness=%.5f'%(rad_sk2),transform=ax3.transAxes)
    ax3.text(0.8,0.96,r'kurtosis=%.5f'%(rad_kt2),transform=ax3.transAxes)
    ax3.text(0.8,0.92,r'mean=%.0f'%(rad_m2),transform=ax3.transAxes)
    ax3.text(0.8,0.88,r'median=%.0f'%(rad_med2),transform=ax3.transAxes)
    ax3.legend(loc=7)
    ax4.set_title(r"$V_{\theta}$ distribution, $0.995<\Psi_N<1.000$")
    ax4.hist(vPOL_2,bins='auto',color='blue')
    ax4.axvline(x=pol_m2,color='r',label='mean')
    ax4.axvline(x=pol_med2,color='k',label='median')
    ax4.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax4.text(0.8,1,r'skewness=%.5f'%(pol_sk2),transform=ax4.transAxes)
    ax4.text(0.8,0.96,r'kurtosis=%.5f'%(pol_kt2),transform=ax4.transAxes)
    ax4.text(0.8,0.92,r'mean=%.0f'%(pol_m2),transform=ax4.transAxes)
    ax4.text(0.8,0.88,r'median=%.0f'%(pol_med2),transform=ax4.transAxes)
    ax4.legend(loc=7)
    ax5.set_title(r"$V_{\Psi}$ distribution, $1.000<\Psi_N<1.005$")
    ax5.hist(vPSI_3,bins='auto',color='blue')
    ax5.axvline(x=rad_m3,color='r',label='mean')
    ax5.axvline(x=rad_med3,color='k',label='median')
    ax5.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax5.text(0.8,1,r'skewness=%.5f'%(rad_sk3),transform=ax5.transAxes)
    ax5.text(0.8,0.96,r'kurtosis=%.5f'%(rad_kt3),transform=ax5.transAxes)
    ax5.text(0.8,0.92,r'mean=%.0f'%(rad_m3),transform=ax5.transAxes)
    ax5.text(0.8,0.88,r'median=%.0f'%(rad_med3),transform=ax5.transAxes)
    ax5.legend(loc=7)
    ax6.set_title(r"$V_{\theta}$ distribution, $1.000<\Psi_N<1.005$")
    ax6.hist(vPOL_3,bins='auto',color='blue')
    ax6.axvline(x=pol_m3,color='r',label='mean')
    ax6.axvline(x=pol_med3,color='k',label='median')
    ax6.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax6.text(0.8,1,r'skewness=%.5f'%(pol_sk3),transform=ax6.transAxes)
    ax6.text(0.8,0.96,r'kurtosis=%.5f'%(pol_kt3),transform=ax6.transAxes)
    ax6.text(0.8,0.92,r'mean=%.0f'%(pol_m3),transform=ax6.transAxes)
    ax6.text(0.8,0.88,r'median=%.0f'%(pol_med3),transform=ax6.transAxes)
    ax6.legend(loc=7)
    ax7.set_title(r"$V_{\Psi}$ distribution, $1.005<\Psi_N<1.010$")
    ax7.hist(vPSI_4,bins='auto',color='blue')
    ax7.axvline(x=rad_m4,color='r',label='mean')
    ax7.axvline(x=rad_med4,color='k',label='median')
    ax7.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax7.text(0.8,1,r'skewness=%.5f'%(rad_sk4),transform=ax7.transAxes)
    ax7.text(0.8,0.96,r'kurtosis=%.5f'%(rad_kt4),transform=ax7.transAxes)
    ax7.text(0.8,0.92,r'mean=%.0f'%(rad_m4),transform=ax7.transAxes)
    ax7.text(0.8,0.88,r'median=%.0f'%(rad_med4),transform=ax7.transAxes)
    ax7.legend(loc=7)
    ax8.set_title(r"$V_{\theta}$ distribution, $1.005<\Psi_N<1.010$")
    ax8.hist(vPOL_4,bins='auto',color='blue')
    ax8.axvline(x=pol_m4,color='r',label='mean')
    ax8.axvline(x=pol_med4,color='k',label='median')
    ax8.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax8.text(0.8,1,r'skewness=%.5f'%(pol_sk4),transform=ax8.transAxes)
    ax8.text(0.8,0.96,r'kurtosis=%.5f'%(pol_kt4),transform=ax8.transAxes)
    ax8.text(0.8,0.92,r'mean=%.0f'%(pol_m4),transform=ax8.transAxes)
    ax8.text(0.8,0.88,r'median=%.0f'%(pol_med4),transform=ax8.transAxes)
    ax8.legend(loc=7)
    ax9.set_title(r"$V_{\Psi}$ distribution, $\Psi_N>1.010$")
    ax9.hist(vPSI_5,bins='auto',color='blue')
    ax9.axvline(x=rad_m5,color='r',label='mean')
    ax9.axvline(x=rad_med5,color='k',label='median')
    ax9.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax9.text(0.8,1,r'skewness=%.5f'%(rad_sk5),transform=ax9.transAxes)
    ax9.text(0.8,0.96,r'kurtosis=%.5f'%(rad_kt5),transform=ax9.transAxes)
    ax9.text(0.8,0.92,r'mean=%.0f'%(rad_m5),transform=ax9.transAxes)
    ax9.text(0.8,0.88,r'median=%.0f'%(rad_med5),transform=ax9.transAxes)
    ax9.legend(loc=7)
    ax10.set_title(r"$V_{\theta}$ distribution, $\Psi_N>1.010$")
    ax10.hist(vPOL_5,bins='auto',color='blue')
    ax10.axvline(x=pol_m5,color='r',label='mean')
    ax10.axvline(x=pol_med5,color='k',label='median')
    ax10.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax10.text(0.8,1,r'skewness=%.5f'%(pol_sk5),transform=ax10.transAxes)
    ax10.text(0.8,0.96,r'kurtosis=%.5f'%(pol_kt5),transform=ax10.transAxes)
    ax10.text(0.8,0.92,r'mean=%.0f'%(pol_m5),transform=ax10.transAxes)
    ax10.text(0.8,0.88,r'median=%.0f'%(pol_med5),transform=ax10.transAxes)
    ax10.legend(loc=7)
    
    plt.show()
    
    
    
def vel_pos(H):
    off=0
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    #break the results in 17 blocks from 0.88-1.05
    
    v1=[];v2=[];v3=[];v4=[];v5=[];v6=[];v7=[];v8=[];v9=[];v10=[];v11=[];v12=[];v13=[];v14=[];v15=[];v16=[];v17=[]
    vPsi=[]
    count=[]
    
    for row in H:
        if np.count_nonzero(row)>=2:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+off)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1+off)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos<=0.89:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v1.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.89 and rad_pos<=0.9:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v2.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.9 and rad_pos<=0.91:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v3.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.91 and rad_pos<=0.92:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v4.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.92 and rad_pos<=0.93:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v5.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.93 and rad_pos<=0.94:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v6.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.94 and rad_pos<=0.95:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v7.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.95 and rad_pos<=0.96:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v8.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.96 and rad_pos<=0.97:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v9.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.97 and rad_pos<=0.98:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v10.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.98 and rad_pos<=0.99:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v11.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>0.99 and rad_pos<=1.0:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v12.append(V_psi)
            
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>1.0 and rad_pos<=1.01:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v13.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>1.01 and rad_pos<=1.02:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v14.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>1.02 and rad_pos<=1.03:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v15.append
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>1.03 and rad_pos<=1.04:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v16.append(V_psi)
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if ma.is_masked(rad_pos) == False:
                        if rad_pos>1.04:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                v17.append(V_psi)
    
    vPsi.append(np.mean(v1))
    vPsi.append(np.mean(v2))
    vPsi.append(np.mean(v3))
    vPsi.append(np.mean(v4))
    vPsi.append(np.mean(v5))
    vPsi.append(np.mean(v6))
    vPsi.append(np.mean(v7))
    vPsi.append(np.mean(v8))
    vPsi.append(np.mean(v9))
    vPsi.append(np.mean(v10))
    vPsi.append(np.mean(v11))
    vPsi.append(np.mean(v12))
    vPsi.append(np.mean(v13))
    vPsi.append(np.mean(v14))
    vPsi.append(np.mean(v15))
    vPsi.append(np.mean(v16))
    vPsi.append(np.mean(v17))
    
    count.append(len(v1))
    count.append(len(v2))
    count.append(len(v3))
    count.append(len(v4))
    count.append(len(v5))
    count.append(len(v6))
    count.append(len(v7))
    count.append(len(v8))
    count.append(len(v9))
    count.append(len(v10))
    count.append(len(v11))
    count.append(len(v12))
    count.append(len(v13))
    count.append(len(v14))
    count.append(len(v15))
    count.append(len(v16))
    count.append(len(v17))
                
    x = np.linspace(0.88,1.05,17)
    
    print(count)
    print(vPsi)
    
    plt.plot(x,vPsi)
    plt.xlabel(r'$\Psi_N$')
    plt.ylabel(r'$V_{\psi}$')
    plt.show()
                
    

def flux_vs_pos(arr_mid,H):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    off = 0
    
    arr_tor = arr_mid.mean(1)
    arr_avg = arr_tor.mean(0)
    arr_avg = np.swapaxes(arr_avg,0,1)
    
    #break the results in 5 blocks: psi<0.995, 0.995-1.000, 1.000-1.005, 1.005-1.010, psi>1.010
    flux_1 = []
    flux_2 = []
    flux_3 = []
    flux_4 = []
    flux_5 = []
    
    for row in H:
        if np.count_nonzero(row)>=2:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+off)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1+off)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                denA = blobA[1]*arr_avg[int(blobA[2]),int(blobA[3])]
                denB = blobB[1]*arr_avg[int(blobB[2]),int(blobB[3])]
                den = (denA+denB)/2.0
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    rad_pos = psi[int(blobA[2]),int(blobA[3])]
                    if np.ma.is_masked(rad_pos) == False:
                        if rad_pos<=0.995:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                flux = den*V_psi
                                
                                flux_1.append(flux)
                                
                        if rad_pos>0.995 and rad_pos<=1.000:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                flux = den*V_psi
                                
                                flux_2.append(flux)
                                
                        if rad_pos>1.000 and rad_pos<=1.005:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                flux = den*V_psi
                                
                                flux_3.append(flux)
                                
                        if rad_pos>1.005 and rad_pos<=1.010:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                flux = den*V_psi
                                
                                flux_4.append(flux)
                        
                        if rad_pos>1.010:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                flux = den*V_psi
                                
                                flux_5.append(flux)
                        
                    else:
                        pass
                else:
                    pass
    
    flux_m1 = np.mean(flux_1)
    flux_m2 = np.mean(flux_2)
    flux_m3 = np.mean(flux_3)
    flux_m4 = np.mean(flux_4)
    flux_m5 = np.mean(flux_5)
    flux_med1 = np.median(flux_1)
    flux_med2 = np.median(flux_2)
    flux_med3 = np.median(flux_3)
    flux_med4 = np.median(flux_4)
    flux_med5 = np.median(flux_5)
    flux_sk1 = skew(flux_1)
    flux_sk2 = skew(flux_2)
    flux_sk3 = skew(flux_3)
    flux_sk4 = skew(flux_4)
    flux_sk5 = skew(flux_5)
    flux_kt1 = kurtosis(flux_1)
    flux_kt2 = kurtosis(flux_2)
    flux_kt3 = kurtosis(flux_3)
    flux_kt4 = kurtosis(flux_4)
    flux_kt5 = kurtosis(flux_5)
                    
    fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2, figsize=(20,40))
    ax1.set_title(r"Flux distribution, $\Psi_N<0.995$")
    ax1.hist(flux_1,bins='auto',color='blue')
    ax1.axvline(x=flux_m1,color='r',label='mean')
    ax1.axvline(x=flux_med1,color='k',label='median')
    ax1.text(0.8,1,r'skewness=%.5f'%(flux_sk1),transform=ax1.transAxes)
    ax1.text(0.8,0.96,r'kurtosis=%.5f'%(flux_kt1),transform=ax1.transAxes)
    ax1.text(0.8,0.92,r'mean=%.2E'%(flux_m1),transform=ax1.transAxes)
    ax1.text(0.8,0.88,r'median=%.2E'%(flux_med1),transform=ax1.transAxes)
    ax1.legend(loc=7)
    ax1.set_xlabel(r"$m^{-2}s^{-1}$",fontsize=20)
    ax2.set_title(r"Flux distribution, $0.995<\Psi_N<1.000$")
    ax2.hist(flux_2,bins='auto',color='blue')
    ax2.axvline(x=flux_m2,color='r',label='mean')
    ax2.axvline(x=flux_med2,color='k',label='median')
    ax2.set_xlabel(r"$m^{-2}s^{-1}$",fontsize=20)
    ax2.text(0.8,1,r'skewness=%.5f'%(flux_sk2),transform=ax2.transAxes)
    ax2.text(0.8,0.96,r'kurtosis=%.5f'%(flux_kt2),transform=ax2.transAxes)
    ax2.text(0.8,0.92,r'mean=%.2E'%(flux_m2),transform=ax2.transAxes)
    ax2.text(0.8,0.88,r'median=%.2E'%(flux_med2),transform=ax2.transAxes)
    ax2.legend(loc=7)
    ax3.set_title(r"Flux distribution, $1.000<\Psi_N<1.005$")
    ax3.hist(flux_3,bins='auto',color='blue')
    ax3.axvline(x=flux_m3,color='r',label='mean')
    ax3.axvline(x=flux_med3,color='k',label='median')
    ax3.set_xlabel(r"$m^{-2}s^{-1}$",fontsize=20)
    ax3.text(0.8,1,r'skewness=%.5f'%(flux_sk3),transform=ax3.transAxes)
    ax3.text(0.8,0.96,r'kurtosis=%.5f'%(flux_kt3),transform=ax3.transAxes)
    ax3.text(0.8,0.92,r'mean=%.2E'%(flux_m3),transform=ax3.transAxes)
    ax3.text(0.8,0.88,r'median=%.2E'%(flux_med3),transform=ax3.transAxes)
    ax3.legend(loc=7)
    ax4.set_title(r"Flux distribution, $1.005<\Psi_N<1.010$")
    ax4.hist(flux_4,bins='auto',color='blue')
    ax4.axvline(x=flux_m4,color='r',label='mean')
    ax4.axvline(x=flux_med4,color='k',label='median')
    ax4.set_xlabel(r"m^{-2}s^{-1}",fontsize=20)
    ax4.text(0.8,1,r'skewness=%.5f'%(flux_sk4),transform=ax4.transAxes)
    ax4.text(0.8,0.96,r'kurtosis=%.5f'%(flux_kt4),transform=ax4.transAxes)
    ax4.text(0.8,0.92,r'mean=%.2E'%(flux_m4),transform=ax4.transAxes)
    ax4.text(0.8,0.88,r'median=%.2E'%(flux_med4),transform=ax4.transAxes)
    ax4.legend(loc=7)
    ax5.set_title(r"Flux distribution, $\Psi_N>1.010$")
    ax5.hist(flux_5,bins='auto',color='blue')
    ax5.axvline(x=flux_m5,color='r',label='mean')
    ax5.axvline(x=flux_med5,color='k',label='median')
    ax5.set_xlabel(r"$m^{-2}s^{-1}$",fontsize=20)
    ax5.text(0.8,1,r'skewness=%.5f'%(flux_sk5),transform=ax5.transAxes)
    ax5.text(0.8,0.96,r'kurtosis=%.5f'%(flux_kt5),transform=ax5.transAxes)
    ax5.text(0.8,0.92,r'mean=%.2E'%(flux_m5),transform=ax5.transAxes)
    ax5.text(0.8,0.88,r'median=%.2E'%(flux_med5),transform=ax5.transAxes)
    ax5.legend(loc=7)
    ax6.set_title(r"Flux distribution, $\Psi_N>1.010$")
    ax6.hist(flux_5,bins='auto',color='blue')
    ax6.axvline(x=flux_m5,color='r',label='mean')
    ax6.axvline(x=flux_med5,color='k',label='median')
    ax6.set_xlabel(r"$m^{-2}s^{-1}$",fontsize=20)
    ax6.text(0.8,1,r'skewness=%.5f'%(flux_sk5),transform=ax5.transAxes)
    ax6.text(0.8,0.96,r'kurtosis=%.5f'%(flux_kt5),transform=ax5.transAxes)
    ax6.text(0.8,0.92,r'mean=%.2E'%(flux_m5),transform=ax5.transAxes)
    ax6.text(0.8,0.88,r'median=%.2E'%(flux_med5),transform=ax5.transAxes)
    ax6.legend(loc=7)
    
    plt.show()
    
def norm_atan(y,x):
    '''Returns the angle formed between a vector with end coordinates x,y and the ray starting at the origin and passing through the
       point(1,0). The angle is between 0 and 360 degrees and the function can also take arrays.
    '''
    angle = np.degrees(np.arctan2(y,x))
    try:
        for index, elem in enumerate(angle):
            if elem<0:
                angle[index] = 360+elem
            else:
                pass
    except TypeError:
        if angle<0:
            angle = 360+angle
        else:
            pass
    
    return angle


def poloidal_range():
    print('top:',norm_atan(Z_points[1],R_points[9]))
    print('bottom:',norm_atan(Z_points[98],R_points[26]))
    
    
    
def vel_vs_pol_pos(H):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    off = 0
    
    #break the results in 2 blocks: theta<360 (below midplane), theta>0 (above midplane)
    vR_1 = []
    vZ_1 = []
    vPSI_1 = []
    vPOL_1 = []
    vR_2 = []
    vZ_2 = []
    vPSI_2 = []
    vPOL_2 = []
    
    for row in H:
        if np.count_nonzero(row)>=2:
            fnz = next((i for i,x in enumerate(row) if x),None)
            ind = 0
            
            while fnz+ind+1<len(row) and row[fnz+ind+1]!=0:
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+off)+" WHERE ROWID="+str(int(row[fnz+ind])))
                blobA = c.fetchall()[0]
                c.execute("SELECT ROWID,* FROM t_"+str(fnz+ind+1+off)+" WHERE ROWID="+str(int(row[fnz+ind+1])))
                blobB = c.fetchall()[0]
                ind = ind+1
                vR = ((blobB[6]-blobA[6])*unitR)/dt
                vZ = ((blobB[7]-blobA[7])*unitZ)/dt
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Br = (BR[int(blobA[3])][int(blobA[2])]+BR[int(blobB[3])][int(blobB[2])])/2
                Bz = (BZ[int(blobA[3])][int(blobA[2])]+BZ[int(blobB[3])][int(blobB[2])])/2
                Bpol = (Bp[int(blobA[3])][int(blobA[2])]+Bp[int(blobB[3])][int(blobB[2])])/2
                
                if ~pd.isnull(blobA[1]) and ~np.isnan(blobA[2]) and ~np.isnan(blobA[3]) and blobA[4] is not None:
                    pol_ang = norm_atan(Z_points[int(blobA[3])],R_points[int(blobA[2])])         
                    if np.ma.is_masked(pol_ang) == False:
                        if pol_ang>350 and pol_ang<=359.999:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                
                                vR_1.append(vR)
                                vZ_1.append(vZ)
                                vPSI_1.append(V_psi)
                                vPOL_1.append(V_pol)
                        if pol_ang>=0.0 and pol_ang<=50.00:
                            if np.isnan(Bpol) or np.isnan(Br) or np.isnan(Bz):
                                pass
                            else:
                                V_psi = vR*(Bz/Bpol) - vZ*(Br/Bpol)
                                V_pol = vR*(Br/Bpol) + vZ*(Bz/Bpol)
                
                                vR_2.append(vR)
                                vZ_2.append(vZ)
                                vPSI_2.append(V_psi)
                                vPOL_2.append(V_pol)
                        
                        
                    else:
                        pass
                else:
                    pass
    
    pol_m1 = np.mean(vPOL_1)
    pol_m2 = np.mean(vPOL_2)
    rad_m1 = np.mean(vPSI_1)
    rad_m2 = np.mean(vPSI_2)
    pol_med1 = np.median(vPOL_1)
    pol_med2 = np.median(vPOL_2)
    rad_med1 = np.median(vPSI_1)
    rad_med2 = np.median(vPSI_2)
    pol_sk1 = skew(vPOL_1)
    pol_sk2 = skew(vPOL_2)
    rad_sk1 = skew(vPSI_1)
    rad_sk2 = skew(vPSI_2)
    pol_kt1 = kurtosis(vPOL_1)
    pol_kt2 = kurtosis(vPOL_2)
    rad_kt1 = kurtosis(vPSI_1)
    rad_kt2 = kurtosis(vPSI_2)    
                    
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize=(20,40))
    ax1.set_title(r"Radial Velocity, below midplane")
    ax1.hist(vPSI_1,bins='auto',color='blue')
    ax1.axvline(x=rad_m1,color='r',label='mean')
    ax1.axvline(x=rad_med1,color='k',label='median')
    ax1.text(0.8,1,r'skewness=%.5f'%(rad_sk1),transform=ax1.transAxes)
    ax1.text(0.8,0.96,r'kurtosis=%.5f'%(rad_kt1),transform=ax1.transAxes)
    ax1.text(0.8,0.92,r'mean=%.0f'%(rad_m1),transform=ax1.transAxes)
    ax1.text(0.8,0.88,r'median=%.0f'%(rad_med1),transform=ax1.transAxes)
    ax1.legend(loc=7)
    ax1.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax2.set_title(r"Poloidal Velocity, below midplane")
    ax2.hist(vPOL_1,bins='auto',color='blue')
    ax2.axvline(x=pol_m1,color='r',label='mean')
    ax2.axvline(x=pol_med1,color='k',label='median')
    ax2.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax2.text(0.8,1,r'skewness=%.5f'%(pol_sk1),transform=ax2.transAxes)
    ax2.text(0.8,0.96,r'kurtosis=%.5f'%(pol_kt1),transform=ax2.transAxes)
    ax2.text(0.8,0.92,r'mean=%.0f'%(pol_m1),transform=ax2.transAxes)
    ax2.text(0.8,0.88,r'median=%.0f'%(pol_med1),transform=ax2.transAxes)
    ax2.legend(loc=6)
    ax3.set_title(r"Radial Velocity, above midplane")
    ax3.hist(vPSI_2,bins='auto',color='blue')
    ax3.axvline(x=rad_m2,color='r',label='mean')
    ax3.axvline(x=rad_med2,color='k',label='median')
    ax3.text(0.8,1,r'skewness=%.5f'%(rad_sk2),transform=ax3.transAxes)
    ax3.text(0.8,0.96,r'kurtosis=%.5f'%(rad_kt2),transform=ax3.transAxes)
    ax3.text(0.8,0.92,r'mean=%.0f'%(rad_m2),transform=ax3.transAxes)
    ax3.text(0.8,0.88,r'median=%.0f'%(rad_med2),transform=ax3.transAxes)
    ax3.legend(loc=7)
    ax3.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax4.set_title(r"Poloidal Velocity, above midplane")
    ax4.hist(vPOL_2,bins='auto',color='blue')
    ax4.axvline(x=pol_m2,color='r',label='mean')
    ax4.axvline(x=pol_med2,color='k',label='median')
    ax4.set_xlabel(r"$\frac{m}{s}$",fontsize=20)
    ax4.text(0.8,1,r'skewness=%.5f'%(pol_sk2),transform=ax4.transAxes)
    ax4.text(0.8,0.96,r'kurtosis=%.5f'%(pol_kt2),transform=ax4.transAxes)
    ax4.text(0.8,0.92,r'mean=%.0f'%(pol_m2),transform=ax4.transAxes)
    ax4.text(0.8,0.88,r'median=%.0f'%(pol_med2),transform=ax4.transAxes)
    ax4.legend(loc=6)
    
    plt.show()
