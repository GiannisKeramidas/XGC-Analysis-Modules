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
def file_read():
    tolerance = 0.5
    crs = []
    for i in range(0,53):
        crs.append(open("The_3_fluxes_vs.angle_%d.txt" %(i),"r"))
    #for i in range(0,24):
        #next(crs[i])
    angle = []
    eq_fl = []
    tur_fl = []
    for i in range(0,53):
        for columns in ( raw.strip().split() for raw in crs[i] ):  
            angle.append(float(columns[0]))
            eq_fl.append(float(columns[1]))
            tur_fl.append(float(columns[2]))
    excl_val = []
    excl_val_2 = []
    for elem in np.where(np.isnan(eq_fl))[0][:]:#filter out nan values from fluxes
        excl_val_2.append(elem)
    for elem in sorted(excl_val_2,reverse=True):
        del angle[elem]
        del eq_fl[elem]
        del tur_fl[elem]
    excl_val_3 = []
    for elem in np.where(np.isnan(tur_fl))[0][:]:
        excl_val_3.append(elem)
    for elem in sorted(excl_val_3,reverse=True):
        del angle[elem]
        del eq_fl[elem]
        del tur_fl[elem]
    for i in range(0,len(angle)-1):#remove duplicates in angle
        for j in range(i+1,len(angle)-1):
            diff = abs(angle[i]-angle[j])
            if (diff<tolerance) or (diff>360-tolerance):
                excl_val.append(j)
        for elem in sorted(excl_val,reverse=True):
            del angle[elem]
            del eq_fl[elem]
            del tur_fl[elem]
        excl_val = []
    print(angle)
    nul_pos = angle.index(min(angle))#start all values from zero angle
    #for i in range(0,len(angle)):
        #if angle[i]<0.8:
            #nul_pos = i
    angle = np.roll(angle,-nul_pos)
    eq_fl = np.roll(eq_fl,-nul_pos)
    tur_fl = np.roll(tur_fl,-nul_pos)
    print(angle)
    plt.xlabel("angle")
    plt.ylabel("fluxes")
    plt.plot(angle[:], eq_fl[:],'b-',angle[:], tur_fl[:],'r-')
    plt.show()



    return angle, eq_fl, tur_fl
file_read()
