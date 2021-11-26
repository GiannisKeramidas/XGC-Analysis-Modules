import random
import numpy as np
import matplotlib.pyplot as plt
import core

def sep_file_read():
    #sep_file = open("Separatrix nodes locations","r")
    sep_file = open("/global/homes/g/giannos/xgc_python_dir/ti262_analysis/ti262_Separatrix.txt","r")
    next(sep_file)
    R=[]
    Z=[]
    for columns in (raw.strip().split() for raw in sep_file):
        R.append(float(columns[0]))
        Z.append(float(columns[1]))
    return R,Z

R,Z = sep_file_read()

Ang=[]
for i in range(0,len(R)):
    R_val = R[i]-core.Rmaj
    Z_val = Z[i]
    angle = np.degrees(np.arctan2(Z_val,R_val))
    if angle<0:
        angle = 360 + angle
    else:
        pass
    Ang.append(angle)
L = np.array([x for x in range(0,len(R))])
#plt.plot(L,Ang,'b-')
fig, ax=plt.subplots()
plt.title('Loading patches (ti262)')
plt.xlabel('R(m)')
plt.ylabel('Z(m)')
for i in range(0,80):#you need 80 processes for ti262 and 60 for ti255
    ax.plot(R[20*i:20*(i+1)],Z[20*i:20*(i+1)],color=(np.random.rand(3)))
plt.grid(True)
plt.show()
