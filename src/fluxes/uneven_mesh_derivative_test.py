import numpy as np
import math
import pylab as py
import matplotlib.pyplot as plt





def fun2():
    x1 = np.linspace(0,100)
    y1 = np.power(x1,3)
    new_list1 = [i*3 for i in x1]
    z1 = np.power(new_list1,2)
    plt.plot(x1[:],y1[:],'go', x1[:],z1[:],'rs')
    plt.show()
def my_lin(lb, ub, steps, spacing=2.6):
    span = (ub-lb)
    dx = 1.0 / (steps-1)
    return [lb + (i*dx)**spacing*span for i in range(steps)]

def fun3(stop):
    x = my_lin(0,stop,40)
    Rp = np.roll(x,1)
    Rm = np.roll(x,-1)
    dRp = x - Rp
    dRm = Rm - x
    y = np.power(x,3)
    derivative = []
    for i in range(1,len(y)-1):
        A = dRp[i+1]/(dRm[i-1]**2 + dRp[i+1]*dRm[i-1]) - dRm[i-1]/(dRp[i+1]**2 + dRp[i+1]*dRm[i-1])
        B = -dRp[i+1]/(dRm[i-1]**2 + dRp[i+1]*dRm[i-1])
        C = dRm[i-1]/(dRp[i+1]**2 + dRp[i+1]*dRm[i-1])
        derivative.append(A*y[i] + B*y[i-1] + C*y[i+1])
    derivative.append((y[-1]-y[-2])/dRp[-1])
    derivative = [(y[1]-y[0])/dRm[0]] + derivative

    new_list = [i*3 for i in x]
    z = 3*np.power(x,2)
    x1 = np.linspace(0,stop)
    y1 = np.power(x1,3)
    y1plus = np.roll(y1,-1)
    y1minus = np.roll(y1,1)
    x1plus = np.roll(x1,-1)
    x1minus = np.roll(x1,1)
    der2 = (y1plus-y1minus)/(x1plus-x1minus)
    der2[0] = (y1[1]-y1[0])/(x1[1]-x1[0])
    der2[-1] = (y1[-1]-y1[-2])/(x1[-1]-x1[-2])
    fig, ax =plt.subplots()
    ax.plot(x[:],y[:],'go',label='$x^3$')
    ax.plot(x[:],z[:],'rs',label='$3 x^2$')
    ax.plot(x[:],derivative[:],'y*',label='uneven der.')
    ax.plot(x1[:],der2[:],'m*',label='even der.')
    legend = ax.legend()
    plt.show()
