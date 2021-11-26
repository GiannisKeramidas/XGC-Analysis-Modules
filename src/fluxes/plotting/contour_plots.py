import numpy as np
import xgc
import math
import pylab as py
import matplotlib.pyplot as plt
import sys
import numpy.ma as ma
import core
from decimal import Decimal
import os
import angle
import glob
import matplotlib.animation as animation
import types

#fileDir = '/scratch1/scratchdirs/giannos/ti255_d3d_Ip_1.5_med'
fileDir='/global/cscratch1/sd/giannos/ti262_cmod_JTM'

trgt = "/global/cscratch1/sd/giannos/"

#DIII-D
#Rmin = 1.97#None#2.15
#Rmax = 2.26#None#2.26
#Zmin = 0.1#None#-0.28
#Zmax = -0.5#None#0.55

#CMOD midplane
Rmin = 0.88
Rmax = 0.908
Zmin = -0.1
Zmax = 0.1

L = np.array([x for x in range(0,100)])
#unitR = (Rmax-Rmin)/len(L)
#R_points = np.array([L[ip]*unitR+Rmin for ip in range(0,len(L))])
#unitZ = (Zmax-Zmin)/len(L)
#Z = np.array([y for y in range(0,100)])
#Z_points = np.array([Z[ip]*unitZ+Zmin for ip in range(0,len(L))])
 
core.getMeshAndCuts(fileDir,Rmin,Rmax,Zmin,Zmax)
def write_sep():
    #file = open("Wedge_Sep-ti255_med.txt",'a')
    file = open("Full_Sep-ti255_med.txt",'a')
    file.write("R" + "\t" + "Z" + "\n")
    Rs_point = core.loader.RZ[core.sepInds,0]
    Zs_point = core.loader.RZ[core.sepInds,1]
    for i in range(0,len(Rs_point)):
        file.write(str(Rs_point[i])+"\t"+str(Zs_point[i])+"\n")
    file.close()

def write_end_points():
    R = core.loader.RZ[:,0]
    Z = core.loader.RZ[:,1]
    Rmin = np.min(R)
    Rmax = np.max(R)
    Zmin = np.min(Z)
    Zmax = np.max(Z)
    file = open("End points for 'None' input.txt",'a')
    file.write("Rmin = "+"\t"+ str(Rmin)+"\n")
    file.write("Rmax = "+"\t"+ str(Rmax)+"\n")
    file.write("Zmin = "+"\t"+ str(Zmin)+"\n")
    file.write("Zmax = "+"\t"+ str(Zmax)+"\n")
    file.close()
#write_end_points()
#write_sep()
#Read in the values
#ne_all = np.load('ne_all_midplane.npy')
#pot_all = np.load('pot_all_midplane.npy')

#writes out all the matrices for postprocessing.
def write_mat():
    B_R = np.array([core.getcutvalue_hor(core.bfield[:,0],cut,1)for cut in range(0,len(L))])
    B_Z = np.array([core.getcutvalue_hor(core.bfield[:,1],cut,1)for cut in range(0,len(L))])
    B_zeta = np.array([core.getcutvalue_hor(core.bfield[:,2],cut,1)for cut in range(0,len(L))])
    ne_all = np.array([[[core.getcutvalue_hor(core.ne[:,iz,it],cut,1) for iz in range(0,core.Nplanes)]for cut in range(0,len(L))] for it in range(0,889)])
    #pot_all = np.array([[[core.getcutvalue_hor(core.pot[:,iz,it],cut,1) for iz in range(0,core.Nplanes)]for cut in range(0,len(L))] for it in range(490,690)])
    np.save(trgt + 'B_R_mid_CMOD', B_R)
    np.save(trgt + 'B_Z_mid_CMOD', B_Z)
    np.save(trgt + 'B_zeta_mid_CMOD', B_zeta)
    np.save(trgt + 'ne_mid_CMOD', ne_all)
    #np.save('pot_all_full', pot_all)

def write_turb():
    pot_all = np.array([[[core.getcutvalue_hor(core.pot[:,iz,it],cut,1) for iz in range(0,core.Nplanes)]for cut in range(0,len(L))] for it in range(490,690)])
    np.save('pot_turb_region', pot_all)

write_mat()

def den_contour_plot(time):
    ne_new = ne_all[time,:,:,:]
    ne_ta = ne_new.mean(axis=1)
    fig, ax =plt.subplots()
    CS = plt.contourf(R_points, Z_points, ne_ta)
    CB = plt.colorbar(CS, extend='both')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Density Contours, t = %d'%(time))
    plt.show()

#TO DO!    
def psi_cont():
    Rmin = core.loader.Rmin
    Rmax = core.loader.Rmax
    Zmin = core.loader.Zmin
    Zmax = core.loader.Zmax
    unitR = (Rmax-Rmin)/100
    unitZ = (Zmax-Zmin)/100
    R_points = [ip*unitR+Rmin for ip in range(100)]
    Z_points = [ip*unitZ+Zmin for ip in range(100)]
    psi=[]
    for cut in range(100):
        psi.append(core.getcutvalue_hor(core.psin[:],cut,1))
    psi = np.reshape(psi,(100,100))
    fig, ax =plt.subplots()
    contour_levels = np.arange(0,1.7,0.05)
    CS = plt.contour(R_points, Z_points, psi,contour_levels)
    CB = plt.colorbar(CS, extend='both')
    #plt.clabel(CS, inline=1, fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.axhline(y=0.0,color='k',linestyle='--')
    plt.axvline(x=core.Rmaj,color='k',linestyle='--')
    plt.title(r'$\Psi$ contours')
    plt.show()
    
    
    
    
def av_den_contour_plot():
    ne_tor = ne_all.mean(axis=2)
    ne_avg = ne_tor.mean(axis=0)
    fig, ax =plt.subplots()
    CS = plt.contourf(R_points, Z_points, ne_avg)
    CB = plt.colorbar(CS, extend='both')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Averaged Density Contours')
    plt.show()

def av_dden_contour_plot():
    ne_tor = ne_all.mean(axis=2)
    ne_avg = ne_tor.mean(axis=0)
    dn_all = np.array(ne_all[:,:,:,:]-ne_avg[:,:,np.newaxis,np.newaxis])
    dn_tor = dn_all.mean(axis=2)
    dn_avg = dn_tor(axis=0)
    fig, ax =plt.subplots()
    CS = plt.contourf(R_points, Z_points, dn_avg)
    CB = plt.colorbar(CS, extend='both')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Time Averaged Density Perturbation Contours')
    plt.show()


#av_den_contour_plot()

def pot_contour_plot(time):
    pot_new = pot_all[time,:,:,:]
    pot_ta = pot_all.mean(axis=1)
    fig, ax =plt.subplots()
    CS = plt.contourf(R_points, Z_points, pot_ta)
    CB = plt.colorbar(CS, extend='both')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Potential Contours, t = %d'%(time))
    plt.show()

def av_pot_contour_plot():
    pot_tor = pot_all.mean(axis=2)
    pot_avg = pot_tor.mean(axis=0)
    fig, ax =plt.subplots()
    CS = plt.contourf(R_points, Z_points, pot_avg)
    CB = plt.colorbar(CS, extend='both')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Averaged Potential Contours')
    plt.show()


#av_pot_contour_plot()


def pot0_contour_plot(time):
    pot0_all = np.array([core.getcutvalue_hor(core.pot0[:,time],cut,1)for cut in range(0,100)])
    fig, ax =plt.subplots()
    CS = plt.contourf(R_points, Z_points, pot0_all)
    CB = plt.colorbar(CS, extend='both')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Pot0 Contours, t=%d'%(time))
    plt.show()

def av_pot0_contour_plot():
    pot0_all = np.array([[core.getcutvalue_hor(core.pot0[:,time],cut,1)for cut in range(0,100)]for time in range(500,650)])
    pot0_av = pot0_all.mean(axis=0)
    fig, ax =plt.subplots()
    CS = plt.contourf(R_points, Z_points, pot0_av)
    CB = plt.colorbar(CS, extend='both')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Averaged Pot0 Contours')
    plt.show()


#av_pot0_contour_plot()


def overlay_contours(time):
    ne_new = ne_all[time,:,:,:]
    pot_new = pot_all[time,:,:,:]
    ne_ta = ne_new.mean(axis=1)
    pot_ta = pot_new.mean(axis=1)
    fig =plt.figure()
    CS = plt.contour(R_points, Z_points, ne_ta,colors='b')
    CS2 = plt.contour(R_points, Z_points, pot_ta,colors='r')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.clabel(CS2, inline=1,fmt='%.2E', fontsize=10)
    labels = ['$n_e$','$\phi$']
    CS.collections[0].set_label('$n_e$')
    CS2.collections[0].set_label('$\phi$')
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    plt.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Overlayed Contours of $n_e$ and $\phi$')
    plt.legend()
    plt.show()


def av_overlay_contours():
    ne_tor = ne_all.mean(axis=2)
    ne_avg = ne_tor.mean(axis=0)
    pot_tor = pot_all.mean(axis=2)
    pot_avg = pot_tor.mean(axis=0)
    fig =plt.figure()
    CS = plt.contour(R_points, Z_points, ne_avg,colors='b')
    CS2 = plt.contour(R_points, Z_points, pot_avg,colors='r')
    plt.clabel(CS, inline=1,fmt='%.2E', fontsize=10)
    plt.clabel(CS2, inline=1,fmt='%.2E', fontsize=10)
    labels = ['$n_e$','$\phi$']
    CS.collections[0].set_label('$n_e$')
    CS2.collections[0].set_label('$\phi$')
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    plt.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Overlayed Contours of averaged $n_e$ and $\phi$')
    plt.legend()
    plt.show()


#av_overlay_contours()

def animation_contour(time_start,time_end):
    fig,ax = plt.subplots()
    ims=[]
    imi=[]
    for time in range(time_start,time_end):
        pot_init = pot_all[time,:,:,:]
        pot_ta = pot_init.mean(axis=1)
        im = plt.contourf(R_points, Z_points, pot_ta)
        imi.append(im)
        #################################################################
    ## Bug fix for Quad Contour set not having attribute 'set_visible', 'set_animated'
        def setvisible(self,vis):
            for c in self.collections: c.set_visible(vis)
        def setanimated(self,ani):
            for c in self.collections: c.set_animated(ani)
        im.set_visible = types.MethodType(setvisible,im)
        im.set_animated = types.MethodType(setanimated,im)
        im.axes = plt.gca()
        im.figure=fig
    ####################################################################
        ims.append([im])
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)    
    
    im_ani = animation.ArtistAnimation(fig, ims, interval=200, repeat_delay=1000,blit=False)
    #plt.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.colorbar(imi[-1], extend='both')
    plt.xlabel('R(m)')
    plt.ylabel("Z(m)")
    plt.title('Potential Contours between t = %d - %d'%(time_start,time_end))
    plt.show()
    im_ani.save('potential.mp4')

#animation_contour(0,200)
