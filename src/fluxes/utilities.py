import numpy as np
import xgc
import math
#import xgcjrm as xgc
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
import sys
from decimal import Decimal


def getcutvalue(arrRZ,j=None):
	"""Calculates the values at Ri along Zcut of arrRZ where arr is a mesh array
	usage:
		netmpi = getcutvalue(ne[:,iz,it]) # where iz is the toroidal plane, and it the time index
	"""
	global triObj,RI,ZI,arrout,tci,jcut,RI,ZI
	# the following construction allows the def to occur before jcut takes a value
	if j==None:
		j = jcut 
	tci=LinearTriInterpolator(triObj,arrRZ)
	arrout=tci(RI,ZI)
	return(arrout[j,:])



def getcutvRvZ_old(iz,it,cut):
    global Zi,Ri,pot,j,bfield,B2
    if cut==None:
	    cut = jcut	
    dZ = Zi[1]-Zi[0]
    dR = Ri[1]-Ri[0]
    theta = (2*math.pi)/Nplanes
    #calculate values at the cut
    br = getcutvalue(bfield[:,0])
    bz = getcutvalue(bfield[:,1])
    bzeta = getcutvalue(bfield[:,2])
    BfMag=[]
    for i in range(len(br)):
        BfMag.append(np.sqrt(br[i]*br[i] + bz[i]*bz[i] + bzeta[i]*bzeta[i]))
#Calculating the potential gradient
# dphi/dZ    
    potim1 = getcutvalue(pot[:,iz,it],cut-1)
    poti = getcutvalue(pot[:,iz,it],cut)
    potip1 = getcutvalue(pot[:,iz,it],cut+1)
    dphidZ = (potip1-potim1)/(2.*dZ)
# dphi/dR
    dphidR=[]
    lengthR = len(getcutvalue(pot[:,iz,it],cut))
    iterlistR = list(range(0,lengthR))
    potrm10 = getcutvalue(pot[:,iz,it],cut)[lengthR-1]
    potrp10 = getcutvalue(pot[:,iz,it],cut)[1]
    dphidR.append((potrp10 - potrm10)/(2.*dR))
    for i in iterlistR[1:-1]:
        potrm1 = getcutvalue(pot[:,iz,it],cut)[i-1]
        potrp1 = getcutvalue(pot[:,iz,it],cut)[i+1]
        dphidR.append((potrp1 - potrm1)/(2.*dR))
    potrm1L = getcutvalue(pot[:,iz,it],cut)[lengthR-2]
    potrp1L = getcutvalue(pot[:,iz,it],cut)[0]
    dphidR.append((potrp1L-potrm1L)/2*dR)
#dphi/dzeta
    if iz < (Nplanes-1):	
        potzp1 = getcutvalue(pot[:,iz+1,it],cut) 
        potz = getcutvalue(pot[:,iz,it],cut)
        dphidzeta=[]
        lengthZeta = len(Ri)
        iterlistZeta = list(range(0,lengthZeta))
        for i in iterlistZeta:
            dphidzeta.append((potzp1[i] - potz[i])/(Ri[i]*theta))
    else:
        potzp1 = getcutvalue(pot[:,0,it],cut)
        potz = getcutvalue(pot[:,iz,it],cut)
        dphidzeta=[]
        lengthZeta = len(Ri)
        iterlistZeta = list(range(0,lengthZeta))
        for i in iterlistZeta:
            dphidzeta.append((potzp1[i] - potz[i])/(Ri[i]*theta))

#Calculating the ExB drift
    vR = []
    lengthdz=len(dphidZ)
    iterlistdz = list(range(0,lengthdz))
    for i in iterlistdz:
        vR.append(dphidZ[i]*bzeta[i] - dphidzeta[i]*bz[i])
    vZ = []
    for i in iterlistR:
        vZ.append(br[i]*dphidzeta[i] - dphidR[i]*bzeta[i])
    return(vR,vZ)


def getcutvRvZ(iz,it,cut):
	global Zi,Ri,pot,jcut,bfield,B2
	
	if cut==None:
		cut = jcut

	dZ = Zi[1]-Zi[0]
	dR = Ri[1]-Ri[0]

	potim1 = getcutvalue(pot[:,iz,it],cut-1)
	poti = getcutvalue(pot[:,iz,it],cut)
	potip1 = getcutvalue(pot[:,iz,it],cut+1)
	dphidZ = (potip1-potim1)/(2.*dZ)

	potiplus = np.roll(poti,-1)
	potiminus = np.roll(poti,1)
	dphidR = (potiplus-potiminus)/(2.*dR)
	dphidR[0] = (poti[1]-poti[0])/dR
	dphidR[-1] = (poti[-1]-poti[-2])/dR

	# calculate B-field = (BR,BZ,Bphi) components at the cut
	BRi = getcutvalue(bfield[:,0],cut)
	BZi = getcutvalue(bfield[:,1],cut)
	Bzetai = getcutvalue(bfield[:,2],cut)
	# save some useful combos
	Bzobzet2 = (BZi/Bzetai)**2
	Brobzet2 = (BRi/Bzetai)**2
	BrBzobzet2 = (BRi*BZi)/(Bzetai**2)
	B2 = BRi**2 + BZi**2 + Bzetai**2

	vRi =  (Bzetai/B2)*(dphidZ*(1.+Bzobzet2)+BrBzobzet2*dphidR)
	vZi = -(Bzetai/B2)*(dphidR*(1.+Brobzet2)+BrBzobzet2*dphidZ)
	return(vRi,vZi)



def getiparrs(psimin,psimax,cut,verbose=False): 
    global chc,psiichc,Richc
    if cut == None:
        cut = jcut
    tci=LinearTriInterpolator(triObj,psin)
    psiout=tci(RI,ZI)
    psii = psiout[cut,:]
    chc1 = np.where(psii > psimin)[0]
    chc2 = np.where(psii < psimax)[0]
    chc = np.intersect1d(chc1,chc2)
    psiichc = psii[chc].filled()
    Richc = Ri[chc]
    if verbose:
	for i in range(len(chc)):
            print ("%d   %4f   %4f" % (chc[i],psiichc[i],Richc[i]))
    return chc

getiparrs(0.95,1.05,None)

def VRandCCplots():
    tfin = pot.shape[2]-1
    time_point = int(input("Enter the time point (Value must between 0 and %s): "%tfin))
    #Define the max and min value of flux
    (psimin,psimax) = (0.95, 1.05)





    VR_all = np.array([getcutvRvZ(iz,time_point,jcut)[0] for iz in range(Nplanes)])
    print('VR_all.shape=', VR_all.shape)

    n_e=loader.calcNeTotal()
    print('n_e.shape=', n_e.shape)
    n_e_all = np.array([getcutvalue(n_e[:,iz,time_point],jcut) for iz in range(Nplanes)])
    print('n_e_all.shape', n_e_all.shape)


    VR_tor_avg = VR_all.mean(axis=0)
    print('VR_tor_avg =', VR_tor_avg.shape)
    VR_in = np.array([[VR_all[iz,ip]-VR_tor_avg[ip] for iz in range(Nplanes)] for ip in getiparrs(psimin,1.,None)])
    print('VR_in=', VR_in.shape)
    n_e_tor_avg = n_e_all.mean(axis=0)
    print('ne_tor_avg=', n_e_tor_avg.shape)
    dn_ov_n_in = np.array([[(n_e_all[iz,ip]-n_e_tor_avg[ip])/n_e_tor_avg[ip] for iz in range(Nplanes)] for ip in getiparrs(psimin,1.,None)])
    print('dn_ov_n_in=', dn_ov_n_in.shape)
    VR_out = np.array([[VR_all[iz,ip]-VR_tor_avg[ip] for iz in range(Nplanes)] for ip in getiparrs(1.,psimax,None)])
    print('VR_out=', VR_out.shape)
    dn_ov_n_out = np.array([[(n_e_all[iz,ip]-n_e_tor_avg[ip])/n_e_tor_avg[ip] for iz in range(Nplanes)] for ip in getiparrs(1.,psimax,None)])
    print('dn_ov_n_out=', dn_ov_n_out.shape)
    Vmaskin = np.isfinite(VR_in)
    Vmaskout = np.isfinite(VR_out)
    dnmaskin = np.isfinite(dn_ov_n_in)
    dnmaskout = np.isfinite(dn_ov_n_out)



    plt.title('$V_r$ $vs.$ $\delta n/n$ $inside$ $separatrix,$ $t=%s$'%(time_point))
    plt.xlabel('$\delta n/n$')
    plt.ylabel('$V_r$ $(m/s)$')
    plt.grid(True)
    plt.plot(dn_ov_n_in[dnmaskin][:],VR_in[Vmaskin][:],'go')
    plt.show()

    plt.title('file = %s \n $V_r$ $vs.$ $\delta n/n$ $outside,$ $t=%s$'%(fileDir,time_point))
    plt.xlabel('$\delta n/n$')
    plt.ylabel('$V_r$')
    plt.grid(True)
    plt.plot(dn_ov_n_out[dnmaskout][:],VR_out[Vmaskout][:],'go')
    plt.show()


def ccornv(ip):
    global dneone,dvRi
    dneone = (n_e_all[:,ip]-n_e_tor_avg[ip])/n_e_tor_avg[ip]
    dvRi = VR_all[:,ip]-VR_tor_avg[ip]
    return np.dot(dneone,dvRi)/(dvRi.shape[0]*np.sqrt(dneone.std()**2*dvRi.std()**2))

chc_list = getiparrs(psimin,psimax,None)
ccornvarr = [ccornv(ip) for ip in chc_list]
plt.title('$Cross-Correlation$ $at$ $t=%s$'%(time_point))
plt.xlabel('R(m)')
plt.ylabel('C')
plt.plot(Richc, ccornvarr,'b.-')
plt.show()



def timehistory():
    n_e=loader.calcNeTotal()
    (psimin,psimax) = (0.95, 1.05)
    tfin = pot.shape[2]-1
    n_e_all_all = np.array([[getcutvalue(n_e[:,iz,it],jcut)for iz in range(Nplanes)]for it in range(ne.shape[2])])
    n_e_all_all_mean = n_e_all_all.mean(axis=1)
    point1 = getiparrs(psimin,1.0,None)[0]
    Sep_point = getiparrs(psimin,1.0,None)[-1]
    pointL = getiparrs(1.0,psimax,None)[-1]
    point = int(input("Select point to plot time history between %s and %s. The Separatrix is located at %s: "%(point1,pointL,Sep_point)))
    pick_plane = int(input("Select a plane for the time history: "))
    n_e_all_2 = n_e_all_all[:,pick_plane,point] - n_e_all_all_mean[:,point]
    times_list = range(0,tfin+1)
    plt.title('file = %s \n $time$ $history$ $of$ $point$ $%s$ $at$ $plane$ $%s$'%(fileDir,point,pick_plane))
    plt.xlabel('time')
    plt.ylabel('density')
    plt.plot(times_list,n_e_all_2)
    plt.show()


def twodtimeplot():
    n_e=loader.calcNeTotal()
    n_e_all_all = np.array([[getcutvalue(n_e[:,iz,it],jcut)for iz in range(Nplanes)]for it in range(ne.shape[2])])
    n_e_all_all_mean = n_e_all_all.mean(axis=1)
    #z = np.array([[((n_e_all_all[it,iz,1:99] - n_e_all_all_mean[it,iz])/n_e_all_all_mean[it,iz]) for iz in range(Nplanes)]for it in range(ne.shape[2]-1)])
    for i in range(Nplanes):
        dmove = 1
        timez, zpace = np.mgrid[slice(0,ne.shape[2]-1,dmove),slice(1,99,dmove)]
        z = (n_e_all_all[timez,i,zpace]-n_e_all_all_mean[timez,zpace])/n_e_all_all_mean[timez,zpace]
       # z_min, z_max = np.abs(z2).min(), np.abs(z2).max()
       # im = plt.pcolor(zpace, timez, z, cmap='RdBu', vmin=z_min, vmax=z_max)
       # plt.title('$\delta n /n$ $vs.$ $time$ $and$ $space$ $for$ $plane$ $%s$'%i)
       # plt.xlabel('space')
       # plt.ylabel('time')
       # plt.colorbar(im)
       # plt.show()
        z2 = n_e_all_all[timez,i,zpace]
        z_min, z_max = np.abs(z2).min(), np.abs(z2).max()
        im2 = plt.pcolor(zpace, timez, z2, cmap='RdBu', vmin=z_min, vmax=z_max)
        plt.title('$Density$ $vs.$ $time$ $and$ $space$ $for$ $plane$ $%s$'%i)
        plt.xlabel('space')
        plt.ylabel('time')
        plt.colorbar(im2)
        plt.show()
       
def RZplot():
    tfin = pot.shape[2]-1
    ud = int(input('Enter how many lines up/down from jcut:'))
    t_start = int(input('Enter starting time:' ))
    t_stop = int(input('Enter stopping time (up to %s):'%tfin))
    choice = int(input('Choose between dn/n (0) and std (1) values to be plotted:'))
    n_e=loader.calcNeTotal()
    ne_RZ = np.array([[[getcutvalue(n_e[:,iz,it],newjay) for iz in range(Nplanes)]for it in range(ne.shape[2])]for newjay in range(jcut-ud,jcut+ud)])
    ne_RZ_P = ne_RZ[:,:,:,3:97]
    ne_RZ_std = ne_RZ_P.std(axis=2)
    ne_RZ_mean = ne_RZ_P.mean(axis=2)
    dmove = 1
    zetas, Rs = np.mgrid[slice(0,ne_RZ.shape[0]-1,dmove),slice(0,94,dmove)]
    for i in range(t_start,t_stop,1):
        if choice == 1:
            z1 = ((ne_RZ_P[zetas, i, 0, Rs]-ne_RZ_mean[zetas,i,Rs])/ne_RZ_mean[zetas,i,Rs])
            z_min, z_max = np.abs(z1).min(), np.abs(z1).max()
            im = plt.pcolor(Rs, zetas, z1, cmap='RdBu', vmin=z_min, vmax=z_max)
        else:
            z_min, z_max = np.abs(ne_RZ_std).min(), np.abs(ne_RZ_std).max()
            im = plt.pcolor(Rs, zetas, ne_RZ_std, cmap='RdBu', vmin=z_min, vmax=z_max)
        plt.colorbar(im)
        plt.title('Density at time=%s'%i)
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.show()
        im.figure.savefig('Density at time=%s'%i)

def RZ_inst_plot():
    tfin = pot.shape[2]-1
    t_point = int(input('Enter the time point of the plot (up to %s):'%tfin))
    choice = int(input('Choose between dn/n (0) and std (1) values to be plotted:'))
    n_e=loader.calcNeTotal()
    ne_RZ = np.array([[getcutvalue(n_e[:,iz,t_point],newjay) for iz in range(Nplanes)]for newjay in range(jcut-120,jcut+140)])
    ne_RZ_P = ne_RZ[:,:,60:150]
    ne_RZ_std = ne_RZ_P.std(axis=1)
    ne_RZ_mean = ne_RZ_P.mean(axis=1)
    dmove = 1
    zetas, Rs = np.mgrid[slice(0,ne_RZ.shape[0],dmove),slice(0,90,dmove)]
    if choice == 0:
        z1 = ((ne_RZ_P[zetas, 0,  Rs]-ne_RZ_mean[zetas, Rs])/ne_RZ_mean[zetas,Rs])
        z_min, z_max = np.abs(z1).min(), np.abs(z1).max()
        im = plt.pcolor(Rs, zetas, z1, cmap='RdBu', vmin=z_min, vmax=z_max)
        plt.title('$\delta n/n$ density at time=%s, whole machine'%t_point)
        plt.colorbar(im)
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.show()
        im.figure.savefig('Dn/n density at time=%s, whole machine'%t_point)
    else:
        z_min, z_max = np.abs(ne_RZ_std).min(), np.abs(ne_RZ_std).max()
        im = plt.pcolor(Rs, zetas, ne_RZ_std, cmap='RdBu', vmin=z_min, vmax=z_max)
        plt.title('RMS density at time=%s, whole machine'%t_point)
        plt.colorbar(im)
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.show()
        im.figure.savefig('RMS density at time=%s, whole machine'%t_point)



def load_flux_data(cut):
    global VR_all, ne_all, VZ_all, br, bz, triObj,RI,ZI,tci,RI,ZI,Sep_ip
    #j_val = int(input('Enter the number of the j-cut value:'%j_val))
    VR_all = np.array([[getcutvRvZ(iz,it,cut)[0] for iz in range(Nplanes)] for it in range(220,230)])
    VZ_all = np.array([[getcutvRvZ(iz,it,cut)[1] for iz in range(Nplanes)] for it in range(220,230)])
    n_e=loader.calcNeTotal()
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],cut) for iz in range(Nplanes)] for it in range(220,230)])
    br = getcutvalue(bfield[:,0],cut)
    bz = getcutvalue(bfield[:,1],cut)
    Sep_ip = getiparrs(0,95,1.,cut)[-1]
   

def flux_plot():
    ne_tor_avg = ne_all.mean(axis=1)
    VR_tor_avg = VR_all.mean(axis=1)
    VZ_tor_avg = VZ_all.mean(axis=1)
    dn = np.array(ne_all[:,:,:] - ne_tor_avg[:,np.newaxis,:])
    dVR = np.array(VR_all[:,:,:] - VR_tor_avg[:,np.newaxis,:])
    dVZ = np.array(VZ_all[:,:,:] - VZ_tor_avg[:,np.newaxis,:])
    R_flux = np.array([[[dn[it,pl,ip]*dVR[it,pl,ip] for it in range(0,10)] for pl in range(Nplanes)] for ip in range(0,100)])
    Z_flux = np.array([[[dn[it,pl,ip]*dVZ[it,pl,ip] for it in range(0,10)] for pl in range(Nplanes)] for ip in range(0,100)])
    R_flux_tor_avg = R_flux.mean(axis=1)
    Z_flux_tor_avg = Z_flux.mean(axis=1)
    R_flux_time_avg = R_flux_tor_avg.mean(axis=1)
    Z_flux_time_avg = Z_flux_tor_avg.mean(axis=1)
    L = [x for x in range(0,100)]
    unit = (Rmax-Rmin)/100
    R = np.array([(L[ip]*unit + Rmin) for ip in range(0,100)])
    fig = plt.figure()
    plt.title('R-flux time average vs. radial position, cut = 46')
    plt.xlabel('R(m)')
    plt.ylabel('R-flux')
    plt.grid(True)
    plt.plot(R[:],R_flux_time_avg[:],'b.-')
    plt.axvline(x=Sep_ip*unit +Rmin, color='b', linestyle='dashed', linewidth=2)
    plt.show()
    fig.savefig('R-flux time average vs. rad. position, cut = 46.png')
    fig = plt.figure()
    plt.title('Z-flux time average vs. radial position, cut = 46')
    plt.xlabel('R (m)')
    plt.ylabel('Z_flux')
    plt.grid(True)
    plt.plot(R[:],Z_flux_time_avg[:],'b.-')
    plt.axvline(x=Sep_ip*unit + Rmin, color='b', linestyle='dashed', linewidth=2)
    plt.show()
    fig.savefig('Z-flux time average vs. rad. position, cut=46.png')
    Bmag=np.array([np.sqrt(br[ip]*br[ip] + bz[ip]*bz[ip])for ip in range(0,100)])
    Gamma = np.array([(-br[ip]*Z_flux_tor_avg[ip,:] + bz[ip]*R_flux_tor_avg[ip,:])/Bmag[ip] for ip in range(0,100)])
    Gamma_time_avg = Gamma.mean(axis=1)
    fig = plt.figure()
    plt.title('Normal flux time average vs. radial position, cut =46')
    plt.xlabel('R (m)')
    plt.ylabel('Normal flux')
    plt.grid(True)
    plt.plot(R[:],Gamma_time_avg[:],'b.-')
    plt.axvline(x=Sep_ip*unit + Rmin, color='b', linestyle='dashed', linewidth=2)
    plt.show()
    fig.savefig('N-flux time average vs. rad. position, jcut=46.png')

def load_flux_data_2():#Change the numbers before running it, Check the thing with the separatrix line in above function
    global VR_all, ne_all, VZ_all, BRi, Bzi, Bzetai,  triObj,RI,ZI,tci,RI,ZI,Sep_ip
    n_e=loader.calcNeTotal()
    angle_Gamma = [] 
    angle_list = []
    file = open("mean flux vs poloidal angle_new_4.txt",'a')
    file.write("pol. angle"+"\t"+"mean flux"+"\t"+"std flux"+"\n")
    unit_R = (Rmax-Rmin)/100
    unit_Z = (Zmax-Zmin)/100
    for cut in range(32,80):
	BRi = getcutvalue(bfield[:,0])
	BZi = getcutvalue(bfield[:,1])
	Bzetai = getcutvalue(bfield[:,2])
        Z_val = (cut-jcut)*unit_Z  
        Flux_ip = getiparrs(0.95,1.,cut)[-2]
	print("Flux ip = ",Flux_ip)
        R_val = Flux_ip*unit_R + Rmin
        tangent = (Z_val/R_val)
        angle = np.degrees(np.arctan(tangent))
        angle_list.append(angle)
	print("angle = ",angle)
        VR_all = np.array([[getcutvRvZ(iz,it,cut)[0][Flux_ip] for iz in range(Nplanes)] for it in range(200,320)])
	print("VR_all done..")
        VZ_all = np.array([[getcutvRvZ(iz,it,cut)[1][Flux_ip] for iz in range(Nplanes)] for it in range(200,320)])        
	print("VZ_all done..")
        ne_all = np.array([[getcutvalue(n_e[:,iz,it],cut)[Flux_ip] for iz in range(Nplanes)] for it in range(200,320)])
	print("ne_all done..")
       # br = getcutvalue(bfield[:,0],cut)
       # bz = getcutvalue(bfield[:,1],cut)
        ne_tor_avg = ne_all.mean(axis=1)
        VR_tor_avg = VR_all.mean(axis=1)
        VZ_tor_avg = VZ_all.mean(axis=1)
        dn = np.array(ne_all[:,:] - ne_tor_avg[:,np.newaxis])
        dVR = np.array(VR_all[:,:] - VR_tor_avg[:,np.newaxis])
        dVZ = np.array(VZ_all[:,:] - VZ_tor_avg[:,np.newaxis])
        R_flux = np.array([[dn[it,pl]*dVR[it,pl] for it in range(0,120)] for pl in range(Nplanes)])
        Z_flux = np.array([[dn[it,pl]*dVZ[it,pl] for it in range(0,120)] for pl in range(Nplanes)])
        R_flux_tor_avg = R_flux.mean(axis=1)
        Z_flux_tor_avg = Z_flux.mean(axis=1)
        R_flux_time_avg = R_flux_tor_avg.mean()
        Z_flux_time_avg = Z_flux_tor_avg.mean()
        Bmag=np.array(np.sqrt(BRi[:]*BRi[:] + BZi[:]*BZi[:]))
        Gamma = np.array((-BRi[Flux_ip]*Z_flux_tor_avg[:] + BZi[Flux_ip]*R_flux_tor_avg[:])/Bmag[Flux_ip])
        Gamma_std = Gamma.std()
        Gamma_time_avg = Gamma.mean()
        angle_Gamma.append(Gamma_time_avg)
        file.write(str(angle)+"\t"+str(Gamma_time_avg)+"\t"+str(Gamma_std)+"\n")
    AN_G = np.asarray(angle_Gamma)
    mean_gamma = AN_G.mean()
    print(mean_gamma)
    file.write("averaged flux:"+ str(mean_gamma)+"\n")
    file.close()
    fig = plt.figure()
    plt.title('mean flux vs. poloidal angle')
    plt.xlabel('poloidal angle')
    plt.ylabel('mean flux')
    plt.grid(True)
    plt.plot(angle_list,angle_Gamma,'ro')    
    plt.show()


def flux_run_av():
    global VR_all, ne_all, VZ_all, BRi, BZi, triObj, RI, ZI, tci, Sep_ip
    n_e=loader.calcNeTotal()
    angle_Gamma = [] 
    angle_list = []
    file = open("running average midplane (10tpts).txt",'a')
    file.write("pol. angle"+"\t"+"mean flux"+"\n")
    unit_R = (Rmax-Rmin)/100
    unit_Z = (Zmax-Zmin)/100
    for cut in range(5,78):
        Z_val = (cut-jcut)*unit_Z #+ Zmin 
        Flux_ip = getiparrs(0.95,1.,cut)[-2]
        print('Flux_ip=',Flux_ip)
        R_val = Flux_ip*unit_R + Rmin
        tangent = (Z_val/R_val)
        angle = np.degrees(np.arctan(tangent))
        angle_list.append(angle)
        VR_all = np.array([[getcutvRvZ(iz,it,cut)[0][Flux_ip] for iz in range(Nplanes)] for it in range(200,320)])
        print("VR_all done...")
        VZ_all = np.array([[getcutvRvZ(iz,it,cut)[1][Flux_ip] for iz in range(Nplanes)] for it in range(200,320)])
        print("VZ_all done...")
        ne_all = np.array([[getcutvalue(n_e[:,iz,it],cut)[Flux_ip] for iz in range(Nplanes)] for it in range(200,320)])
        print("ne_all done...")
        BRi = getcutvalue(bfield[:,0],cut)
        BZi = getcutvalue(bfield[:,1],cut)
        ne_tor_avg = ne_all.mean(axis=1)
        VR_tor_avg = VR_all.mean(axis=1)
        VZ_tor_avg = VZ_all.mean(axis=1)
        ne_red = ne_all[10:-10,:]#reduce by number of values left out in the for loop of the time average
        VR_red = VR_all[10:-10,:]
        VZ_red = VZ_all[10:-10,:]
        ne_m = []
        VR_m = []
        VZ_m = []
        for i in range(10,110):#range(x,y), ne_tor_avg[i-a,i+b] => x=a, y+b=tot_time_points
            ne_m.append(np.mean(ne_tor_avg[i-10:i+10]))
            VR_m.append(np.mean(VR_tor_avg[i-10:i+10]))
            VZ_m.append(np.mean(VZ_tor_avg[i-10:i+10]))
        ne_avg = np.asarray(ne_m)
        VR_avg = np.asarray(VR_m)
        VZ_avg = np.asarray(VZ_m)
        dn = np.array(ne_red[:,:] - ne_avg[:,np.newaxis])
        dVR = np.array(VR_red[:,:] - VR_avg[:,np.newaxis])
        dVZ = np.array(VZ_red[:,:] - VZ_avg[:,np.newaxis])
        R_flux = np.array([[dn[it,pl]*dVR[it,pl] for it in range(0,100)] for pl in range(Nplanes)])#put number of reduced time elements
        Z_flux = np.array([[dn[it,pl]*dVZ[it,pl] for it in range(0,100)] for pl in range(Nplanes)])
        R_flux_tor_avg = R_flux.mean(axis=1)
        Z_flux_tor_avg = Z_flux.mean(axis=1)
        R_flux_time_avg = R_flux_tor_avg.mean()
        Z_flux_time_avg = Z_flux_tor_avg.mean()
        Bmag=np.array(np.sqrt(BRi[:]*BRi[:] + BZi[:]*BZi[:]))
        Gamma = np.array((-BRi[Flux_ip]*Z_flux_tor_avg[:] + BZi[Flux_ip]*R_flux_tor_avg[:])/Bmag[Flux_ip])
        Gamma_time_avg = Gamma.mean()
        angle_Gamma.append(Gamma_time_avg)
        file.write(str(angle)+"\t"+str(Gamma_time_avg)+"\n")
        print("cut %s done..."%cut)
    AN_G = np.asarray(angle_Gamma)
    mean_gamma = AN_G.mean()
    print(mean_gamma)
    file.write("averaged flux:"+ str(mean_gamma)+"\n")
    file.close()
    fig = plt.figure()
    plt.title('mean flux vs. poloidal angle')
    plt.xlabel('poloidal angle (degrees)')
    plt.ylabel('mean flux')
    plt.grid(True)
    plt.plot(angle_list,angle_Gamma,'ro')    
    plt.show()
    #fig.savefig('C:\\Users\\giannis\\Desktop\\mean_flux_vs_pol_angle')rint(Sep_ip)

def density_profile():
    n_e=loader.calcNeTotal()
    Sep_ip = getiparrs(0.95,1.,jcut)[-1]
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],jcut) for iz in range(Nplanes)] for it in range(0,320)])
    ne_tor_avg=ne_all.mean(axis=1)
    ne_avg = ne_tor_avg.mean(axis=0)
    data = ne_avg[Sep_ip:80]
    mean = data.mean()
    new_data = data/mean
    L = np.array([x for x in range(0,100)])
    unit = (Rmax-Rmin)/100
    R = np.array([L[ip]*unit+Rmin for ip in range(0,80-Sep_ip)])
    plt.title('Density')
    plt.xlabel('R')
    plt.ylabel('ne')
    plt.plot(L[:],ne_avg[:],'go')
    plt.show()

def fit_width():
    n_e=loader.calcNeTotal()
    Sep_ip = getiparrs(0.95,1.,jcut)[-1]
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],jcut) for iz in range(Nplanes)] for it in range(210,320)])
    ne_tor_avg=ne_all.mean(axis=1)
    ne_avg = ne_tor_avg.mean(axis=0)
    data = ne_avg[Sep_ip:78]
    mean = data.mean()
    new_data = data/mean
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    L = np.array([x for x in range(0,78-Sep_ip)])
    unit = (Rmax-Rmin)/100
    R = np.array([L[ip]*unit+Rmin for ip in range(0,78-Sep_ip)])
    popt, pcov = curve_fit(func, L, new_data, p0=(1,1,1), maxfev = 20000)
    new_popt = [mean*popt[0]*np.exp((popt[1]*Rmin)/unit), popt[1]/unit, popt[2]*mean]
    plt.title('Density Profile')
    plt.xlabel('R(m)')
    plt.ylabel('$n_e$')
    plt.text(2.211,1.3e+19,'SOL width = 4 mm')
    plt.plot(R[:], data[:],'go',label='Original Data')
    plt.plot(R[:], func(R, *new_popt), 'r-', label="Fitted Curve")
    plt.legend()
    plt.show()
    print('fit values =',new_popt)
    print('SOL width = ',1/new_popt[1])


def speed_of_sound():
    m_i = 1.6723e-27
    k = 1.38e-23
    joule_conv = 1.602e-19
    Flux_ip = getiparrs(0.95,1.,jcut)[-2]
    psi_val = getcutvalue(psin[:],jcut)[Flux_ip]
    oned_location = (np.where(psin1d>psi_val))[0][0]
    Te_avg = Te1d.mean(axis=0)
    T_e = Te_avg[oned_location-1]
    speed = np.sqrt((joule_conv*T_e)/(2*m_i))
    return speed   

def sol_width():
    n_e=loader.calcNeTotal()
    Flux_ip = getiparrs(0.95,1.,jcut)[-2]
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],jcut)[Flux_ip] for iz in range(Nplanes)] for it in range(210,320)])
    ne_tor_avg=ne_all.mean(axis=1) 
    ne_avg = ne_tor_avg.mean(axis=0)

    flux = 1.01856781613e+20
    tau = 13/speed_of_sound()
    print('speed of sound:',speed_of_sound())
    print('tau:',tau)
    print('ne_avg',ne_avg)
    w = (flux*tau)/ne_avg
    print('width:',w)


def angle_transform():
    #old_anglelist = []
    #new_anglelist = []
    file = open("new_angle_list.txt",'a')
    file.write("cut"+"\t"+"old angle"+"\t"+"new angle"+"\n")
    unit_R = (Rmax-Rmin)/100
    unit_Z = (Zmax-Zmin)/100
    Rmaj=1.67
    for cut in range(5,95):
        Z_val_new = (cut-jcut)*unit_Z #+ Zmin 
        #Z_val_old = (cut-jcut)*unit_Z + Zmin
        Flux_ip = getiparrs(0.95,1.,cut)[-2]
        #print('Flux_ip=',Flux_ip)
        R_val_old = Flux_ip*unit_R + Rmin
        R_val_new = Flux_ip*unit_R - Rmaj
        tangent_new = (Z_val_new/R_val_new)
        tangent_old = (Z_val_new/R_val_old)
        angle_new = np.degrees(np.arctan(tangent_new))
        #new_anglelist.append(angle_new)
        angle_old = np.degrees(np.arctan(tangent_old))
        #old_anglelist.append(angle_old)
        print('cut=',cut)
        print('old=',angle_old)
        print('new=',angle_new)                           
        file.write(str(cut)+"\t"+str(angle_old)+"\t"+str(angle_new)+"\n")
    file.close()

def density_map():
    n_e=loader.calcNeTotal()
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],jcut) for iz in range(Nplanes)] for it in range(0,320)])
    ne_std=ne_all.std(axis=1)
    dmove = 1
    Sep_ip = getiparrs(0.95,1.,jcut)[-1]
    unit = (Rmax-Rmin)/100
    timez, zpace = np.mgrid[slice(0,320,dmove),slice(1,99,dmove)]
    R = np.array(zpace*unit+Rmin)
    #print(R.shape)
    z = ne_std[timez,zpace]
    mz = np.ma.masked_invalid(z)
    #print(z[100,86:])
    z_min, z_max = np.abs(mz).min(), np.abs(mz).max()
    im = plt.pcolor(R[0][:], timez, mz, cmap='RdBu', vmin=z_min, vmax=z_max)
    plt.title('$Time$ $evolution$ $of$ $rms$ $density$ $at$ $midplane$')
    plt.axvline(x=Sep_ip*unit + Rmin, color='b', linestyle='dashed', linewidth=1)
    plt.gca().set_xlim(2.205,2.28)
    plt.gca().set_ylim(0,320)
    plt.xlabel('R(m)')
    plt.ylabel('time step')
    plt.colorbar(im)
    plt.show()



def fit_temp():
    Flux_ip = getiparrs(0.95,1.,jcut)[-2]
    psi_val = getcutvalue(psin[:],jcut)[Flux_ip]
    oned_location = (np.where(psin1d>psi_val))[0][0]
    Te_avg = Te1d.mean(axis=0)
    T_e = Te_avg[oned_location:]
    #mean = T_e.mean()
    #print('Te:',T_e.shape)
    #new_data = T_e/mean
    #def func(x, a, b, c):
        #return  a*np.exp(-(1/b) * x) + c
    L = np.array([x for x in range(0,11)])
    unit = (Rmax-Rmin)/100
    R = np.array([L[ip]*unit for ip in range(0,11)])
    #popt, pcov = curve_fit(func, R, T_e, p0=(200,750,100), maxfev = 20000)
    #new_popt = [300,popt[1]/unit, popt[2]*mean]
    #print('fit values =',popt)

    plt.title('Temperature Profile')
    plt.xlabel('R(m)')
    plt.ylabel('$n_e$')
    plt.plot(R[:],T_e[:],'go',label='Original Data')
    #plt.plot(R[:], func(R, *popt), 'r-', label="Fitted Curve")
    plt.plot(R[:],100*np.exp(-1800* R)+5,'r-',label="Fitted Curve")
    plt.legend()
    plt.show()


def av_vel():
    Flux_ip = getiparrs(0.95,1.,jcut)[-1]
    velocities = np.array([[getcutvRvZ(iz,it,jcut)[0][Flux_ip] for iz in range(Nplanes)] for it in range(220,320)])
    vel_tor_avg = velocities.mean(axis=1)
    vel_avg = vel_tor_avg.mean(axis=0)
    print('vel_avg:',vel_avg)


def velocity_prof():
    v_all = np.array([[getcutvRvZ(iz,it,jcut)[0] for iz in range(Nplanes)] for it in range(220,315)])
    v_time_avg = v_all.mean(axis=0)
    v_avg = v_time_avg.mean(axis=0)
    Sep_ip = getiparrs(0.95,1.,jcut)[-1]
    L = np.array([x for x in range(0,100)])
    unit = (Rmax-Rmin)/100
    R = np.array([L[ip]*unit+Rmin for ip in range(0,100)])
    plt.title('Velocity Profile')
    plt.xlabel('R(m)')
    plt.ylabel('$V_r$ (m/s)')
    plt.plot(R[:], v_avg[:],'go')
    plt.axvline(x=Sep_ip*unit + Rmin, color='b', linestyle='dashed', linewidth=2)
    plt.legend()
    plt.show()


def sep_boundaries():
    Sep_ip_L = []
    Sep_ip_R = []
    #Maximum and Minimum Values for when loading data as None
    Rmint = 1.0009998
    Rmaxt = 2.3769998
    Zmint = -1.363
    Zmaxt = 1.348000
    for cut in range(0,100):
        if getiparrs(0.95,1.,cut).shape[0]>0:
	    Sep_ip_L.append(getiparrs(0.95,1.,cut)[-1])
            Sep_ip_R.append(getiparrs(0.95,1.,cut)[0])
	else:
	    Sep_ip_L.append(0)
            Sep_ip_R.append(0)
    L = np.array([x for x in range(0,100)])
    unit_x = (Rmaxt-Rmint)/100
    unit_y = (Zmaxt-Zmint)/100
    R_L = np.array([Sep_ip_L[ip]*unit_x+Rmint for ip in range(0,100)])
    R_R = np.array([Sep_ip_R[ip]*unit_x+Rmint for ip in range(0,100)])
    cut = np.array([x for x in range(0,100)])
    Z = np.array([cut[ip]*unit_y+Zmint for ip in range(0,100)])
    file = open("Separatrix.txt",'a')
    for R in R_L:
        file.write(str(R)+"\n")
    file.write("R_L done..."+"\n")
    for R in R_R:
        file.write(str(R)+"\n")
    file.write("R_R done..."+"\n")
    for Zet in Z:
        file.write(str(Zet)+"\n")
    file.write("Z done....")
    file.close()
    plt.title('Separatrix Profile')
    plt.xlabel('$R(m)$')
    plt.ylabel('$Z(m)$')
    plt.grid(True)
    plt.plot(R_L[:], Z[:],'go', R_R[:], Z[:], 'bo')
    plt.show()

def time_avg_test():
    global VR_all, ne_all, VZ_all, br, bz, triObj,RI,ZI,tci,RI,ZI,Sep_ip,start
    #Inputs:
    time_range=85
    start=240 #start where linear phase ends
    ts_max=30
    file = open("mean flux vs time averaging jcut test.txt",'a')
    file.write("time interval"+"\t"+"mean flux"+"\n")
    Flux_ip = getiparrs(0.95,1.,jcut)[-2]
    print('Flux_ip=',Flux_ip)
    n_e=loader.calcNeTotal()
    temp = np.array([[getcutvRvZ(iz,it,cut) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    VR_all = temp[:,:,0,Flux_ip]
    VZ_all = temp[:,:,1,Flux_ip]

   # VR_all = np.array([[getcutvRvZ(iz,it,jcut)[0][Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    print("VR_all and VZ_all done...")
    #VZ_all = np.array([[getcutvRvZ(iz,it,jcut)[1][Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    #print("VZ_all done...")
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],jcut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    print("ne_all done...")
    br = getcutvalue(bfield[:,0],jcut)
    bz = getcutvalue(bfield[:,1],jcut)
    ne_tor_avg = ne_all.mean(axis=1)
    VR_tor_avg = VR_all.mean(axis=1)
    VZ_tor_avg = VZ_all.mean(axis=1)
    Gamma_list=[]
    ts_list=[]
    for ts in range(1,ts_max): 
    	l=1
    	m=0
        print('ts:',ts)
    	ne_m = []
    	VR_m = []
    	VZ_m = []
    	for i in range(0,ts_max):#left side
	    ne_m.append(np.mean(ne_tor_avg[i-m:i+ts+m]))
            VR_m.append(np.mean(VR_tor_avg[i-m:i+ts+m]))
            VZ_m.append(np.mean(VZ_tor_avg[i-m:i+ts+m]))
            m=m+1
    	for i in range(ts_max,time_range-ts_max):#center
            ne_m.append(np.mean(ne_tor_avg[i-ts:i+ts]))
            VR_m.append(np.mean(VR_tor_avg[i-ts:i+ts]))
            VZ_m.append(np.mean(VZ_tor_avg[i-ts:i+ts]))
    	for i in range(time_range-ts_max,time_range):#right side
	    ne_m.append(np.mean(ne_tor_avg[i-ts-l:i+ts-l]))
            VR_m.append(np.mean(VR_tor_avg[i-ts-l:i+ts-l]))
            VZ_m.append(np.mean(VZ_tor_avg[i-ts-l:i+ts-l]))
            l=l+1
    	ne_avg = np.asarray(ne_m)
    	VR_avg = np.asarray(VR_m)
    	VZ_avg = np.asarray(VZ_m)
    	dn = np.array(ne_all[:,:] - ne_avg[:,np.newaxis])
    	dVR = np.array(VR_all[:,:] - VR_avg[:,np.newaxis])
    	dVZ = np.array(VZ_all[:,:] - VZ_avg[:,np.newaxis])
    	R_flux = np.array([[dn[it,pl]*dVR[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    	Z_flux = np.array([[dn[it,pl]*dVZ[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    	R_flux_tor_avg = R_flux.mean(axis=1)
    	Z_flux_tor_avg = Z_flux.mean(axis=1)
    	R_flux_time_avg = R_flux_tor_avg.mean()
    	Z_flux_time_avg = Z_flux_tor_avg.mean()
    	Bmag=np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:]))
    	Gamma = np.array((-br[Flux_ip]*Z_flux_tor_avg[:] + bz[Flux_ip]*R_flux_tor_avg[:])/Bmag[Flux_ip])
    	Gamma_time_avg = Gamma.mean()
    	Gamma_list.append(Gamma_time_avg)
        ts_list.append(ts)
    	file.write(str(ts)+"\t"+str(Gamma_time_avg)+"\n")
    file.close()
    fig = plt.figure()
    plt.title('mean flux vs. time interval')
    plt.xlabel('time interval')
    plt.ylabel('mean flux')
    plt.grid(True)
    plt.plot(ts_list,Gamma_list,'b.-')    
    plt.show()

def total_flux():
    global VR_all, ne_all, VZ_all, br, bz, triObj,RI,ZI,tci,Sep_ip,start    
    Rmaj = 1.67
    R_list = []
    Total_Flux = [] 
    Equilibrium = []
    Turbulent = []
    Magnetic = []
    angle_list = []
    file = open("Total of Three Fluxes_3.txt",'a')
    file.write("Coarse grid:50X50. Starting point at t=240")
    file.write("pol. angle"+"\t"+"Equilibrium flux"+"\t"+"Turbulent Flux"+"\t"+"Magnetic Flux"+"\t"+"Total flux"+"\n")
    unit_R = (Rmax-Rmin)/100
    unit_Z = (Zmax-Zmin)/100
    for cut in range(10,15):
        Flux_ip_2 = getiparrs(0.95,1.,jcut)[-2]
        Z_val = (cut-jcut)*unit_Z #+ Zmin 
        R_val = Flux_ip_2*unit_R + Rmin-Rmaj
        tangent = (Z_val/R_val)
        angle = np.degrees(np.arctan(tangent))
        angle_list.append(angle)
        R_list.append(R_val)
        Eq_flux, Tur_flux, Mag_flux = three_fluxes(cut) 
        tot_flux = Eq_flux + Tur_flux + Mag_flux
        Total_Flux.append(Eq_flux + Tur_flux + Mag_flux)
        Equilibrium.append(Eq_flux)
        Turbulent.append(Tur_flux)
        Magnetic.append(Mag_flux)
        file.write(str(angle)+"\t"+str(Eq_flux)+"\t"+str(Tur_flux)+"\t"+str(Mag_flux)+"\t"+str(tot_flux)+"\n")
        print("cut done...")
    flux_tot = np.asarray(Total_Flux)
    mean_gamma = flux_tot.mean()
    file.write("averaged flux:"+ str(mean_gamma)+"\n")
    file.close()
    #fig = plt.figure()
    #plt.title('total flux vs. poloidal angle')
    #plt.xlabel('poloidal angle (degrees)')
    #plt.ylabel('total flux')
    #plt.grid(True)
    #plt.plot(angle_list,Total_Flux,'r.-')    
    #plt.show()
    #fig = plt.figure()
    #plt.title('equilibrium flux vs. poloidal angle')
    #plt.xlabel('poloidal angle (degrees)')
    #plt.ylabel('equilibrium flux')
    #plt.grid(True)
    #plt.plot(angle_list,Equilibrium,'r.-')    
    #plt.show()
    #fig = plt.figure()
    #plt.title('turbulent flux vs. poloidal angle')
    #plt.xlabel('poloidal angle (degrees)')
    #plt.ylabel('turbulent flux')
    #plt.grid(True)
    #plt.plot(angle_list,Turbulent,'r.-')    
    #plt.show()
    #fig = plt.figure()
    #plt.title('magnetic flux vs. poloidal angle')
    #plt.xlabel('poloidal angle (degrees)')
    #plt.ylabel('magnetic flux')
    #plt.grid(True)
    #plt.plot(angle_list,Magnetic,'r.-')    
    #plt.show()
    return np.asarray(Equilibrium),np.asarray(Turbulent),np.asarray(Magnetic),np.asarray(Total_Flux),np.asarray(R_list), angle_list

def den_plot():#numbers for ti255 Cut them down significantly so it can run
    Flux_ip = getiparrs(0.95,1.,jcut)[-2]
    print('Flux_ip=',Flux_ip)
    file = open("ti255_density.txt",'a')
    n_e=loader.calcNeTotal()
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],jcut)[Flux_ip] for iz in range(Nplanes)] for it in range(0,693)])
    ne_tor_avg = ne_all.mean(axis=1)
    print('n_e:',ne_tor_avg)
    for i in range(0,693):
        file.write(str(ne_tor_avg[i])+"\n")
    file.close()
    time_list = np.array([x for x in range(0,693)])
    fig = plt.figure()
    plt.title('Toroidally averaged density vs. time at Flux point')
    plt.xlabel('time step')
    plt.ylabel('density')
    plt.grid(True)
    plt.plot(time_list,ne_tor_avg,'b.-')    
    plt.show()

def dden_plot():#not complete yet. try to find the correct numbers for dn
    Flux_ip = getiparrs(0.95,1.,jcut)[-2]
    print('Flux_ip=',Flux_ip)
    file = open("ti255_den_pert.txt",'a')
    n_e=loader.calcNeTotal()
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],jcut)[Flux_ip] for iz in range(Nplanes)] for it in range(0,693)])
    ne_tor_avg = ne_all.mean(axis=1)
    dn = np.array(ne_all[:,:]-ne_tor_avg[:,np.newaxis])
   # dn_tor_avg = dn.mean(axis=1)
    for i in range(0,693):
        file.write(str(dn[i][0])+"\n")
    file.close()
   # time_list = np.array([x for x in range(0,693)])
   # fig = plt.figure()
   # plt.title('Density perturbation vs. time at Flux point')
   # plt.xlabel('time step')
   # plt.ylabel('density perturbation')
   # plt.grid(True)
   # plt.plot(time_list,dn_tor_avg,'b.-')    
   # plt.show()



def  three_fluxes(cut):
    global VR_all, ne_all, VZ_all, br, bz, triObj,RI,ZI,tci,Sep_ip,start
    #Inputs:
    time_range=85
    start=240 #start where linear phase ends
    ts_max=30
    file = open("The three fluxes.txt",'a')
    file.write("Equilibrium flux"+"\t"+"Turbulent Flux"+"\t"+"Magnetic Flux"+"\n")
    Flux_ip = getiparrs(0.95,1.,cut)[-2]
    print('Flux_ip=',Flux_ip)
    n_e=loader.calcNeTotal()
    temp = np.array([[getcutvRvZ(iz,it,cut) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    VR_all = temp[:,:,0,Flux_ip]
    VZ_all = temp[:,:,1,Flux_ip]
    #VR_all = np.array([[getcutvRvZ(iz,it,cut)[0][Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    print("VR_all and VZ_all done...")
    #VZ_all = np.array([[getcutvRvZ(iz,it,cut)[1][Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    #print("VZ_all done...")
    ne_all = np.array([[getcutvalue(n_e[:,iz,it],cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    print("ne_all done...")
    gradBR,gradBZ,_,curvR,curvZ,_ = magnetic_drifts(cut,Flux_ip,start)
    print("magnetic drifts done...")
    br = getcutvalue(bfield[:,0],cut)
    bz = getcutvalue(bfield[:,1],cut)
    ne_tor_avg = ne_all.mean(axis=1)
    VR_tor_avg = VR_all.mean(axis=1)
    VZ_tor_avg = VZ_all.mean(axis=1)
    Tur_Gamma_list=[]
    l=1
    m=0
    ne_m = []
    VR_m = []
    VZ_m = []
    for i in range(0,ts_max):#left side
	ne_m.append(np.mean(ne_tor_avg[i-m:i+ts_max+m]))
        VR_m.append(np.mean(VR_tor_avg[i-m:i+ts_max+m]))
        VZ_m.append(np.mean(VZ_tor_avg[i-m:i+ts_max+m]))
        m=m+1
    for i in range(ts_max,time_range-ts_max):#center
        ne_m.append(np.mean(ne_tor_avg[i-ts_max:i+ts_max]))
        VR_m.append(np.mean(VR_tor_avg[i-ts_max:i+ts_max]))
        VZ_m.append(np.mean(VZ_tor_avg[i-ts_max:i+ts_max]))
    for i in range(time_range-ts_max,time_range):#right side
	ne_m.append(np.mean(ne_tor_avg[i-ts_max-l:i+ts_max-l]))
        VR_m.append(np.mean(VR_tor_avg[i-ts_max-l:i+ts_max-l]))
        VZ_m.append(np.mean(VZ_tor_avg[i-ts_max-l:i+ts_max-l]))
        l=l+1
    ne_avg = np.asarray(ne_m)
    VR_avg = np.asarray(VR_m)
    VZ_avg = np.asarray(VZ_m)
    dn = np.array(ne_all[:,:] - ne_avg[:,np.newaxis])
    dVR = np.array(VR_all[:,:] - VR_avg[:,np.newaxis])
    dVZ = np.array(VZ_all[:,:] - VZ_avg[:,np.newaxis])
    Bmag=np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:]))
    #Equilibrium Flux
    Eq_R_flux = np.array(ne_tor_avg[:]*VR_tor_avg[:])
    Eq_Z_flux = np.array(ne_tor_avg[:]*VZ_tor_avg[:])        
    Eq_Gamma = np.array((-br[Flux_ip]*Eq_Z_flux[:] + bz[Flux_ip]*Eq_R_flux[:])/Bmag[Flux_ip])
    Eq_Gamma_time_avg = Eq_Gamma.mean()
    #Magnetic Flux
    Mag_R_flux = np.array(ne_tor_avg[:]*(curvR+gradBR))
    Mag_Z_flux = np.array(ne_tor_avg[:]*(curvZ+gradBZ))
    Mag_Gamma = np.array((-br[Flux_ip]*Mag_Z_flux[:] + bz[Flux_ip]*Mag_R_flux[:])/Bmag[Flux_ip])
    Mag_Gamma_time_avg = Mag_Gamma.mean()
    #Turbulent Flux
    Tur_R_flux = np.array([[dn[it,pl]*dVR[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    Tur_Z_flux = np.array([[dn[it,pl]*dVZ[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    Tur_R_flux_tor_avg = Tur_R_flux.mean(axis=1)
    Tur_Z_flux_tor_avg = Tur_Z_flux.mean(axis=1)
    Tur_R_flux_time_avg = Tur_R_flux_tor_avg.mean()
    Tur_Z_flux_time_avg = Tur_Z_flux_tor_avg.mean()
    Tur_Gamma = np.array((-br[Flux_ip]*Tur_Z_flux_tor_avg[:] + bz[Flux_ip]*Tur_R_flux_tor_avg[:])/Bmag[Flux_ip])
    Tur_Gamma_time_avg = Tur_Gamma.mean()
    Tur_Gamma_list.append(Tur_Gamma_time_avg)
    file.write(str(Eq_Gamma_time_avg)+"\t"+str(Tur_Gamma_time_avg)+"\t"+str(Mag_Gamma_time_avg)+"\n")
    file.close()
    return Eq_Gamma_time_avg, Tur_Gamma_time_avg, Mag_Gamma_time_avg
    #print("Equilibrium Flux=,",Eq_Gamma_time_avg)
    #print("Turbulent Flux=,",Tur_Gamma_time_avg)
    #print("Magnetic Flux=,",Mag_Gamma_time_avg)
    #fig = plt.figure()
    #plt.title('mean flux vs. time interval')
    #plt.xlabel('time interval')
    #plt.ylabel('mean flux')
    #plt.grid(True)
    #plt.plot(ts_list,Gamma_list,'b.-')    
    #plt.show()

def Line_int_along_curve(rarray,thetalist,weight):#function that returns the line integral of a weight along a curve
    #index of starting and ending angles. Assume that all three arrays have the same size
    theta_start = 0
    theta_end = len(thetalist)-1
    prefactor = (thetalist[theta_end]-thetalist[theta_start])/(len(thetalist[theta_start:theta_end])-1)
    #dr/dtheta
    drdtheta = []
    for i in range(1,theta_end):
        drdtheta.append((rarray[i+1]-rarray[i-1])/(thetalist[i+1]-thetalist[i-1]))
    #linear correction for first and last element of derivative list
    dr2dtheta_start = (rarray[2]-2*rarray[1] + rarray[0])/(thetalist[2]-thetalist[0])**2
    dr2dtheta_end = (rarray[theta_end]-2*rarray[theta_end-1] + rarray[theta_end-2])/(thetalist[theta_end]-thetalist[theta_end-2])**2
    last_elem = drdtheta[-1] + dr2dtheta_end
    drdtheta.append(last_elem)
    first_elem = drdtheta[0] - dr2dtheta_start
    a = first_elem
    drdtheta_new = [a]+drdtheta
    der = np.asarray(drdtheta_new)
    #Integration scheme of the full square root 
    temp = weight*rarray*np.sqrt(1 + (1/np.power(rarray,2))*np.power(der,2))
    Summation = (1/2)*weight[theta_start]*rarray[theta_start]*np.sqrt(1+(1/rarray[theta_start]**2)*der[theta_start]**2) +\
    (1/2)*weight[theta_end]*rarray[theta_end]*np.sqrt(1+(1/rarray[theta_end]**2)*der[theta_end]**2) +\
    np.sum(temp[theta_start+1:theta_end-2])
    Integral = prefactor*Summation
    print("Integral = ",Integral)
    return Integral



















def magnetic_drifts(cut,ip_value,start):
    global Zi,Ri,pot,bfield,B2
    if cut==None:
        cut = jcut
    dZ = Zi[1]-Zi[0]
    dR = Ri[1]-Ri[0]
    theta = (2*math.pi)/Nplanes

    #Define flux points, constants, average temperatures and velocities on them
    m_i = 2*1.6723e-27
    k = 1.38e-23
    joule_conv = 1.602e-19
    charge = 1.602176e-19 
    permeability = 1.2566e-6
    #c = 3e8
    Flux_ip = getiparrs(0.95,1.,cut)[-2]
    psi_val = getcutvalue(psin[:],cut)[Flux_ip]
    oned_location = (np.where(psin1d>psi_val))[0][0]
    Tigc1d = loader.Ti1d
    Tigcperp = loader.Tiperp
    Tigcpar = loader.Tipar
    Ti_prp_avg = np.mean(Tigcperp[start:,oned_location]) 
    Ti_prl_avg = 2*np.mean(Tigcpar[start:,oned_location]) 
    Ti_avg = np.mean(Tigc1d[start:,oned_location])
    V_th = np.sqrt((joule_conv*Ti_avg)/m_i)
   # print('V_th:', '%.2E' % Decimal(V_th))
   # print('Temp_ratio:','%.2E' % Decimal(Ti_prp_avg/Ti_prl_avg))
    V_perp_2 = (2*joule_conv*Ti_prp_avg)/m_i
    V_par_2 = (joule_conv*Ti_prl_avg)/m_i
   # print('V_perp_sq=','%.2E' % Decimal(V_perp_2))
   # print('V_par_sq=','%.2E' % Decimal(V_par_2))
    #Loading B-field 
    #bfield=loader.loadBfield()
    #calculate values at the cut and find the gyrofrequency
    Br = getcutvalue(bfield[:,0],cut)
    Bz = getcutvalue(bfield[:,1],cut)
    Bzeta = getcutvalue(bfield[:,2],cut)
    BfMag=[]
    for i in range(len(Br)):
        BfMag.append(np.sqrt(Br[i]*Br[i] + Bz[i]*Bz[i] + Bzeta[i]*Bzeta[i]))
    BfMag = np.asarray(BfMag)
    #print('B:','%.2E' % Decimal(BfMag[Flux_ip]))
    Omega = (charge*BfMag/(m_i))
    #print('Omega:',np.ma.masked_invalid(Omega).mean(),np.std(np.ma.masked_invalid(Omega)))
   # print('Omega:','%.2E' % Decimal(Omega[Flux_ip]))
   # print('rho:','%.2E' % Decimal(V_th/Omega[Flux_ip]))

    mu = (0.5*m_i*V_perp_2)/BfMag
    #print('B=',BfMag)
    #print('mu=',mu)
    #print('Omega=',Omega)
    #finding J_perp/B
   # i_pol_flow = loader.i_pol_f
   # e_pol_flow = loader.e_pol_f
   # i_pol_flow_av = i_pol_flow.mean(axis=0)
   # e_pol_flow_av = e_pol_flow.mean(axis=0)
   # n_e=loader.calcNeTotal()
   # ne_all = np.array([[getcutvalue(n_e[:,iz,it],cut)[Flux_ip] for iz in range(Nplanes)] for it in range(0,130)])
   # ne_int = ne_all.mean(axis=0)
   # ne = ne_int.mean(axis=0)
    #print('ne',ne_all.shape)
   # J_perp = ne*charge*(i_pol_flow_av[oned_location]-e_pol_flow_av[oned_location])*permeability
   # print('i_poloidal:','%.2E' % Decimal(i_pol_flow_av[oned_location]))
   # print('e_poloidal:','%.2E' % Decimal(e_pol_flow_av[oned_location]))
   # print('J_perp:','%.2E' % Decimal(J_perp))
   # print('J_perp/B:','%.2E' % Decimal(J_perp/BfMag[Flux_ip]))


    #Calculating the gradient of scalar |B| and dbr/dZ, dbr/dR, dbz/dZ, dbz/dR, dbzeta/dZ, dbzeta/dR
# dB/dZ
    Brm1 = getcutvalue(bfield[:,0],cut-1)
    Bzm1 = getcutvalue(bfield[:,1],cut-1)
    Bzetam1 = getcutvalue(bfield[:,2],cut-1)
    BfMagm1=[]
    for i in range(len(Brm1)):
        BfMagm1.append(np.sqrt(Brm1[i]*Brm1[i] + Bzm1[i]*Bzm1[i] + Bzetam1[i]*Bzetam1[i]))
    BMm1 = np.asarray(BfMagm1)
    brm1 = Brm1/BMm1
    bzm1 = Bzm1/BMm1
    bzetam1 = Bzetam1/BMm1
    Br = getcutvalue(bfield[:,0],cut)
    Bz = getcutvalue(bfield[:,1],cut)
    Bzeta = getcutvalue(bfield[:,2],cut)
    BfMag=[]
    for i in range(len(Br)):
        BfMag.append(np.sqrt(Br[i]*Br[i] + Bz[i]*Bz[i] + Bzeta[i]*Bzeta[i]))
    BM = np.asarray(BfMag)
    br = Br/BM
    bz = Bz/BM
    bzeta = Bzeta/BM
    Brp1 = getcutvalue(bfield[:,0],cut+1)
    Bzp1 = getcutvalue(bfield[:,1],cut+1)
    Bzetap1 = getcutvalue(bfield[:,2],cut+1)
    BfMagp1=[]
    for i in range(len(Brp1)):
        BfMagp1.append(np.sqrt(Brp1[i]*Brp1[i] + Bzp1[i]*Bzp1[i] + Bzetap1[i]*Bzetap1[i]))
    BMp1 = np.asarray(BfMagp1)
    brp1 = Brp1/BMp1
    bzp1 = Bzp1/BMp1
    bzetap1 = Bzetap1/BMp1
    dBdZ = (BMp1-BMm1)/(2.*dZ)
    dbrdZ = (brp1-brm1)/(2.*dZ)
    dbzdZ = (bzp1-bzm1)/(2.*dZ)
    dbzetadZ = (bzetap1-bzetam1)/(2.*dZ)
    #print('dB/dZ=',dBdZ) #print('dbr/dZ=',dbrdZ)
    #print('dbz/dZ=',dbzdZ)
    #print('dbzeta/dZ=',dbzetadZ)
# dB/dR
    dBdR=[]
    dbrdR=[]
    dbzdR=[]
    dbzetadR=[]
    lengthR = len(getcutvalue(bfield[:,0],cut))
    iterlistR = list(range(0,lengthR))
    Brrm10 = getcutvalue(bfield[:,0],cut)[lengthR-1]
    Bzrm10 = getcutvalue(bfield[:,1],cut)[lengthR-1]
    Bzetarm10 = getcutvalue(bfield[:,2],cut)[lengthR-1]
    BfMagrm10=np.sqrt(Brrm10*Brrm10 + Bzrm10*Bzrm10 + Bzetarm10*Bzetarm10)
    brrm10 = Brrm10/BfMagrm10
    bzrm10 = Bzrm10/BfMagrm10
    bzetarm10 = Bzetarm10/BfMagrm10
    Brrp10 = getcutvalue(bfield[:,0],cut)[1]
    Bzrp10 = getcutvalue(bfield[:,1],cut)[1]
    Bzetarp10 = getcutvalue(bfield[:,2],cut)[1]
    BfMagrp10=np.sqrt(Brrp10*Brrp10 + Bzrp10*Bzrp10 + Bzetarp10*Bzetarp10)
    brrp10 = Brrp10/BfMagrp10
    bzrp10 = Bzrp10/BfMagrp10
    bzetarp10 = Bzetarp10/BfMagrp10
    dbrdR.append((brrp10-brrm10)/(2.*dR))
    dbzdR.append((bzrp10-bzrm10)/(2.*dR))
    dbzetadR.append((bzetarp10-bzetarm10)/(2.*dR))
    dBdR.append((BfMagrp10-BfMagrm10)/(2.*dR))
    for i in iterlistR[1:-1]:
        Brrm1 = getcutvalue(bfield[:,0],cut)[i-1]
        Bzrm1 = getcutvalue(bfield[:,1],cut)[i-1]
        Bzetarm1 = getcutvalue(bfield[:,2],cut)[i-1]
        BMrm1 = np.sqrt(Brrm1*Brrm1 + Bzrm1*Bzrm1 + Bzetarm1*Bzetarm1)
        brrm1 = Brrm1/BMrm1
        bzrm1 = Bzrm1/BMrm1
        bzetarm1 = Bzetarm1/BMrm1
        Brrp1 = getcutvalue(bfield[:,0],cut)[i+1]
        Bzrp1 = getcutvalue(bfield[:,1],cut)[i+1]
        Bzetarp1 = getcutvalue(bfield[:,2],cut)[i+1]
        BMrp1=np.sqrt(Brrp1*Brrp1 + Bzrp1*Bzrp1 + Bzetarp1*Bzetarp1)
        brrp1 = Brrp1/BMrp1
        bzrp1 = Bzrp1/BMrp1
        bzetarp1 = Bzetarp1/BMrp1
        dbrdR.append((brrp1-brrm1)/(2.*dR))
        dbzdR.append((bzrp1-bzrm1)/(2.*dR))
        dbzetadR.append((bzetarp1-bzetarm1)/(2.*dR))
        dBdR.append((BMrp1-BMrm1)/(2.*dR))
    Brrm1L = getcutvalue(bfield[:,0],cut)[lengthR-2]
    Bzrm1L = getcutvalue(bfield[:,1],cut)[lengthR-2]
    Bzetarm1L = getcutvalue(bfield[:,2],cut)[lengthR-2]
    BfMagrm1L=np.sqrt(Brrm1L*Brrm1L + Bzrm1L*Bzrm1L + Bzetarm1L*Bzetarm1L)
    brrm1L = Brrm1L/BfMagrm1L
    bzrm1L = Bzrm1L/BfMagrm1L
    bzetarm1L = Bzetarm1L/BfMagrm1L
    Brrp1L = getcutvalue(bfield[:,0],cut)[0]
    Bzrp1L = getcutvalue(bfield[:,1],cut)[0]
    Bzetarp1L = getcutvalue(bfield[:,2],cut)[0]
    BfMagrp1L=np.sqrt(Brrp1L*Brrp1L + Bzrp1L*Bzrp1L + Bzetarp1L*Bzetarp1L)
    brrp1L = Brrp1L/BfMagrp1L
    bzrp1L = Bzrp1L/BfMagrp1L
    bzetarp1L = Bzetarp1L/BfMagrp1L
    dbrdR.append((brrp1L-brrm1L)/(2.*dR))
    dbzdR.append((bzrp1L-bzrm1L)/(2.*dR))
    dbzetadR.append((bzetarp1L-bzetarm1L)/(2.*dR))
    dBdR.append((BfMagrp1L-BfMagrm1L)/(2.*dR))
    dbrdR = np.asarray(dbrdR)
    dbzdR = np.asarray(dbzdR)   
    dbzetadR = np.asarray(dbzetadR)
    dBdR = np.asarray(dBdR)
    #print('dbr/dR=',dbrdR)
    #print('dbz/dR=',dbzdR)
    #print('dbzeta/dR=',dbzetadR)
    #dBdzeta
    dBdzeta = 0
    dbrdzeta = 0
    dbzdzeta = 0
    dbzetadzeta = 0
    L = np.array([x for x in range(0,100)])
    unit = (Rmax-Rmin)/100
    R = np.array([L[ip]*unit+Rmin for ip in range(0,len(br))])

    
    #calculation of gradBr, gradBz, gradBzeta drifts
    gradBR = (bz*(dBdzeta/R) - bzeta*dBdZ)*(mu/(m_i*Omega))
    gradBZ = (dBdR*bzeta - br*(dBdzeta/R))*(mu/(m_i*Omega))
    gradBzeta = (dBdR*bz - dBdZ*br)*(mu/(m_i*Omega))
    gradBZ_geom = (dBdR*bzeta - br*(dBdzeta/R))
   # print('gradBR:',np.ma.masked_invalid(gradBR).mean(),np.std(np.ma.masked_invalid(gradBR)))
   # print('gradBZ:',np.ma.masked_invalid(gradBZ).mean(),np.std(np.ma.masked_invalid(gradBZ)))
   # print('gradBzeta:',np.ma.masked_invalid(gradBzeta).mean(),np.std(np.ma.masked_invalid(gradBzeta)))
   # print('gradBR:','%.2E' % Decimal(gradBR[Flux_ip]))
   # print('gradBZ:','%.2E' % Decimal(gradBZ[Flux_ip]))
   # print('gradBzeta:','%.2E' % Decimal(gradBzeta[Flux_ip]))
   # print('gradBZ_geom:','%.2E' % Decimal(gradBZ_geom[Flux_ip]))
    
    #calculation of curvature drifts
    #b*delb
    #b*Del = br*d/dr+bz*d/dz+bzeta*d/dzeta
    #(b*Del)(br,bz,bzeta) = 
    #R= br*dbrdR+bz*dbrdZ+bzeta*dbrdzeta
    #Z =  br*dbzdR+bz*dbzdZ+bzeta*dbzdzeta
    #zeta = br*dbzetadR+bz*dbzetadZ+bzeta*dbzetadzeta
    #bX(b*Del)b=
    curvR = (bz*(br*dbzetadR+bz*dbzetadZ+(bzeta/R)*dbzetadzeta+(bzeta*br)/R)-bzeta*(br*dbzdR+bz*dbzdZ+(bzeta/R)*dbzdzeta))*(V_par_2/Omega)
    curvZ = (bzeta*(br*dbrdR+bz*dbrdZ+(bzeta/R)*dbrdzeta - (bzeta*bzeta)/R)-br*(br*dbzetadR+bz*dbzetadZ+(bzeta/R)*dbzetadzeta+(bzeta*br)/R))*(V_par_2/Omega)
    curvzeta = (br*(br*dbzdR+bz*dbzdZ+(bzeta/R)*dbzdzeta)-bz*(br*dbrdR+bz*dbrdZ+(bzeta/R)*dbrdzeta - (bzeta*bzeta)/R))*(V_par_2/Omega)
    curvZ_geom = (bzeta*(br*dbrdR+bz*dbrdZ+(bzeta/R)*dbrdzeta - (bzeta*bzeta)/R)-br*(br*dbzetadR+bz*dbzetadZ+(bzeta/R)*dbzetadzeta+(bzeta*br)/R))

   # print('curvR:',np.ma.masked_invalid(curvR).mean(),np.std(np.ma.masked_invalid(curvR)))
   # print('curvZ:',np.ma.masked_invalid(curvZ).mean(),np.std(np.ma.masked_invalid(curvZ)))
   # print('curvzeta:',np.ma.masked_invalid(curvzeta).mean(),np.std(np.ma.masked_invalid(curvzeta)))
    return gradBR[ip_value],gradBZ[ip_value],gradBzeta[ip_value],curvR[ip_value],curvZ[ip_value],curvzeta[ip_value]
   # print('curvR:','%.2E' % Decimal(curvR[Flux_ip]))
   # print('curvZ:','%.2E' % Decimal(curvZ[Flux_ip]))
   # print('curvzeta:','%.2E' % Decimal(curvzeta[Flux_ip]))
   # print('curvZ_geom','%.2E' % Decimal(curvZ_geom[Flux_ip]))
   # print('ratio of geometric factors:','%.2E' % Decimal(gradBZ_geom[Flux_ip]/curvZ_geom[Flux_ip]))
   # print('ratio of drifts(gradBZ/curvZ):','%.2E' % Decimal(gradBZ[Flux_ip]/curvZ[Flux_ip]))

def test_function():
    #subprocess.call("load_data.py")
    Eq, Tur, Mag, Tot, R, angle = total_flux()
    Eq_f = Line_int_along_curve(R,angle,Eq)
    Tur_f = Line_int_along_curve(R,angle,Tur)
    Mag_f = Line_int_along_curve(R,angle,Mag)
    Tot_f = Line_int_along_curve(R,angle,Tot)
    Tot2_f = Line_int_along_curve(R,angle,Eq+Tur+Mag)
    print("Equilibrium Flux Integral:", Eq_f)
    print("Turbulent Flux Integral:", Tur_f)
    print("Magnetic Flux Integral:", Mag_f)
    print("Total Flux Integral:", Tot_f)
    print("Total2 Flux Integral:", Tot2_f)
