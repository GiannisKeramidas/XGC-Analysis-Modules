import numpy as np
import xgc
import math
import pylab as py
#import xgcjrm as xgc
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
import sys
import cProfile
import glob
import angle
import numpy.ma as ma
from decimal import Decimal


#limit data to the [Rmin,Rmax,Zmin,Zmax] box, and read only the first two toroidal planes
#Rmin=None
Rmin=0.4401#0.88#2.15#1.8#2.15#1.9#2.2
#Rmax=None
Rmax=0.9099#0.9#2.26#2.19#2.26#2.2#2.31
#Zmin=None
Zmin=-0.5958#-0.1#-0.28#-0.8#-0.28#-0.9#-0.25
#Zmax=None
Zmax=0.4328#0.1#0.55#-0.28#0.55#-0.3#0.4
phi_start=0
phi_end=None
#Rmaj=1.67
Rmaj = 0.68

#fileDirec='/scratch1/scratchdirs/giannos/particle_pinch/particle_pinch_idl/hdf5'
#fileDir= 'D:\particle_pinch_idl\hdf5'
#fileDir= '/global/homes/g/giannos/ti253_d3d_Ip_1.5_large'
#fileDir='/scratch1/scratchdirs/giannos/ti255_d3d_Ip_1.5_med'
fileDir='/global/cscratch1/sd/giannos/ti255_d3d_Ip_1.5_med'
#fileDir='/global/cscratch1/sd/giannos/ti262_cmod_JTM'

def getMeshAndCuts(fileDir,Rmin,Rmax,Zmin,Zmax,nRi=None,nZi=None):

    global loader,ne,pot,psin,RZ,tri,time,psin1d,psin001d,Te1d,ne1d,pot001d,bfield,pot0,dpot,dpot,pot,Tigc1d,nigc1d,i_T_perp,i_E_para,e_T_perp,e_E_para,i_u_para, e_u_para
    global Nplanes,Ntimes,dt
    global Ri,Zi,RI,ZI,triObj,unitR,unitZ
    global sepInds,psinvalues,jcut,Zcut,Rcut,rcut,fsInds#,psii

    
    loader=xgc.load(fileDir,Rmin=Rmin,Rmax=Rmax,Zmin=Zmin,Zmax=Zmax,phi_start=0,phi_end=None)
    ne = loader.calcNeTotal()
    neAd = loader.calcNeAdiabatic()
    pot = loader.calcPotential()
    print('Rmin = ',loader.Rmin)
    print('Rmax = ',loader.Rmax)
    print('Zmin = ',loader.Zmin)
    print('Zmax = ',loader.Zmax)
    print('ne.shape = ',ne.shape)

    # (could have imported * from xgc but this way we see explicitly what we are getting)

    # arrays
    #rdeflux = loader.e_radial_pflux
    psin = loader.psin
    RZ = loader.RZ
    tri = loader.tri
    time = loader.time
    psin1d = loader.psin1d
    psin001d = loader.psin001d
    Te1d = loader.Te1d
    ne1d = loader.ne1d
    pot001d = loader.pot001d
    bfield = loader.bfield
    pot0 = loader.pot0
    dpot = loader.dpot
    pot = loader.pot
    eden = loader.eden#non-adiabatic part of the density
    #e_rad_fl = loader.e_rad_fl
    #e_rad_exb_fl = loader.e_rad_exb_fl
    #if hasattr(loader,'i_T_perp'):
    i_u_para = loader.i_u_para
    i_T_perp = loader.i_T_perp
    i_E_para = loader.i_E_para
    e_T_perp = loader.e_T_perp
    e_E_para = loader.e_E_para
    e_u_para = loader.e_u_para
    
    # scalars
    Zcut = None
    Rcut = None
    Nplanes = loader.Nplanes
    Ntimes = loader.Ntimes
    dt = loader.dt

    #setup mesh grid
    if nRi == None:
        nRi = 100
    if nZi == None:
        nZi = 100
            
    Ri = np.linspace(loader.Rmin,loader.Rmax,nRi)    
    Zi = np.linspace(loader.Zmin,loader.Zmax,nZi)
    (RI,ZI) = np.meshgrid(Ri,Zi)
    # set up to interpolate using the TriInterpolator class. Should be what tricontourf() uses
    triObj = Triangulation(RZ[:,0],RZ[:,1],tri)

    # get indices of RZ vertices on the separatrix
    sepInds = np.where(np.abs(loader.psin-1.0)<1e-4)[0]

    # flux surface values on the triangular mesh
    psinvalues = np.sort(np.array(list(set(np.round(psin,4))))) # 4 round to 4 decimal places; adjust?
    psinvalues = list(psinvalues)
    #ind_of_sep = psinvalues.index(1.)#find the index of separatrix flux
    #inner_fs_ind = 177#index of flux of flux surface before the separatrix
    
    # get indices of RZ vertices on the inner flux surface
    #fs_val = 0.973#psinvalues[inner_fs_ind]
    #fsInds = np.where(np.abs(loader.psin-fs_val)<1e-2)[0]


    # locate j value where Z = Zcut
    if Zcut == None:
        Zcut = Zi.mean()
        jcut = int(round(np.mean(np.where(np.abs(Zi-Zcut)<1.e-1))))
        #print "Zi[0:5] = ",Zi[0:5]
        #print "Zi[-5:0] = ",Zi[-5:-1]
        #print "Zcut = ",Zcut
        print("jcut = ",jcut)

    if Rcut == None:
        Rcut = Ri.mean()
        rcut = int(round(np.mean(np.where(np.abs(Ri-Rcut)<1.e-1))))
        #print "Zi[0:5] = ",Zi[0:5]
        #print "Zi[-5:0] = ",Zi[-5:-1]
        #print "Zcut = ",Zcut
        print("rcut = ",rcut)


        # calculate psi(Ri) = psii at the cut
        tci=LinearTriInterpolator(triObj,psin)
        psiout=tci(RI,ZI)
        psii = psiout[jcut,:] # a 1D array of psin values at each Ri along the cut

        print(loader.Rmin)
        #unitR = (float(loader.Rmax)-float(loader.Rmin))/len(Ri)
        #unitZ = (float(loader.Zmax)-float(loader.Zmin))/len(Zi)
        
        
def getcutvalue_hor(arrRZ,j,option):
    #"""Calculates the values at Ri along Zcut of arrRZ where arr is a mesh array
    #usage:
     #   netmpi = getcutvalue_hor(ne[:,iz,it]) # where iz is the toroidal plane, and it the time index
    #"""
    #option=1 is for evaluating on the rectangular grid
    #option=2 is for evaluating on the grid by spline-generated flux points, equidistant in poloidal angle.
    global triObj,RI,ZI,arrout,tci,jcut
    if option == 1:
    # the following construction allows the def to occur before jcut takes a value
        if j==None:
            j = jcut 
        tci=LinearTriInterpolator(triObj,arrRZ)
        arrout=tci(RI,ZI)
        return(arrout[j,:])
    if option == 2:
        if j==None:
            j = int(round(np.mean(angle.ZA))) 
        tci=LinearTriInterpolator(triObj,arrRZ)
        arrout=tci(angle.RA,angle.ZA)
        return(arrout[j,:])
    
def getcutvalue_ver(arrRZ,r,option):
    """Calculates the values at Zi along Rcut of arrRZ where arr is a mesh array
    usage:
        netmpi = getcutvalue_ver(ne[:,iz,it]) # where iz is the toroidal plane, and it the time index
    """
    global triObj,RI,ZI,arrout,tci,jcut,RI,ZI
    # the following construction allows the def to occur before jcut takes a value
    if option==1:
        if r==None:
            r = rcut 
        tci=LinearTriInterpolator(triObj,arrRZ)
        arrout=tci(RI,ZI)
        return(arrout[:,r])
    if option == 2:
        if r==None:
            r = int(round(np.mean(angle.ZA))) 
        tci=LinearTriInterpolator(triObj,arrRZ)
        arrout=tci(angle.RA,angle.ZA)
        return(arrout[:,r])
    
def getcutvRvZ_old(iz,it,cut):
    global Zi,Ri,pot,bfield,B2
    option=1
    if cut==None:
        cut = jcut
    dZ = Zi[1]-Zi[0]
    dR = Ri[1]-Ri[0]
    theta = (2*math.pi)/Nplanes
    #Loading B-field 
    #bfield=loader.loadBfield()
    #calculate values at the cut
    br = getcutvalue_hor(bfield[:,0],cut,option)
    bz = getcutvalue_hor(bfield[:,1],cut,option)
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    BfMag=[]
    for i in range(len(br)):
        BfMag.append(np.sqrt(br[i]*br[i] + bz[i]*bz[i] + bzeta[i]*bzeta[i]))
    #print(BfMag)
#Calculating the potential gradient
# dphi/dZ    
    potim1 = getcutvalue_hor(pot[:,iz,it],cut-1,option)
    poti = getcutvalue_hor(pot[:,iz,it],cut,option)
    potip1 = getcutvalue_hor(pot[:,iz,it],cut+1,option)
    dphidZ = (potip1-potim1)/(2.*dZ)
# dphi/dR
    dphidR=[]
    lengthR = len(getcutvalue_hor(pot[:,iz,it],cut,option))
    iterlistR = list(range(0,lengthR))
    potrm10 = getcutvalue_hor(pot[:,iz,it],cut,option)[lengthR-1]
    potrp10 = getcutvalue_hor(pot[:,iz,it],cut,option)[1]
    dphidR.append((potrp10 - potrm10)/(2.*dR))
    for i in iterlistR[1:-1]:
        potrm1 = getcutvalue_hor(pot[:,iz,it],cut,option)[i-1]
        potrp1 = getcutvalue_hor(pot[:,iz,it],cut,option)[i+1]
        dphidR.append((potrp1 - potrm1)/(2.*dR))
    potrm1L = getcutvalue_hor(pot[:,iz,it],cut,option)[lengthR-2]
    potrp1L = getcutvalue_hor(pot[:,iz,it],cut,option)[0]
    dphidR.append((potrp1L-potrm1L)/2*dR)
#dphi/dzeta
    if iz < (Nplanes-1):
        potzp1 = getcutvalue_hor(pot[:,iz+1,it],cut,option) 
        potz = getcutvalue_hor(pot[:,iz,it],cut,option)
        dphidzeta=[]
        lengthZeta = len(Ri)
        iterlistZeta = list(range(0,lengthZeta))
        for i in iterlistZeta:
            dphidzeta.append((potzp1[i] - potz[i])/(Ri[i]*theta))
    else:
        potzp1 = getcutvalue_hor(pot[:,0,it],cut,option)
        potz = getcutvalue_hor(pot[:,iz,it],cut,option)
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
        vR.append((dphidZ[i]*bzeta[i] - dphidzeta[i]*bz[i])/BfMag[i]**2)
    vZ = []
    for i in iterlistR:
        vZ.append((br[i]*dphidzeta[i] - dphidR[i]*bzeta[i])/BfMag[i]**2)
    return(vR,vZ,dphidR)
    
    
def getcutvRvZ(iz,it,cut):
    global Zi,Ri,pot,jcut,bfield,B2
    option=1
    if cut==None:
        cut = jcut
        
    dZ = Zi[1] - Zi[0]
    dR = Ri[1] - Ri[0]
    
    potim1 = getcutvalue_hor(pot[:,iz,it],(cut-1)%len(Zi),option)
    poti = getcutvalue_hor(pot[:,iz,it],cut,option)
    potip1 = getcutvalue_hor(pot[:,iz,it],(cut+1)%len(Zi),option)
    dphidZ = (potip1-potim1)/(2.*dZ)
    
    potiplus = np.roll(poti,-1)
    potiminus = np.roll(poti,1)
    dphidR = (potiplus-potiminus)/(2.*dR)
    dphidR[0] = (poti[1]-poti[0])/dR
    dphidR[-1] = (poti[-1]-poti[-2])/dR
    
    BRi = getcutvalue_hor(bfield[:,0],cut,option)
    BZi = getcutvalue_hor(bfield[:,1],cut,option)
    Bzetai = getcutvalue_hor(bfield[:,2],cut,option)
    
    Bzobzet2 = (BZi/Bzetai)**2
    Brobzet2 = (BRi/Bzetai)**2
    BrBzobzet2 = (BRi*BZi)/(Bzetai**2)
    B2 = BRi**2 + BZi**2 + Bzetai**2
    #print(B2)
    
    vRi = (Bzetai/B2)*(dphidZ*(1.+Bzobzet2)+BrBzobzet2*dphidR)
    vZi = -(Bzetai/B2)*(dphidR*(1.+Brobzet2)+BrBzobzet2*dphidZ)
    return(vRi,vZi)

#calculates EXB velocity matrix
def EXB_matr():
    global temp
    getMeshAndCuts(fileDir,Rmin,Rmax,Zmin,Zmax)
    temp = np.array([[[getcutvRvZ(iz,it,cut) for cut in range(0,100)] for iz in range(0,Nplanes)] for it in range(750,903)])
    np.save('temp',temp)





def av_flows(cut,it):
    option = 1
    vr,vz = zip(*[getcutvRvZ(iz,it,cut) for iz in range(0,Nplanes)])
    BRi = getcutvalue_hor(bfield[:,0],cut,option)
    BZi = getcutvalue_hor(bfield[:,1],cut,option)
    Bzetai = getcutvalue_hor(bfield[:,2],cut,option)
    Bp = np.sqrt(np.square(BRi)+np.square(BZi))
    B = np.sqrt(np.square(BRi)+np.square(BZi)+np.square(Bzetai))
    v_theta = (1/Bp)*(BRi*vr + BZi*vz)
    return v_theta

def flow_in_time(cut):
    t_start = 650
    t_end = 660
    pol_flow = [av_flows(cut,it) for it in range(t_start,t_end)]
    poloidal = np.asarray(pol_flow)
    print(poloidal.shape)
    return poloidal
    
    
    
def two_d_pol_flow():
    pol_flow = [flow_in_time(cut) for cut in range(0,len(Zi))]
    poloidal_flow = np.asarray(pol_flow)
    return poloidal_flow# flow[Z,time,planes,R]

def tor_avg(it,cut):
    global Zi,Ri,pot,jcut,bfield,B2
    option=1
    arrim1 = [getcutvalue_hor(pot[:,iz,it],cut-1,option) for iz in range(0,Nplanes)]
    return np.asarray(arrim1)
    #arri = getcutvalue_hor(arr[:,iz,it],cut,option)
    #arrip1 = getcutvalue_hor(arr[:,iz,it],cut+1,option)

#FUNCTIONS FOR THE CALCULATION AND IMAGING OF SHEAR
def psi_der(arr,it,cut):#takes the d/dpsi on the potential and divides by RBp
    global Zi,Ri,pot,jcut,bfield,B2,unitR
    option=1
    if cut==None:
        cut = jcut
        
    dZ = Zi[1] - Zi[0]
    dR = Ri[1] - Ri[0]
    
    arrim1 = np.asarray([getcutvalue_hor(arr[:,iz,it],(cut-1)%len(Zi),option) for iz in range(0,Nplanes)]).mean(axis=0)
    arri = np.asarray([getcutvalue_hor(arr[:,iz,it],cut,option) for iz in range(0,Nplanes)]).mean(axis=0)
    arrip1 = np.asarray([getcutvalue_hor(arr[:,iz,it],(cut+1)%len(Zi),option) for iz in range(0,Nplanes)]).mean(axis=0)
    darrdZ = (arrip1-arrim1)/(2.*dZ)
    
    arriplus = np.roll(arri,-1)
    arriminus = np.roll(arri,1)
    darrdR = (arriplus-arriminus)/(2.*dR)
    darrdR[0] = (arri[1]-arri[0])/dR
    darrdR[-1] = (arri[-1]-arri[-2])/dR
    
    BRi = getcutvalue_hor(bfield[:,0],cut,option)
    BZi = getcutvalue_hor(bfield[:,1],cut,option)
    Bzetai = getcutvalue_hor(bfield[:,2],cut,option)
    Bp = np.sqrt(np.square(BRi)+np.square(BZi))
    L = np.array([x for x in range(0,len(Ri))])
    R = np.array([L[ip]*unitR+Rmin for ip in range(0,len(Ri))])+Rmaj

    
    darrdpsi = (1/Bp)*(BRi*darrdZ - BZi*darrdR)
    return -darrdpsi/(R*Bp)

def two_d_psi_der(it):#produces 2D E_psi/RBp matrix
    global Zi,Ri,pot,jcut,bfield,B2
    E_rad = [psi_der(pot,it,cut) for cut in range(0,len(Zi))]
    E_psi = np.asarray(E_rad)
    return E_psi
    
def shear(it,cut):#takes the d/dpsi of E_psi/RBp and multiplies with (RBp)^2/B to find the shear.
    global Zi,Ri,pot,jcut,bfield,B2
    option = 1
    c = 3.0e8
    dZ = Zi[1] - Zi[0]
    dR = Ri[1] - Ri[0]
    
    dphidpsi = two_d_psi_der(it)
    
    arrim1 = dphidpsi[(cut-1)%len(Zi),:]
    arri = dphidpsi[cut,:]
    arrip1 = dphidpsi[(cut+1)%len(Zi),:]
    darrdZ = (arrip1-arrim1)/(2.*dZ)
    
    arriplus = np.roll(arri,-1)
    arriminus = np.roll(arri,1)
    darrdR = (arriplus-arriminus)/(2.*dR)
    darrdR[0] = (arri[1]-arri[0])/dR
    darrdR[-1] = (arri[-1]-arri[-2])/dR
    
    BRi = getcutvalue_hor(bfield[:,0],cut,option)
    BZi = getcutvalue_hor(bfield[:,1],cut,option)
    Bzetai = getcutvalue_hor(bfield[:,2],cut,option)
    Bp = np.sqrt(np.square(BRi)+np.square(BZi))
    Bmag = np.sqrt(np.square(BRi)+np.square(BZi)+np.square(Bzetai))
    L = np.array([x for x in range(0,len(Ri))])
    R = np.array([L[ip]*unitR+Rmin for ip in range(0,len(Ri))])+Rmaj

    doubledpsi = (1/Bp)*(BRi*darrdZ - BZi*darrdR)
    Shear = ((np.square(R)*np.square(Bp))/Bmag)*doubledpsi
    return Shear
    
def two_d_shear(it): #produces a 2D matrix of Shear
    global Zi,Ri,pot,jcut,bfield,B2
    Shear = [shear(it,cut) for cut in range(0,len(Zi))]
    Bur_shear = np.asarray(Shear)
    return Bur_shear



def shear_call_n_save(r):
    sh = two_d_shear(650)
    sh2 = np.nan_to_num(sh)
    np.save("shear_%s" %(r),sh2)

    
def sep_save(r):
    global RZ
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    np.save("Rsep_%s" %(r),Rpol)
    np.save("Zsep_%s" %(r),Zpol)
    
def two_d_plot(arr):
    global RZ
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    fig,ax = plt.subplots()
    ax.plot(Rpol,Zpol,'k--')
    dmove = 1
    Zs, Rs = np.mgrid[slice(0,arr.shape[0],dmove),slice(0,arr.shape[1],dmove)]
    z = arr[Zs,Rs]
    im = ax.pcolor(Rs*unitR+Rmin, Zs*unitZ+Zmin, z, cmap='RdBu')
    plt.title('Shear Rate')
    #plt.title(r'Reynolds Tur. Stress: $\Pi_{\theta\theta}$')
    #plt.title(r'Reynolds Tur. Poloidal Force: $-\partial_{\psi}\tilde{\Pi_{\psi\theta}}-\partial_{\theta}\tilde{\Pi_{\theta\theta}}$')
    #plt.title(r'Reynolds Radial Force: $\partial_{\psi}\Pi_{\theta\psi}$')
    #plt.title(r'Poloidal Flow: $u_{\theta}$')
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.colorbar(im,format='%.0e')
    plt.show()

def sh_on_sep(arr):
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    Zs = np.array([x for x in range(0,arr.shape[0])])
    Rs = np.array([x for x in range(0,arr.shape[1])])
    file = open("shear_file.txt",'a')
    #file = open("Reynolds_turpf_mid.txt",'a')
    #file = open("Reynolds_force_midplane.txt",'a')
    #file = open("pol_flow_mid.txt",'a')
    #file = open("Re_rad_f_mid.txt",'a')
    #file = open("Re_rad_f_mid.txt",'a')
    Septrix = list(zip(Rpol,Zpol))
    Z_coord = []
    Z_loc = []
    R_loc = []
    R_coord = []
    for Z in Zpol:
        diff=[]
        for zeta in Zs:
            diff.append(abs(Z-(zeta*unitZ+Zmin)))
        tmp = np.amin(diff)
        Sep_loc_z = diff.index(tmp)
        Z_loc.append(Sep_loc_z)
        Z_coord.append(Zs[Sep_loc_z]*unitZ+Zmin)
    for R in Rpol:
        diff=[]
        for rho in Rs:
            diff.append(abs(R-(rho*unitR+Rmin)))
        Sep_loc_r = np.argmin(diff)
        R_loc.append(Sep_loc_r)
        R_coord.append(Rs[Sep_loc_r]*unitR+Rmin)
    #plt.plot(R_coord,Z_coord,'r--',Rpol,Zpol,'b*')
    #plt.show()
    #file.write("R"+"\t"+"Z"+"\t"+"Shear"+"\n")
    #file.write("R"+"\t"+"Z"+"\t"+"Re_S"+"\n")
    file.write("R"+"\t"+"Z"+"\t"+"Re_F"+"\n")
    coordinates = list(zip(R_coord,Z_coord))
    for i in range(0,len(coordinates)):
        file.write(str(coordinates[i][0])+"\t"+str(coordinates[i][1])+"\t"+str(arr[Z_loc[i],R_loc[i]])+"\n")
    file.close()
    
def shear_read():
    sh_file = open("/global/homes/g/giannos/xgc_python_dir/shear_file.txt","r")
    #sh_file = open("/global/homes/g/giannos/xgc_python_dir/Reynolds_turpf_mid.txt","r")
    #sh_file = open("/global/homes/g/giannos/xgc_python_dir/Reynolds_force_midplane.txt","r")
    #sh_file = open("/global/homes/g/giannos/xgc_python_dir/pol_flow_mid.txt","r")
    #sh_file = open("/global/homes/g/giannos/xgc_python_dir/Re_rad_f_mid.txt","r")
    next(sh_file)
    R=[]
    Z=[]
    S=[]
    for columns in (raw.strip().split() for raw in sh_file):
        R.append(columns[0])
        Z.append(columns[1])
        S.append(columns[2])
    ang = []
    for i in range(0,len(R)):
        ang.append(np.degrees(np.arctan2(float(Z[i]),float(R[i]))))
    plt.title("Shear on Separatrix")
    #plt.title("Reynolds Stress on Separatrix")
    #plt.title("Reynolds Turbulent Poloidal Force on Separatrix")
    #plt.title("Poloidal Flow on Separatrix")
    plt.xlabel("poloidal angle")
    plt.ylabel("Shear")
    #plt.ylabel('Stress')
    #plt.ylabel(r'$u_{\theta}$')
    #plt.ylabel("Force")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid()
    plt.plot(ang,S)
    plt.show()
    
    
    
    
    
    
# utility routine to locate psii indices (i.e. ip values) inside and outside the separatrix
def getiparrs_hor(psimin,psimax,cut,verbose=False):
	global chc,psiichc,Richc
	if cut == None:
		cut = jcut
    # calculate psi(Ri) = psii at the cut
	tci=LinearTriInterpolator(triObj,psin)
	psiout=tci(RI,ZI)
	psii = psiout[cut,:] # a 1D array of psin values at each Ri along the cut

	chc1 = np.where(psii > psimin)[0]
	chc2 = np.where(psii <= psimax)[0]
	chc = np.intersect1d(chc1,chc2) # to invoke both conditions
	psiichc = psii[chc].filled()
	Richc = Ri[chc]
	if verbose:
		for i in range(len(chc)):
			print ("%d   %4f   %4f" % (chc[i],psiichc[i],Richc[i]))
	return chc

# utility routine to locate psii indices (i.e. ip values) inside and outside the separatrix/on vertical direction
def getiparrs_ver(psimin,psimax,cut,verbose=False):
	global chc,psiichc,Richc
	if cut == None:
		cut = rcut
    # calculate psi(Ri) = psii at the cut
	tci=LinearTriInterpolator(triObj,psin)
	psiout=tci(RI,ZI)
	psiz = psiout[:,cut] # a 1D array of psin values at each Ri along the cut

	chc1 = np.where(psiz > psimin)[0]
	chc2 = np.where(psiz <= psimax)[0]
	chc = np.intersect1d(chc1,chc2) # to invoke both conditions
	psiichc = psiz[chc].filled()
	Zichc = Zi[chc]
	if verbose:
		for i in range(len(chc)):
			print ("%d   %4f   %4f" % (chc[i],psiichc[i],Zichc[i]))
	return chc

def magnetic_drifts(cut,ip_value,start):
    global Zi,Ri,pot,bfield,B2
    option = 1
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
    #Flux_ip = getiparrs_hor(0.95,1.,cut)[-2]
    #psi_val = getcutvalue_hor(psin[:],cut,option)[Flux_ip]
    oned_location = (np.where(psin1d>=0.99))[0][0]
    Tigc1d = loader.Ti1d
    Tigcperp = loader.Tiperp
    Tigcpar = loader.Tipar
    Ti_prp_avg = np.mean(Tigcperp[start:,oned_location]) 
    Ti_prl_avg = 2*np.mean(Tigcpar[start:,oned_location]) 
    Ti_avg = np.mean(Tigc1d[start:,oned_location])
    V_th = np.sqrt((joule_conv*Ti_avg)/m_i)
    print('V_th:', '%.2E' % Decimal(V_th))
    print('Temp_ratio:','%.2E' % Decimal(Ti_prp_avg/Ti_prl_avg))
    V_perp_2 = (2*joule_conv*Ti_prp_avg)/m_i
    V_par_2 = (joule_conv*Ti_prl_avg)/m_i
    print('V_perp_sq=','%.2E' % Decimal(V_perp_2))
    print('V_par_sq=','%.2E' % Decimal(V_par_2))
    #Loading B-field 
    #bfield=loader.loadBfield()
    #calculate values at the cut and find the gyrofrequency
    Br = getcutvalue_hor(bfield[:,0],cut,option)
    Bz = getcutvalue_hor(bfield[:,1],cut,option)
    Bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    BfMag=[]
    for i in range(len(Br)):
        BfMag.append(np.sqrt(Br[i]*Br[i] + Bz[i]*Bz[i] + Bzeta[i]*Bzeta[i]))
    BfMag = np.asarray(BfMag)
    print('B:','%.2E' % Decimal(BfMag[ip_value]))
    Omega = (charge*BfMag/(m_i))
    print('Omega:',np.ma.masked_invalid(Omega).mean(),np.std(np.ma.masked_invalid(Omega)))
    print('Omega:','%.2E' % Decimal(Omega[ip_value]))
    print('rho:','%.2E' % Decimal(V_th/Omega[ip_value]))

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
    Brm1 = getcutvalue_hor(bfield[:,0],cut-1,option)
    Bzm1 = getcutvalue_hor(bfield[:,1],cut-1,option)
    Bzetam1 = getcutvalue_hor(bfield[:,2],cut-1,option)
    BfMagm1=[]
    for i in range(len(Brm1)):
        BfMagm1.append(np.sqrt(Brm1[i]*Brm1[i] + Bzm1[i]*Bzm1[i] + Bzetam1[i]*Bzetam1[i]))
    BMm1 = np.asarray(BfMagm1)
    brm1 = Brm1/BMm1
    bzm1 = Bzm1/BMm1
    bzetam1 = Bzetam1/BMm1
    Br = getcutvalue_hor(bfield[:,0],cut,option)
    Bz = getcutvalue_hor(bfield[:,1],cut,option)
    Bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    BfMag=[]
    for i in range(len(Br)):
        BfMag.append(np.sqrt(Br[i]*Br[i] + Bz[i]*Bz[i] + Bzeta[i]*Bzeta[i]))
    BM = np.asarray(BfMag)
    br = Br/BM
    bz = Bz/BM
    bzeta = Bzeta/BM
    Brp1 = getcutvalue_hor(bfield[:,0],cut+1,option)
    Bzp1 = getcutvalue_hor(bfield[:,1],cut+1,option)
    Bzetap1 = getcutvalue_hor(bfield[:,2],cut+1,option)
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
    lengthR = len(getcutvalue_hor(bfield[:,0],cut,option))
    iterlistR = list(range(0,lengthR))
    Brrm10 = getcutvalue_hor(bfield[:,0],cut,option)[lengthR-1]
    Bzrm10 = getcutvalue_hor(bfield[:,1],cut,option)[lengthR-1]
    Bzetarm10 = getcutvalue_hor(bfield[:,2],cut,option)[lengthR-1]
    BfMagrm10=np.sqrt(Brrm10*Brrm10 + Bzrm10*Bzrm10 + Bzetarm10*Bzetarm10)
    brrm10 = Brrm10/BfMagrm10
    bzrm10 = Bzrm10/BfMagrm10
    bzetarm10 = Bzetarm10/BfMagrm10
    Brrp10 = getcutvalue_hor(bfield[:,0],cut,option)[1]
    Bzrp10 = getcutvalue_hor(bfield[:,1],cut,option)[1]
    Bzetarp10 = getcutvalue_hor(bfield[:,2],cut,option)[1]
    BfMagrp10=np.sqrt(Brrp10*Brrp10 + Bzrp10*Bzrp10 + Bzetarp10*Bzetarp10)
    brrp10 = Brrp10/BfMagrp10
    bzrp10 = Bzrp10/BfMagrp10
    bzetarp10 = Bzetarp10/BfMagrp10
    dbrdR.append((brrp10-brrm10)/(2.*dR))
    dbzdR.append((bzrp10-bzrm10)/(2.*dR))
    dbzetadR.append((bzetarp10-bzetarm10)/(2.*dR))
    dBdR.append((BfMagrp10-BfMagrm10)/(2.*dR))
    for i in iterlistR[1:-1]:
        Brrm1 = getcutvalue_hor(bfield[:,0],cut,option)[i-1]
        Bzrm1 = getcutvalue_hor(bfield[:,1],cut,option)[i-1]
        Bzetarm1 = getcutvalue_hor(bfield[:,2],cut,option)[i-1]
        BMrm1 = np.sqrt(Brrm1*Brrm1 + Bzrm1*Bzrm1 + Bzetarm1*Bzetarm1)
        brrm1 = Brrm1/BMrm1
        bzrm1 = Bzrm1/BMrm1
        bzetarm1 = Bzetarm1/BMrm1
        Brrp1 = getcutvalue_hor(bfield[:,0],cut,option)[i+1]
        Bzrp1 = getcutvalue_hor(bfield[:,1],cut,option)[i+1]
        Bzetarp1 = getcutvalue_hor(bfield[:,2],cut,option)[i+1]
        BMrp1=np.sqrt(Brrp1*Brrp1 + Bzrp1*Bzrp1 + Bzetarp1*Bzetarp1)
        brrp1 = Brrp1/BMrp1
        bzrp1 = Bzrp1/BMrp1
        bzetarp1 = Bzetarp1/BMrp1
        dbrdR.append((brrp1-brrm1)/(2.*dR))
        dbzdR.append((bzrp1-bzrm1)/(2.*dR))
        dbzetadR.append((bzetarp1-bzetarm1)/(2.*dR))
        dBdR.append((BMrp1-BMrm1)/(2.*dR))
    Brrm1L = getcutvalue_hor(bfield[:,0],cut,option)[lengthR-2]
    Bzrm1L = getcutvalue_hor(bfield[:,1],cut,option)[lengthR-2]
    Bzetarm1L = getcutvalue_hor(bfield[:,2],cut,option)[lengthR-2]
    BfMagrm1L=np.sqrt(Brrm1L*Brrm1L + Bzrm1L*Bzrm1L + Bzetarm1L*Bzetarm1L)
    brrm1L = Brrm1L/BfMagrm1L
    bzrm1L = Bzrm1L/BfMagrm1L
    bzetarm1L = Bzetarm1L/BfMagrm1L
    Brrp1L = getcutvalue_hor(bfield[:,0],cut,option)[0]
    Bzrp1L = getcutvalue_hor(bfield[:,1],cut,option)[0]
    Bzetarp1L = getcutvalue_hor(bfield[:,2],cut,option)[0]
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





#routines for sensitivity to grid dimension
def R_grid_ch(pn):
    VZ_all = np.asarray(getcutvRvZ_old(0,120,jcut)[1])
    Vmask = np.isfinite(VZ_all)
    L = np.array([x for x in range(0,pn-1)])
    unit = ((Rmax-Rmin)/pn)
    R = np.array([L[ip]*unit+Rmin for ip in range(0,pn-1)])
    plt.title("VZ vs R evaluated on a %dX100 mesh" % (pn))
    plt.xlabel('$R(m)$')
    plt.ylabel('$V_Z$ $(m/s)$')
    plt.plot(R[Vmask],VZ_all[Vmask],'go')
    plt.show()

def Z_grid_ch(pn):
    VR = []
    for i in range(0,pn-1):
        VR_all = np.asarray(getcutvRvZ_old(0,120,i)[0])
        VR.append(VR_all[rcut])
    #print(VR)
    VR = np.asarray(VR)
    Vmask = np.isfinite(VR)
    #print(Vmask)
    L = np.array([x for x in range(0,pn-1)])
    unit = ((Zmax-Zmin)/pn)
    Z = np.array([L[ip]*unit+Zmin for ip in range(0,pn-1)])
    plt.title("VR vs Z evaluated on a 100X%d mesh" % (pn))
    plt.xlabel('$Z(m)$')
    plt.ylabel('$V_R$ $(m/s)$')
    plt.plot(Z[Vmask],VR[Vmask],'go')
    plt.show()



#diamagnetic frequency

def diamagnetic_frequency():
    Rp = np.roll(psin1d,1)
    Rm = np.roll(psin1d,-1)
    dRp = psin1d-Rp
    dRm = Rm - psin1d
    
    ne_avg = ne1d.mean(axis=0)
    
    dndpsi = []
    for i in range(1,len(ne_avg)-1):
        A = dRp[i+1]/(dRm[i-1]**2 + dRp[i+1]*dRm[i-1]) - dRm[i-1]/(dRp[i+1]**2 + dRp[i+1]*dRm[i-1])
        B = -dRp[i+1]/(dRm[i-1]**2 + dRp[i+1]*dRm[i-1])
        C = dRm[i-1]/(dRp[i+1]**2 + dRp[i+1]*dRm[i-1])
        dndpsi.append(A*ne_avg[i] + B*ne_avg[i-1] + C*ne_avg[i+1])
    dndpsi.append((ne_avg[-1]-ne_avg[-2])/dRp[-1])
    dndpsi = [(ne_avg[1]-ne_avg[0])/dRm[0]] + dndpsi
    
    dx = 1
    psiplus = np.roll(psin1d,-1)
    psiminus = np.roll(psin1d,1)
    dpsidx = (psiplus-psiminus)/(2.*dx)
    dpsidx[0] = (psin1d[1]-psin1d[0])/dx
    dpsidx[-1] = (psin1d[-1]-psin1d[-2])/dx

    R_list = RZ[:,0]#call it with Rmax=None
    Rmax = np.max(R_list)
    unitR = Rmax/len(psin1d)
    
    dndR = dndpsi*dpsidx*(1/unitR)
    
    L = [x for x in range(0,len(psin1d))]
    
    Ln = -(ne_avg)*(1/dndR)
    
    m_i = 2*1.6723e-27
    k = 1.38e-23
    joule_conv = 1.602e-19
    Te_avg = Te1d.mean(axis=0)
    speed = np.sqrt((joule_conv*Te_avg)/(2*m_i))
    
    
    
    diamagnetic_fr = (0.2)*(speed/abs(Ln))*(1./2.*math.pi)/1000 #divide to make it khz
    
    '''fig, ax1 = plt.subplots()
    plt.title('Diamagnetic Frequency')
    ax1.plot(L,diamagnetic_fr,'r-')
    ax1.set_xlabel('x')
    ax1.set_ylabel('frequency (kHz)',color='r')
    ax2 = ax1.twinx()
    ax2.plot(L,psin1d,'b-')
    ax2.set_ylabel('$\psi$',color='b')
    fig.tight_layout()
    plt.show()'''
    
    plt.title('Diamagnetic Frequency vs. $\psi$')
    plt.plot(psin1d,diamagnetic_fr,'r-')
    plt.xlabel('$\psi$',fontsize=18)
    plt.ylabel('$\omega_{\star}$ (kHz)',fontsize=18)
    plt.show()


def save_on_cut(arr):
    arr_all = np.array([[getcutvalue_hor(arr[:,iz,it],jcut,1) for iz in range(Nplanes)] for it in range(490,690)])
    arr_tor = arr_all.mean(axis=1)
    arr_avg = arr_tor.mean(axis=0)
    np.save('argument',arr_avg)
    

def potential_around_sep():
    Rpol = RZ[sepInds[91:1271],0]
    Zpol = RZ[sepInds[91:1271],1]
    dRpol = Rpol-Rpol.mean()
    dZpol = Zpol-Zpol.mean()
    thetapol = angle.norm_atan(dZpol,dRpol)#np.arctan2(dZpol,dRpol)#do it with angle.norm_atan
    phipol = pot[sepInds[91:1271],:,:]
    phi_tor = phipol.mean(axis=1)
    phi_tt_avg = phi_tor.mean(axis=1)
    phi_avg = phi_tt_avg.mean()
    phinorm = phi_tt_avg/phi_avg
    n_e=loader.calcNeTotal()
    npol = n_e[sepInds[91:1271],:,:]
    n_tor = npol.mean(axis=1)
    n_tt_avg = n_tor.mean(axis=1)
    n_avg = n_tt_avg.mean()
    nnorm = n_tt_avg/n_avg
    #plt.plot(thetapol,phinorm,'r-',thetapol,nnorm,'b-')
    fig, ax1 = plt.subplots()
    ax1.set_xlim([0,360])
    ax1.plot(thetapol[0:10],n_tt_avg[0:10],'b.')
    ax1.plot(thetapol[10:],n_tt_avg[10:],'b.')
    ax1.set_xlabel('angle')
    ax1.set_ylabel('$n_e$',color='b')
    ax2 = ax1.twinx()
    ax2.plot(thetapol,phi_tt_avg,'r.')
    ax2.set_xlim([0,360])
    ax2.set_ylabel('$\Phi$',color='r')
    fig.tight_layout()
    plt.show()

def sep_plot():
    Rpol = RZ[sepInds[91:1271],0]
    Zpol = RZ[sepInds[91:1271],1]
    plt.plot(Rpol,Zpol)
    plt.show()
    
def parallel_velocity():
    Rpol = RZ[sepInds[91:1271],0]
    Zpol = RZ[sepInds[91:1271],1]
    dRpol = Rpol-Rpol.mean()
    dZpol = Zpol-Zpol.mean()
    thetapol = angle.norm_atan(dZpol,dRpol)#np.arctan2(dZpol,dRpol)#do it with angle.norm_atan
    phipol = pot[sepInds[91:1271],:,:]
    phi_tor = phipol.mean(axis=1)
    phi_tt_avg = phi_tor.mean(axis=1)
    phi_avg = phi_tt_avg.mean()
    phinorm = phi_tt_avg/phi_avg
    npol = i_u_para[sepInds[91:1271],:,:]
    n_tor = npol.mean(axis=1)
    n_tt_avg = n_tor.mean(axis=1)
    n_avg = n_tt_avg.mean()
    nnorm = n_tt_avg/n_avg
    #plt.plot(thetapol,phinorm,'r-',thetapol,nnorm,'b-')
    fig, ax1 = plt.subplots()
    ax1.set_xlim([0,360])
    ax1.plot(thetapol[0:10],n_tt_avg[0:10],'b.')
    ax1.plot(thetapol[10:],n_tt_avg[10:],'b.')
    ax1.set_xlabel('angle')
    ax1.set_ylabel('$u_\parallel$',color='b')
    ax2 = ax1.twinx()
    ax2.plot(thetapol,phi_tt_avg,'r.')
    ax2.set_xlim([0,360])
    ax2.set_ylabel('$\Phi$',color='r')
    fig.tight_layout()
    plt.show()
    
def separatrix():
    Rpol = RZ[sepInds[91:1271],0]
    Zpol = RZ[sepInds[91:1271],1]
    plt.plot(Rpol,Zpol)
    plt.show()
    
#reads two files with data and plots them with common x-axis    
def double_plot():
    file1=open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/sep_flux/short_run/ti255_Short_fluxes_total.txt","r")
    file2=open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/Reynolds/Re_tur_pol_f_total.txt","r")
    next(file1)
    next(file2)
    x1=[]
    x2=[]
    y1=[]
    y2=[]
    for columns in (raw.strip().split() for raw in file1):  
            x1.append(float(columns[0]))
            y1.append(float(columns[2]))
    for columns in (raw.strip().split() for raw in file2):  
            x2.append(float(columns[0]))
            y2.append(float(columns[1]))
    fig, ax1 = plt.subplots()
    ax1.set_xlim([0,360])
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax1.plot(x1,y1,'b.')
    #ax1.plot(thetapol[10:],n_tt_avg[10:],'b.')
    ax1.set_xlabel('angle')
    ax1.set_ylabel('Turbulent Flux',color='b')
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    ax2.plot(x2,y2,'r.')
    ax2.set_xlim([0,360])
    ax2.set_ylabel('Turbulent Poloidal Reynolds Force',color='r')
    fig.tight_layout()
    ax1.grid()
    plt.show()
    fig.savefig("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/figures/Re_tur_pol_vs_tur_flux.png")

#function that calculates the electron ExB and parallel heat fluxes, equilibrium and turbulent. Also returns Z-components.    
def e_heat_flux(cut,Flux_ip,Angle,id_num,Rmin):
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start#add temperatures
    #Inputs:
    option=1
    time_range=200 #ti255_shorter:149,ti255_short:194, ti255_long:259, ti253:125
    start=0 #start where linear phase ends:ti255_shorter:490, ti255_short:445, ti255_long:380, ti253:200
    ts_max=30 #30 looks like the optimum
    n_e=loader.calcNeTotal()
    temp = np.array([[getcutvRvZ(iz,it,cut) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    VR_all = temp[:,:,0,Flux_ip]
    VZ_all = temp[:,:,1,Flux_ip]
    print("VR_all done...")
    print("VZ_all done...")
    ne_all = np.array([[getcutvalue_hor(n_e[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Te_perp_all = np.array([[angle.map_array_rect(e_T_perp,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ee_para_all = np.array([[angle.map_array_rect(e_E_para,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Vpar_all = np.array([[angle.map_array_rect(e_u_para,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Temperature = 2/3.*(Te_perp_all + Ee_para_all) 
    Pe_all = np.multiply(ne_all,Temperature)
    br = getcutvalue_hor(bfield[:,0],cut,option)
    bz = getcutvalue_hor(bfield[:,1],cut,option)
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    Pe_tor_avg = Pe_all.mean(axis=1)
    VR_tor_avg = VR_all.mean(axis=1)
    VZ_tor_avg = VZ_all.mean(axis=1)
    Vpar_tor_avg = Vpar_all.mean(axis=1)
    l=1
    m=0
    Pe_m = []
    VR_m = []
    VZ_m = []
    Vpar_m = []
    for i in range(0,ts_max):
        Pe_m.append(np.mean(Pe_tor_avg[i-m:i+ts_max]))
        VR_m.append(np.mean(VR_tor_avg[i-m:i+ts_max]))
        VZ_m.append(np.mean(VZ_tor_avg[i-m:i+ts_max]))
        Vpar_m.append(np.mean(Vpar_tor_avg[i-m:i+ts_max]))
        m=m+1
    for i in range(ts_max,time_range-ts_max):#center
        Pe_m.append(np.mean(Pe_tor_avg[i-ts_max:i+ts_max]))
        VR_m.append(np.mean(VR_tor_avg[i-ts_max:i+ts_max]))
        VZ_m.append(np.mean(VZ_tor_avg[i-ts_max:i+ts_max]))
        Vpar_m.append(np.mean(Vpar_tor_avg[i-ts_max:i+ts_max]))
    for i in range(time_range-ts_max,time_range):#right side
        Pe_m.append(np.mean(Pe_tor_avg[i-ts_max-l:i+ts_max-l]))
        VR_m.append(np.mean(VR_tor_avg[i-ts_max-l:i+ts_max-l]))
        VZ_m.append(np.mean(VZ_tor_avg[i-ts_max-l:i+ts_max-l]))
        Vpar_m.append(np.mean(Vpar_tor_avg[i-ts_max-l:i+ts_max-l]))
        l=l+1
    Pe_avg = np.asarray(Pe_m)
    VR_avg = np.asarray(VR_m)
    VZ_avg = np.asarray(VZ_m)
    Vpar_avg = np.asarray(Vpar_m)
    dP = np.array(Pe_all[:,:] - Pe_avg[:,np.newaxis])
    dVR = np.array(VR_all[:,:] - VR_avg[:,np.newaxis])
    dVZ = np.array(VZ_all[:,:] - VZ_avg[:,np.newaxis])
    dVpar = np.array(Vpar_all[:,:] - Vpar_avg[:,np.newaxis])
    Bpol=np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:]))
    Bmag=np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:] + bzeta[:]*bzeta[:]))
    #Equilibrium ExB Flux
    Eq_R_flux = np.array(Pe_tor_avg[:]*VR_tor_avg[:])
    Eq_Z_flux = np.array(Pe_tor_avg[:]*VZ_tor_avg[:])        
    Eq_Gamma = np.array((-br[Flux_ip]*Eq_Z_flux[:] + bz[Flux_ip]*Eq_R_flux[:])/Bpol[Flux_ip])
    Eq_Gamma_time_avg = Eq_Gamma.mean()
    Eq_Z = Eq_Z_flux.mean()
    #Turbulent ExB Flux
    Tur_R_flux = np.array([[dP[it,pl]*dVR[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    Tur_Z_flux = np.array([[dP[it,pl]*dVZ[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    Tur_R_flux_tor_avg = Tur_R_flux.mean(axis=1)
    Tur_Z_flux_tor_avg = Tur_Z_flux.mean(axis=1)
    Tur_R_flux_time_avg = Tur_R_flux_tor_avg.mean()
    Tur_Z_flux_time_avg = Tur_Z_flux_tor_avg.mean()
    Tur_Gamma = np.array((-br[Flux_ip]*Tur_Z_flux_tor_avg[:] + bz[Flux_ip]*Tur_R_flux_tor_avg[:])/Bpol[Flux_ip])
    Tur_Gamma_time_avg = Tur_Gamma.mean()
    #Equilibrium parallel flux
    Eq_par = np.array(Pe_tor_avg[:]*Vpar_tor_avg[:])
    Eq_par_Z = Eq_par * (bz[Flux_ip]/Bmag[Flux_ip])
    Eq_par_ta = Eq_par.mean()
    Eq_par_Z_ta = Eq_par_Z.mean()
    #Turbulent parallel flux
    Tur_par_flux = np.array([[dP[it,pl]*dVpar[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    Tur_par_tor = Tur_par_flux.mean(axis=1)
    Tur_par_Z = Tur_par_tor*(bz[Flux_ip]/Bmag[Flux_ip])
    Tur_par = Tur_par_tor.mean()
    Tur_par_Z_ta = Tur_par_Z.mean()
    return Eq_Gamma_time_avg, Tur_Gamma_time_avg, Eq_Z, Tur_Z_flux_time_avg, Eq_par_ta, Eq_par_Z_ta, Tur_par, Tur_par_Z_ta

#function that calculates the ion ExB and parallel heat fluxes, equilibrium and turbulent. Also returns Z-components.    
def i_heat_flux(cut,Flux_ip,Angle,id_num,Rmin):
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start#add temperatures
    #Inputs:
    option=1
    time_range=200 #ti255_shorter:149,ti255_short:194, ti255_long:259, ti253:125
    start=0 #start where linear phase ends:ti255_shorter:490, ti255_short:445, ti255_long:380, ti253:200
    ts_max=30 #30 looks like the optimum
    n_e=loader.calcNeTotal()
    temp = np.array([[getcutvRvZ(iz,it,cut) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    VR_all = temp[:,:,0,Flux_ip]
    VZ_all = temp[:,:,1,Flux_ip]
    print("VR_all done...")
    print("VZ_all done...")
    ne_all = np.array([[getcutvalue_hor(n_e[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ti_perp_all = np.array([[angle.map_array_rect(i_T_perp,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ei_para_all = np.array([[angle.map_array_rect(i_E_para,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Vpar_all = np.array([[angle.map_array_rect(i_u_para,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Temperature = 2/3.*(Ti_perp_all + Ei_para_all) 
    Pi_all = np.multiply(ne_all,Temperature)
    br = getcutvalue_hor(bfield[:,0],cut,option)
    bz = getcutvalue_hor(bfield[:,1],cut,option)
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    Pi_tor_avg = Pi_all.mean(axis=1)
    VR_tor_avg = VR_all.mean(axis=1)
    VZ_tor_avg = VZ_all.mean(axis=1)
    Vpar_tor_avg = Vpar_all.mean(axis=1)
    l=1
    m=0
    Pi_m = []
    VR_m = []
    VZ_m = []
    Vpar_m = []
    for i in range(0,ts_max):
        Pi_m.append(np.mean(Pi_tor_avg[i-m:i+ts_max]))
        VR_m.append(np.mean(VR_tor_avg[i-m:i+ts_max]))
        VZ_m.append(np.mean(VZ_tor_avg[i-m:i+ts_max]))
        Vpar_m.append(np.mean(Vpar_tor_avg[i-m:i+ts_max]))
        m=m+1
    for i in range(ts_max,time_range-ts_max):#center
        Pi_m.append(np.mean(Pi_tor_avg[i-ts_max:i+ts_max]))
        VR_m.append(np.mean(VR_tor_avg[i-ts_max:i+ts_max]))
        VZ_m.append(np.mean(VZ_tor_avg[i-ts_max:i+ts_max]))
        Vpar_m.append(np.mean(Vpar_tor_avg[i-ts_max:i+ts_max]))
    for i in range(time_range-ts_max,time_range):#right side
        Pi_m.append(np.mean(Pi_tor_avg[i-ts_max-l:i+ts_max-l]))
        VR_m.append(np.mean(VR_tor_avg[i-ts_max-l:i+ts_max-l]))
        VZ_m.append(np.mean(VZ_tor_avg[i-ts_max-l:i+ts_max-l]))
        Vpar_m.append(np.mean(Vpar_tor_avg[i-ts_max-l:i+ts_max-l]))
        l=l+1
    Pi_avg = np.asarray(Pi_m)
    VR_avg = np.asarray(VR_m)
    VZ_avg = np.asarray(VZ_m)
    Vpar_avg = np.asarray(Vpar_m)
    dP = np.array(Pi_all[:,:] - Pi_avg[:,np.newaxis])
    dVR = np.array(VR_all[:,:] - VR_avg[:,np.newaxis])
    dVZ = np.array(VZ_all[:,:] - VZ_avg[:,np.newaxis])
    dVpar = np.array(Vpar_all[:,:] - Vpar_avg[:,np.newaxis])
    Bpol=np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:]))
    Bmag=np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:] + bzeta[:]*bzeta[:]))
    #Equilibrium ExB Flux
    Eq_R_flux = np.array(Pi_tor_avg[:]*VR_tor_avg[:])
    Eq_Z_flux = np.array(Pi_tor_avg[:]*VZ_tor_avg[:])        
    Eq_Gamma = np.array((-br[Flux_ip]*Eq_Z_flux[:] + bz[Flux_ip]*Eq_R_flux[:])/Bpol[Flux_ip])
    Eq_Gamma_time_avg = Eq_Gamma.mean()
    Eq_Z = Eq_Z_flux.mean()
    #Turbulent ExB Flux
    Tur_R_flux = np.array([[dP[it,pl]*dVR[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    Tur_Z_flux = np.array([[dP[it,pl]*dVZ[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    Tur_R_flux_tor_avg = Tur_R_flux.mean(axis=1)
    Tur_Z_flux_tor_avg = Tur_Z_flux.mean(axis=1)
    Tur_R_flux_time_avg = Tur_R_flux_tor_avg.mean()
    Tur_Z_flux_time_avg = Tur_Z_flux_tor_avg.mean()
    Tur_Gamma = np.array((-br[Flux_ip]*Tur_Z_flux_tor_avg[:] + bz[Flux_ip]*Tur_R_flux_tor_avg[:])/Bpol[Flux_ip])
    Tur_Gamma_time_avg = Tur_Gamma.mean()
    #Equilibrium parallel flux
    Eq_par = np.array(Pi_tor_avg[:]*Vpar_tor_avg[:])
    Eq_par_Z = Eq_par * (bz[Flux_ip]/Bmag[Flux_ip])
    Eq_par_ta = Eq_par.mean()
    Eq_par_Z_ta = Eq_par_Z.mean()
    #Turbulent parallel flux
    Tur_par_flux = np.array([[dP[it,pl]*dVpar[it,pl] for it in range(0,time_range)] for pl in range(Nplanes)])
    Tur_par_tor = Tur_par_flux.mean(axis=1)
    Tur_par_Z = Tur_par_tor*(bz[Flux_ip]/Bmag[Flux_ip])
    Tur_par = Tur_par_tor.mean()
    Tur_par_Z_ta = Tur_par_Z.mean()
    return Eq_Gamma_time_avg, Tur_Gamma_time_avg, Eq_Z, Tur_Z_flux_time_avg, Eq_par_ta, Eq_par_Z_ta, Tur_par, Tur_par_Z_ta

#function that calculates the ion magnetic heat fluxes. Returns Z-component.    
def mag_i_heat_flux(cut,Flux_ip,Angle,id_num,Rmin):
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start#add temperatures
    #Inputs:
    option=1
    time_range=200 #ti255_shorter:149,ti255_short:194, ti255_long:259, ti253:125
    start=0 #start where linear phase ends:ti255_shorter:490, ti255_short:445, ti255_long:380, ti253:200
    ts_max=30 #30 looks like the optimum
    n_e=loader.calcNeTotal()
    gradBR,gradBZ,_,curvR,curvZ,_ = magnetic_drifts(cut,Flux_ip,start)
    ne_all = np.array([[getcutvalue_hor(n_e[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ti_perp_all = np.array([[angle.map_array_rect(i_T_perp,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ei_para_all = np.array([[angle.map_array_rect(i_E_para,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Temperature = 2/3.*(Ti_perp_all + Ei_para_all) 
    Pi_all = np.multiply(ne_all,Temperature)    
    Pi_tor_avg = Pi_all.mean(axis=1)
    l=1
    m=0
    Pi_m = []
    for i in range(0,ts_max):
        Pi_m.append(np.mean(Pi_tor_avg[i-m:i+ts_max]))
        m=m+1
    for i in range(ts_max,time_range-ts_max):#center
        Pi_m.append(np.mean(Pi_tor_avg[i-ts_max:i+ts_max]))
    for i in range(time_range-ts_max,time_range):#right side
        Pi_m.append(np.mean(Pi_tor_avg[i-ts_max-l:i+ts_max-l]))
        l=l+1
    Pi_avg = np.asarray(Pi_m)
    #Magnetic flux
    Mag_R_flux = np.array(Pi_tor_avg[:]*(gradBR+curvR))
    Mag_Z_flux = np.array(Pi_tor_avg[:]*(gradBZ+curvZ))
    Mag_R = Mag_R_flux.mean()
    Mag_Z = Mag_Z_flux.mean()
    return Mag_R, Mag_Z

#function that calculates the ion magnetic heat fluxes. Returns Z-component.    
def mag_e_heat_flux(cut,Flux_ip,Angle,id_num,Rmin):
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start#add temperatures
    #Inputs:
    option=1
    time_range=200 #ti255_shorter:149,ti255_short:194, ti255_long:259, ti253:125
    start=0 #start where linear phase ends:ti255_shorter:490, ti255_short:445, ti255_long:380, ti253:200
    ts_max=30 #30 looks like the optimum
    n_e=loader.calcNeTotal()
    gradBR,gradBZ,_,curvR,curvZ,_ = magnetic_drifts(cut,Flux_ip,start)
    ne_all = np.array([[getcutvalue_hor(n_e[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Te_perp_all = np.array([[angle.map_array_rect(e_T_perp,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ee_para_all = np.array([[angle.map_array_rect(e_E_para,(iz-1)%Nplanes,it,cut)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Temperature = 2/3.*(Te_perp_all + Ee_para_all) 
    Pe_all = np.multiply(ne_all,Temperature)    
    Pe_tor_avg = Pe_all.mean(axis=1)
    l=1
    m=0
    Pe_m = []
    for i in range(0,ts_max):
        Pe_m.append(np.mean(Pe_tor_avg[i-m:i+ts_max]))
        m=m+1
    for i in range(ts_max,time_range-ts_max):#center
        Pe_m.append(np.mean(Pe_tor_avg[i-ts_max:i+ts_max]))
    for i in range(time_range-ts_max,time_range):#right side
        Pe_m.append(np.mean(Pe_tor_avg[i-ts_max-l:i+ts_max-l]))
        l=l+1
    Pe_avg = np.asarray(Pe_m)
    #Magnetic flux
    Mag_R_flux = np.array(Pe_tor_avg[:]*(gradBR+curvR))
    Mag_Z_flux = np.array(Pe_tor_avg[:]*(gradBZ+curvZ))
    Mag_R = Mag_R_flux.mean()
    Mag_Z = Mag_Z_flux.mean()
    return Mag_R, Mag_Z






def Par_quiver():
    option = 1
    minR = 1.001
    maxR = 2.377
    minZ = -1.363
    maxZ = 1.348
    
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    
    #getMeshAndCuts(fileDir,None,None,None,None)
    
    unitR = (maxR-minR)/100
    unitZ = (maxZ-minZ)/99
    
    V = np.load('/global/cscratch1/sd/giannos/ti255_par_vel.npy')
    V_tor = V.mean(axis=2)
    Vpar = V_tor.mean(axis = 1)
    R, Z = np.meshgrid(np.linspace(minR, maxR, num=100), np.linspace(minZ, maxZ, num=99))
    br = []
    bz = []
    bzeta = []
    for i in range(0,99):
        br.append(getcutvalue_hor(bfield[:,0],jcut,option))
        bz.append(getcutvalue_hor(bfield[:,1],jcut,option))
        bzeta.append(getcutvalue_hor(bfield[:,2],jcut,option))
    br = np.asarray(br)
    bz = np.asarray(bz)
    bzeta = np.asarray(bzeta)
    Bpol = []
    for i in range(0,99):
        Bpol.append(np.array(np.sqrt(br[i][:]*br[i][:] + bz[i][:]*bz[i][:])))
    Bpol = np.asarray(Bpol)
    VparR = []
    VparZ = []
    for z in range(0,99):
        for r in range(0,100):
            VparR.append(Vpar[z,r]*(br[z,r]/Bpol[z,r]))
            VparZ.append(Vpar[z,r]*(bz[z,r]/Bpol[z,r]))
    VparR = np.asarray(VparR)
    VparZ = np.asarray(VparZ)
    VparR = np.reshape(VparR,(99,100))
    VparZ = np.reshape(VparZ,(99,100))
    Vpol = []
    for z in range(0,99):
        for r in range(0,100):
            Vpol.append(np.sqrt(VparR[z,r]*VparR[z,r] + VparZ[z,r]*VparZ[z,r]))
    Vpol = np.asarray(Vpol)
    Vpol = np.reshape(Vpol,(99,100))
    
    fig,ax = plt.subplots()
    plt.xlabel('R(m)')
    plt.ylabel('Z(m)')
    plt.title("Projection of $V_{\parallel}$ on the poloidal plane")
    Q = ax.quiver(R,Z,VparR,VparZ,Vpol, scale = 1e6)
    ax.plot(Rpol,Zpol,'k--')
    plt.show()
    
def exb_quiv():
    minR = 1.001
    maxR = 2.377
    minZ = -1.363
    maxZ = 1.348
    
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    
    #getMeshAndCuts(fileDir,None,None,None,None)
    
    unitR = (maxR-minR)/100
    unitZ = (maxZ-minZ)/100
    
    V = np.load('/global/cscratch1/sd/giannos/temp.npy')
    VR = V[:,:,:,0,:]
    VZ = V[:,:,:,0,:]
    VR_tor = VR.mean(axis=1)
    VZ_tor = VZ.mean(axis=1)
    VR_ta = VR_tor.mean(axis=0)
    VZ_ta = VZ_tor.mean(axis=0)
    R, Z = np.meshgrid(np.linspace(minR, maxR, num=100), np.linspace(minZ, maxZ, num=100))
    Vpol = []
    for z in range(0,100):
        for r in range(0,100):
            Vpol.append(np.sqrt(VR_ta[z,r]*VR_ta[z,r] + VZ_ta[z,r]*VZ_ta[z,r]))
    Vpol = np.asarray(Vpol)
    Vpol = np.reshape(Vpol,(100,100))
    
    
    fig,ax = plt.subplots()
    plt.xlabel('R(m)')
    plt.ylabel('Z(m)')
    plt.title(" ExB Flow")
    Q = ax.quiver(R,Z,VR_ta,VZ_ta, Vpol, scale = 5e5)
    ax.plot(Rpol,Zpol,'k--')
    plt.show()
    
    
def temp_plot():
    global  triObj,tci,Sep_ip,start
    #Inputs:
    option=1
    cut=50
    time_range=2 #ti255_shorter:149,ti255_short:194, ti255_long:259, ti253:125
    start=0 #start where linear phase ends:ti255_shorter:490, ti255_short:445, ti255_long:380, ti253:200
    br = getcutvalue_hor(bfield[:,0],cut,1)
    temp_list = []
    psi_list = []
    for Flux_ip in range(0,100):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else: 
            Te_perp_all = np.array([[getcutvalue_hor(e_T_perp[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
            Ee_para_all = np.array([[getcutvalue_hor(e_E_para[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
            psi = np.array(getcutvalue_hor(psin[:],cut,option)[Flux_ip])
            Temperature = 2/3.*(Te_perp_all + Ee_para_all) 
            Temp_tor = Temperature.mean(axis=1)
            Temp = Temp_tor.mean()
            temp_list.append(Temp)
            psi_list.append(psi)
    fig,ax = plt.subplots()
    plt.xlabel('$\psi$')
    plt.ylabel('$T_e$')
    ax.plot(psi_list,temp_list,'r--')
    plt.show()
    
    
    