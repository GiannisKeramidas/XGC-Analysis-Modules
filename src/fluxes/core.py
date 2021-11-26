import numpy as np
import xgc
import math
import pylab as py
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.tri as mtri
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
import sys
import cProfile
import glob
import angle
import numpy.ma as ma
from decimal import Decimal
import os
#import ti255
import ti262
#import ti344
#plt.switch_backend('agg')

#set up parameters for plotting:
plt.style.use('ggplot')
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'heavy'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12
plt.rcParams['text.usetex']=True


'''List of things to check before each run:
   1) Import the right module of the discharge
   2) Change the options directly below. 
'''

Rmin = ti262.Rmin
Rmax = ti262.Rmax
Zmin = ti262.Zmin
Zmax = ti262.Zmax
phi_start = ti262.phi_start
phi_end = ti262.phi_end
Rmaj = ti262.Rmaj
fileDir = ti262.fileDir


def getMeshAndCuts(fileDir,Rmin,Rmax,Zmin,Zmax,nRi=None,nZi=None):

    global loader,ne,pot,psin,RZ,tri,time,psin1d,psin001d,Te1d,ne1d,pot001d,bfield,pot0,dpot,dpot,pot,Ti1d,Tigc1d,nigc1d#,i_T_perp,i_E_para,e_T_perp,e_E_para,i_u_para, e_u_para#, i_rad_fl, i_rad_exb_fl
    global Nplanes,Ntimes,dt
    global Ri,Zi,RI,ZI,triObj,unitR,unitZ
    global sepInds,psinvalues,jcut,Zcut,Rcut,rcut,fsInds,psii

    
    loader=xgc.load(fileDir,Rmin=Rmin,Rmax=Rmax,Zmin=Zmin,Zmax=Zmax,phi_start=phi_start,phi_end=phi_end)
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
    Ti1d = loader.Ti1d
    pot001d = loader.pot001d
    bfield = loader.bfield
    pot0 = loader.pot0
    dpot = loader.dpot
    pot = loader.pot
    eden = loader.eden#non-adiabatic part of the density
    #i_rad_fl = loader.i_rad_fl
    #i_rad_exb_fl = loader.i_rad_exb_fl
    #if hasattr(loader,'i_T_perp'):
    #i_u_para = loader.i_u_para
    #i_T_perp = loader.i_T_perp
    #i_E_para = loader.i_E_para
    #e_T_perp = loader.e_T_perp
    #e_E_para = loader.e_E_para
    #e_u_para = loader.e_u_para
    
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
    # get indices of RZ vertices on a flux surface
    fsInds = np.where(np.abs(loader.psin-1.01)<1e-4)[0]

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
        psii = psiout[jcut,:]# a 1D array of psin values at each Ri along the cut
        
        
        print(loader.Rmin)
        unitR = (loader.Rmax-loader.Rmin)/len(Ri)
        unitZ = (loader.Zmax-loader.Zmin)/len(Zi)
        
        
        
def getcutvalue_hor(arrRZ,j,option):
    """Calculates the values at Ri along Zcut of arrRZ where arr is a mesh array
    usage:
        netmpi = getcutvalue_hor(ne[:,iz,it]) # where iz is the toroidal plane, and it the time index
    """
    '''option=1 is for evaluating on the rectangular grid
       option=2 is for evaluating on the grid by spline-generated flux points, equidistant in poloidal angle.'''
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
    
def getframe(arrRZ):
    '''Returns a full frame at the preselected resolution without the need to interpolate at each new cut.'''
    global triObj,RI,ZI,arrout,tci,jcut,RI,ZI
    tci=LinearTriInterpolator(triObj,arrRZ)
    arrout=tci(RI,ZI)
    return(arrout[:,:])    
       
def getcutvRvZ(iz,it,cut):
    '''Calculates the R and Z components of the equilibrium ExB drift'''
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
    
    #vRi = (Bzetai/B2)*(dphidZ*(1.+Bzobzet2)+BrBzobzet2*dphidR)
    #vZi = -(Bzetai/B2)*(dphidR*(1.+Brobzet2)+BrBzobzet2*dphidZ)
    #actual equilibrium velocities
    vRi = (1/B2)*(Bzetai*dphidZ)
    vZi = -(1/B2)*(dphidR*Bzetai)
    
    return(vRi,vZi)

def map_array(array,iz,it,cut):
    '''Takes an array at time = it, on a plane = iz and returns it's value at the half plane forward assuming that the values stay
       constant along field lines. Used to map the temperature values of the f3d files. Works on the rectangular grid.
    '''
    global RA,ZA,pot,bfield,B2
    option = 1
    dzeta = (2*np.pi)/Nplanes
    
    dZ = Zi[1] - Zi[0]
    dR = Ri[1] - Ri[0]
    
    arrim1 = getcutvalue_hor(array[:,iz,it],cut-1,option)
    arri = getcutvalue_hor(array[:,iz,it],cut,option)
    arrip1 = getcutvalue_hor(array[:,iz,it],cut+1,option)
    darridZ = (arrip1-arrim1)/(2.*dZ)
    
    arriplus = np.roll(arri,-1)
    arriminus = np.roll(arri,1)
    darrdR = (arriplus-arriminus)/(2.*dR)
    darrdR[0] = (arri[1]-arri[0])/dR
    darrdR[-1] = (arri[-1]-arri[-2])/dR
    
    arriF = getcutvalue_hor(array[:,(iz+1)%Nplanes,it],cut,option)
    
    darrdRF = []
    arriFplus = np.roll(arriF,-1)
    arriFminus = np.roll(arriF,1)
    darrdRF = (arriFplus-arriFminus)/(2.*dR)
    darrdRF[0] = (arriF[1]-arriF[0])/dR
    darrdRF[-1] = (arriF[-1]-arriF[-2])/dR
    
    arrim1F = getcutvalue_hor(array[:,(iz+1)%Nplanes,it],cut-1,option)
    arrip1F = getcutvalue_hor(array[:,(iz+1)%Nplanes,it],cut+1,option)
    darridZF = (arrip1F-arrim1F)/(2.*dZ)
    
    BRi = getcutvalue_hor(bfield[:,0],cut,option)
    BZi = getcutvalue_hor(bfield[:,1],cut,option)
    Bzetai = getcutvalue_hor(bfield[:,2],cut,option)
    
    R = np.asarray([loader.Rmin + unitR*x for x in range(len(Ri))])
    
    Sigma = -R*((BRi/Bzetai)*darrdR - (BZi/Bzetai)*darridZ)
    SigmaF = -R*((BRi/Bzetai)*darrdRF - (BZi/Bzetai)*darridZF)
    arr_half = (arri+arriF)/2. + (1./4.)*(Sigma - SigmaF)*dzeta
    
    return arr_half


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
    t_start = 0
    t_end = 200
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

'''FUNCTIONS FOR THE CALCULATION AND IMAGING OF SHEAR. Important to make sure at which time step shear needs to be calculated.'''

def psi_der(arr,it,cut):#takes the d/dpsi on the potential and divides by RBp
    '''Takes a single derivative wrt psi along a cut.'''
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
    '''Produces a 2D matrix of psi derivatives.'''
    global Zi,Ri,pot,jcut,bfield,B2
    E_rad = [psi_der(pot,it,cut) for cut in range(0,len(Zi))]
    E_psi = np.asarray(E_rad)
    return E_psi

def shear(it,cut):#takes the d/dpsi of E_psi/RBp and multiplies with (RBp)^2/B to find the shear.
    '''Calculates the Burrell shear along a cut at a particular time.'''
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

def two_d_shear(it):
    '''Calculates the Burell shear on an RZ plane. Time point needs to be specified.'''
    global Zi,Ri,pot,jcut,bfield,B2
    Shear = [shear(it,cut) for cut in range(0,len(Zi))]
    Bur_shear = np.asarray(Shear)
    return Bur_shear

def shear_call_n_save(r):
    '''Function that saves the 2D shear values on an array. Gets called from parallel program.'''
    sh = two_d_shear(0)
    sh2 = np.nan_to_num(sh)
    np.save("shear_%s" %(r),sh2)

def sep_save(r):
    '''Function that saves the separatrix locations on arrays. Gets called from parallel program.''' 
    global RZ
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    np.save("Rsep_%s" %(r),Rpol)
    np.save("Zsep_%s" %(r),Zpol)
    
def two_d_plot(arr):
    '''Function that makes a 2D plot of an array. Usually shear but also can be used for any other 2D array to visualize on a plane.'''
    global RZ
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    arr = np.nan_to_num(arr)
    fig,ax = plt.subplots()
    ax.plot(Rpol,Zpol,'k--')
    dmove = 1
    Zs, Rs = np.mgrid[slice(0,arr.shape[0],dmove),slice(0,arr.shape[1],dmove)]
    z = arr[Zs,Rs]/(2*np.pi*1000)
    im = ax.pcolor(Rs*unitR+Rmin, Zs*unitZ+Zmin, z, cmap='RdBu_r')
    plt.title(r'$\Omega_{E\times B}$')
    #plt.title(r'Reynolds Tur. Stress: $\Pi_{\theta\theta}$')
    #plt.title(r'Reynolds Tur. Poloidal Force: $-\partial_{\psi}\tilde{\Pi_{\psi\theta}}-\partial_{\theta}\tilde{\Pi_{\theta\theta}}$')
    #plt.title(r'Reynolds Radial Force: $\partial_{\psi}\Pi_{\theta\psi}$')
    #plt.title(r'Poloidal Flow: $u_{\theta}$')
    plt.xlabel(r'$R (m)$')
    plt.ylabel(r'$Z (m)$')
    cbar = plt.colorbar(im,format='%.0e')
    cbar.set_label(r'kHz')
    plt.savefig("snew_d3d.png", bbox_inches='tight')
    plt.show()

def sh_on_sep2(arr):
    global RZ
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    arr = np.nan_to_num(arr)
    sep_loc = zip(int(Rpol))
    
def sh_on_sep(arr):
    '''Takes the shear array and extracts the shear values on the separatrix. Usually used from postpro module.'''
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
    '''Plotting from file.'''
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


#returns components of magnetic drifts. They need to be multiplied by squared velocities to be made into real drifts.
def magnetic_drifts(cut,ip_value, species):  
    '''Calculates all components of the magnetic drifts at a point. Select between ions and electrons.'''
    global Zi,Ri,pot,bfield,B2
    option = 1
    if cut==None:
        cut = jcut
    dZ = Zi[1]-Zi[0]
    dR = Ri[1]-Ri[0]
    theta = (2*math.pi)/Nplanes

    #Define constants, average temperatures and velocities
    k = 1.38e-23
    joule_conv = 1.602e-19
    permeability = 1.2566e-6
    
    #Loading B-field 
    #calculate values at the cut and find the gyrofrequency, gyroradius and magnetic moment.
    Br = getcutvalue_hor(bfield[:,0],cut,option)
    Bz = getcutvalue_hor(bfield[:,1],cut,option)
    Bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    BfMag=[]
    for i in range(len(Br)):
        BfMag.append(np.sqrt(Br[i]*Br[i] + Bz[i]*Bz[i] + Bzeta[i]*Bzeta[i]))
    BfMag = np.asarray(BfMag)
    #print('B:','%.2E' % Decimal(BfMag[ip_value]))

    if species == 'ion':
        #Constants
        m_i = 2*1.6723e-27
        charge = 1.602176e-19 
        #Temperatures are constant so, it suffices to take avg. temp. at a single time point.
        Ti_perp_all = np.array([map_array(i_T_perp,(iz-1)%Nplanes, 0,cut)[ip_value] for iz in range(Nplanes)])
        Ei_para_all = np.array([map_array(i_E_para,(iz-1)%Nplanes, 0,cut)[ip_value] for iz in range(Nplanes)])
        Temperature = 2/3.*(Ti_perp_all + Ei_para_all)
        Ti_avg = Temperature.mean()
        Ti_prp_avg = Ti_perp_all.mean()
        Ti_prl_avg = 2*Ei_para_all.mean()
    
        V_th = np.sqrt((joule_conv*Ti_avg)/m_i)
        #print('V_th:', '%.2E' % Decimal(V_th))
        #print('Temp_ratio:','%.2E' % Decimal(Ti_prp_avg/Ti_prl_avg))
        V_perp_2 = (joule_conv*Ti_prp_avg)/m_i
        V_par_2 = (joule_conv*Ti_prl_avg)/m_i
        #print('V_perp_sq=','%.2E' % Decimal(V_perp_2))
        #print('V_par_sq=','%.2E' % Decimal(V_par_2))
    
        Omega = np.asarray((charge*BfMag/(m_i)))
        #print('Omega:','%.2E' % Decimal(Omega[ip_value]))
        #print('rho:','%.2E' % Decimal(V_th/Omega[ip_value]))
        
    if species == 'electron':
        #Constants
        m_e = 9.1094e-31
        charge = -1.602176e-19 
        #Temperature profiles are constant so, it suffices to take avg. temp. at a single time point.
        Te_perp_all = np.array([map_array(e_T_perp,(iz-1)%Nplanes, 0,cut)[ip_value] for iz in range(Nplanes)])
        Ee_para_all = np.array([map_array(e_E_para,(iz-1)%Nplanes, 0,cut)[ip_value] for iz in range(Nplanes)])
        Temperature = 2/3.*(Te_perp_all + Ee_para_all)
        Te_avg = Temperature.mean()
        Te_prp_avg = Te_perp_all.mean()
        Te_prl_avg = 2*Ee_para_all.mean()
    
        V_th = np.sqrt((joule_conv*Te_avg)/m_e)
        #print('V_th:', '%.2E' % Decimal(V_th))
        #print('Temp_ratio:','%.2E' % Decimal(Te_prp_avg/Te_prl_avg))
        V_perp_2 = (joule_conv*Te_prp_avg)/m_e
        V_par_2 = (joule_conv*Te_prl_avg)/m_e
        #print('V_perp_sq=','%.2E' % Decimal(V_perp_2))
        #print('V_par_sq=','%.2E' % Decimal(V_par_2))
    
        Omega = np.asarray((charge*BfMag/(m_e)))
        #print('Omega:','%.2E' % Decimal(Omega[ip_value]))
        #print('rho:','%.2E' % Decimal(V_th/Omega[ip_value]))
        
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


    #Calculating the gradients of scalar |B| and dbr/dZ, dbr/dR, dbz/dZ, dbz/dR, dbzeta/dZ, dbzeta/dR

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
    br = np.asarray(Br/BM)
    bz = np.asarray(Bz/BM)
    bzeta = np.asarray(Bzeta/BM)
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
    dBdZ = np.asarray(dBdZ)
    dbrdZ = np.asarray(dbrdZ)
    dbzdZ = np.asarray(dbzdZ)
    dbzetadZ = np.asarray(dbzetadZ)

    # dB/dR
    dBdR=[]
    dbrdR=[]
    dbzdR=[]
    dbzetadR=[]
    
    brplus = np.roll(br,-1)
    bzplus = np.roll(bz,-1)
    bzetaplus = np.roll(bzeta,-1)
    BMplus = np.roll(BM,-1)
    
    brminus = np.roll(br,1)
    bzminus = np.roll(bz,1)
    bzetaminus = np.roll(bzeta,1)
    BMminus = np.roll(BM,1)
    
    dbrdR = (brplus-brminus)/(2.*dR)
    dbzdR = (bzplus-bzminus)/(2.*dR)
    dbzetadR = (bzetaplus-bzetaminus)/(2.*dR)
    dBdR = (BMplus-BMminus)/(2.*dR)
    
    dbrdR[0] = (br[1]-br[0])/dR
    dbzdR[0] = (bz[1]-bz[0])/dR
    dbzetadR[0] = (bzeta[1]-bzeta[0])/dR
    dBdR[0] = (BM[1]-BM[0])/dR
    
    dbrdR[-1] = (br[-1]-br[-2])/dR
    dbzdR[-1] = (bz[-1]-bz[-2])/dR
    dbzetadR[-1] = (bzeta[-1]-bzeta[-2])/dR
    dBdR[-1] = (BM[-1]-BM[-2])/dR
   
    dbrdR = np.asarray(dbrdR)
    dbzdR = np.asarray(dbzdR)   
    dbzetadR = np.asarray(dbzetadR)
    dBdR = np.asarray(dBdR)
    
    #dBdzeta
    dBdzeta = 0
    dbrdzeta = 0
    dbzdzeta = 0
    dbzetadzeta = 0
    L = np.array([x for x in range(0,200)])
    unit = (Rmax-Rmin)/200
    R = np.array([L[ip]*unit+Rmin for ip in range(0,len(br))])
    
    #calculation of b X gradB
    crossR = (-bz*(dBdzeta/R) + bzeta*dBdZ)
    crossZ = (-dBdR*bzeta + br*(dBdzeta/R))
    crosszeta = (dBdR*bz - dBdZ*br)
    
    #calculation of geometric portion of gradB drifts (they need to be multiplied by v_perp^2 to be made drifts), 1/Omega (b X gradB)/B
    gradBR = (1/Omega)*(crossR/BM)
    gradBZ = (1/Omega)*(crossZ/BM)
    gradBzeta = (1/Omega)*(crosszeta/BM)
     
    #real drifts. Disable if you only want geometric effect.
    gradBR = gradBR*V_perp_2
    gradBZ = gradBZ*V_perp_2
    gradBzeta = gradBzeta*V_perp_2
    
    
    #calculation of curvature drifts
    #b*delb
    #b*Del = br*d/dr+bz*d/dz+bzeta*d/dzeta
    #(b*Del)(br,bz,bzeta) = 
    #R= br*dbrdR+bz*dbrdZ+(bzeta/R)*dbrdzeta - (bzeta*bzeta)/R
    #Z =  br*dbzdR+bz*dbzdZ+(bzeta/R)*dbzdzeta
    #zeta = br*dbzetadR+bz*dbzetadZ+(bzeta/R)*dbzetadzeta +(bzeta*br)/R
    
    #(1/Omega)*bX(b*Del)b= (they need to be multiplied by v_par^2 to be made into drifts)
    curvR = (-bz*(br*dbzetadR+bz*dbzetadZ+(bzeta/R)*dbzetadzeta+(bzeta*br)/R)+bzeta*(br*dbzdR+bz*dbzdZ+(bzeta/R)*dbzdzeta))*(1/Omega)
    curvZ = (-bzeta*(br*dbrdR+bz*dbrdZ+(bzeta/R)*dbrdzeta - (bzeta*bzeta)/R)+br*(br*dbzetadR+bz*dbzetadZ+(bzeta/R)*dbzetadzeta+(bzeta*br)/R))*(1/Omega)
    curvzeta = (-br*(br*dbzdR+bz*dbzdZ+(bzeta/R)*dbzdzeta)+bz*(br*dbrdR+bz*dbrdZ+(bzeta/R)*dbrdzeta - (bzeta*bzeta)/R))*(1/Omega)
    
    #Real drifts. Disable if you only want geometric effect.
    curvR = curvR*V_par_2
    curvZ = curvZ*V_par_2
    curvzeta = curvzeta*V_par_2
    
    return gradBR[ip_value],gradBZ[ip_value],gradBzeta[ip_value],curvR[ip_value],curvZ[ip_value],curvzeta[ip_value]
   

def core_der(cut):
    global Zi,Ri,pot,bfield,B2
    #derivatives on rectangular mesh
    option = 1
    if cut==None:
        cut = jcut
    dZ = Zi[1]-Zi[0]
    dR = Ri[1]-Ri[0]
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
    br = np.asarray(Br/BM)
    bz = np.asarray(Bz/BM)
    bzeta = np.asarray(Bzeta/BM)
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
    dBdZ_c = np.asarray(dBdZ)
    dbrdZ_c = np.asarray(dbrdZ)
    dbzdZ_c = np.asarray(dbzdZ)
    dbzetadZ_c = np.asarray(dbzetadZ)

    # dB/dR
    dBdR=[]
    dbrdR=[]
    dbzdR=[]
    dbzetadR=[]
    
    brplus = np.roll(br,-1)
    bzplus = np.roll(bz,-1)
    bzetaplus = np.roll(bzeta,-1)
    BMplus = np.roll(BM,-1)
    
    brminus = np.roll(br,1)
    bzminus = np.roll(bz,1)
    bzetaminus = np.roll(bzeta,1)
    BMminus = np.roll(BM,1)
    
    dbrdR = (brplus-brminus)/(2.*dR)
    dbzdR = (bzplus-bzminus)/(2.*dR)
    dbzetadR = (bzetaplus-bzetaminus)/(2.*dR)
    dBdR = (BMplus-BMminus)/(2.*dR)
    
    dbrdR[0] = (br[1]-br[0])/dR
    dbzdR[0] = (bz[1]-bz[0])/dR
    dbzetadR[0] = (bzeta[1]-bzeta[0])/dR
    dBdR[0] = (BM[1]-BM[0])/dR
    
    dbrdR[-1] = (br[-1]-br[-2])/dR
    dbzdR[-1] = (bz[-1]-bz[-2])/dR
    dbzetadR[-1] = (bzeta[-1]-bzeta[-2])/dR
    dBdR[-1] = (BM[-1]-BM[-2])/dR
   
    dbrdR_c = np.asarray(dbrdR)
    dbzdR_c = np.asarray(dbzdR)   
    dbzetadR_c = np.asarray(dbzetadR)
    dBdR_c = np.asarray(dBdR)
    
    #Electric field
    potim1 = getcutvalue_hor(pot[:,0,0],(cut-1)%len(Zi),option)
    poti = getcutvalue_hor(pot[:,0,0],cut,option)
    potip1 = getcutvalue_hor(pot[:,0,0],(cut+1)%len(Zi),option)
    dphidZ = (potip1-potim1)/(2.*dZ)
    dphidZ_c = np.asarray(dphidZ)
    
    potiplus = np.roll(poti,-1)
    potiminus = np.roll(poti,1)
    dphidR = (potiplus-potiminus)/(2.*dR)
    dphidR[0] = (poti[1]-poti[0])/dR
    dphidR[-1] = (poti[-1]-poti[-2])/dR
    dphidR_c = np.asarray(dphidR)

    return dBdZ_c, dbrdZ_c, dbzdZ_c, dbzetadZ_c, dbrdR_c, dbzdR_c, dbzetadR_c, dBdR_c, dphidZ_c, dphidR_c 

    


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


    
def flux_surfaces():
    '''Function that plots flux surfaces.'''
    psi_list = []
    for cut in range(100):
        psi_list.append(getcutvalue_hor(psin[:],cut,1))
    psi_mat = np.reshape(psi_list,(100,100))
    R = np.array([ip*unitR+Rmin for ip in range(0,len(Ri))])
    Z = np.array([ip*unitZ+Zmin for ip in range(0,len(Zi))])
    fig, ax =plt.subplots(figsize=(10,10))
    contour_levels = np.arange(0.7,1.4,0.05)
    CS = plt.contour(R, Z, psi_mat,contour_levels)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(loader.RZ[sepInds,0],loader.RZ[sepInds,1],'k--')
    ax.set_aspect('equal')
    plt.savefig("d3d_flux_surf.png", bbox_inches='tight')
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
    
    
    Ln= -(ne_avg)*(1/dndR)
    
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

#function that saves on an npy matrix 
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
    Rpol = RZ[sepInds[206:1757],0]
    Zpol = RZ[sepInds[206:1757],1]
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
    
#functions that finds the (cut, ip_value) pairs of grid points that are closest to the separatrix nodes.     
def sep_list():
    Rsep = RZ[sepInds[:],0]
    Zsep = RZ[sepInds[:],1]
    sep = zip(Rsep,Zsep)
    
    LR = [x for x in range(len(Ri))]
    LZ = [x for x in range(len(Zi))]
    
    ip_values = []
    for i in range(len(sep)):
        diff_R = []
        for ip in range(len(LR)):
            diff_R.append(abs(sep[i][0] - (ip*unitR+Rmin)))
        tmp = np.amin(diff_R)
        R_loc = diff_R.index(tmp)
        ip_values.append(R_loc)
        
    cuts = []
    for i in range(len(sep)):
        diff_Z = []
        for cut in range(len(LZ)):
            diff_Z.append(abs(sep[i][1] - (cut*unitZ+Zmin)))
        tmp = np.amin(diff_Z)
        Z_loc = diff_Z.index(tmp)
        cuts.append(Z_loc)
        
    flux_points = zip(cuts,ip_values)
    return flux_points
    
    
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

    
def main_fluxes():
    flux_points = sep_list()#come in the form [(cuts, ip_values)]
    print("Number of points:", len(flux_points))
    file = open("ti255_Rect_Particle_fluxes_vs.angle_Mid.txt", 'a')
    for i in range(len(flux_points)):
        R_val = flux_points[i][1]*unitR+Rmin
        Z_val = flux_points[i][0]*unitZ+Zmin
        Ang = angle.norm_atan(Z_val,R_val)
        Eq,Tur,_,_ = particle_fluxes(flux_points[i][0],flux_points[i][1])
        file.write(str(Ang)+"\t"+str(Eq)+"\t"+str(Tur)+"\t"+str(R_val) + "\t"+ str(Z_val)+ "\n")
        print("Point %s Done!" %(i))
    file.close()
        
    
'''Functions that calculate the Normal-to-divertors fluxes.'''       
def particle_fluxes_norm(cut,Flux_ip,N,VZ,Vpar):
    '''Calculates the ion ExB and parallel particle fluxes. These are total fluxes, no distinction between equilibrium and turbulent.
       Returns only Z-components.
    '''
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start
    #Inputs:
    option=1
    
    VZ_all = VZ[:,:,Flux_ip]
    ne_all = N[:,:,Flux_ip]
    Vpar_all = Vpar[:,:,Flux_ip]
        
    br = getcutvalue_hor(bfield[:,0],cut,option)
    bz = getcutvalue_hor(bfield[:,1],cut,option)
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    B = np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:] + bzeta[:]*bzeta[:]))
    
    #Calculation of fluxes
    E_flux_Z = ne_all*VZ_all
    Par_flux = ne_all*Vpar_all
    
    E_Z_tor = E_flux_Z.mean(axis=1)
    E_Z = E_Z_tor.mean()
    
    Par_tor = Par_flux.mean(axis=1)
    Par = Par_tor.mean()
    Par_Z = Par * (bz[Flux_ip]/B[Flux_ip])
    
    return E_Z,Par_Z

def mta(arr, ts_max, time_range):
    '''Performs a moving time average of an array that is already toroidally averaged. Returns an array of averages the same size.'''
    l=1
    m=0
    arr_m = []
    for i in range(0,ts_max):
        arr_m.append(np.mean(arr[i-m:i+ts_max]))
        m=m+1
    for i in range(ts_max,time_range-ts_max):#center
        arr_m.append(np.mean(arr[i-ts_max:i+ts_max]))
    for i in range(time_range-ts_max,time_range):#right side
        arr_m.append(np.mean(arr[i-ts_max-l:i+ts_max-l]))        
        l=l+1
    arr_avg = np.asarray(arr_m)
    return arr_avg

def magnetic_flux_norm(cut,Flux_ip,N):
    '''Calculates the magnetic particle flux. Returns Z-component.'''
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start
    #Inputs:
    gradBR,gradBZ,_,curvR,curvZ,_ = magnetic_drifts(cut,Flux_ip,'electron')
    ne_all = N[:,:,Flux_ip]
        
    #Calculation of flux
    Mag_Z_flux = ne_all*(gradBZ+curvZ)
    Mag_Z_tor = Mag_Z_flux.mean(axis=1)
    Mag_Z = Mag_Z_tor.mean()
    
    return  Mag_Z
            
def e_heat_flux_norm(cut,Flux_ip,VZ,N,Te,Ee,Vpar):
    '''Calculates the electron ExB and parallel heat fluxes. Returns Z-components.'''
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start#add temperatures
    #Inputs:
    option=1
    
    VZ_all = VZ[:,:,Flux_ip]
    ne_all = N[:,:,Flux_ip]
    Te_perp_all = Te[:,:,Flux_ip]
    Ee_para_all = Ee[:,:,Flux_ip]
    Vpar_all = Vpar[:,:,Flux_ip]
    
    Temperature = 2/3.*(Te_perp_all + Ee_para_all) 
    Pe_all = np.multiply(ne_all,Temperature)
    
    br = getcutvalue_hor(bfield[:,0],cut,option)
    bz = getcutvalue_hor(bfield[:,1],cut,option)
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    B = np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:] + bzeta[:]*bzeta[:]))
    
    #Calculation of fluxes
    E_flux_Z = Pe_all*VZ_all
    Par_flux = Pe_all*Vpar_all
    
    E_Z_tor = E_flux_Z.mean(axis=1)
    E_Z = E_Z_tor.mean()
    
    Par_tor = Par_flux.mean(axis=1)
    Par = Par_tor.mean()
    Par_Z = Par * (bz[Flux_ip]/B[Flux_ip])
    
    return E_Z,Par_Z
    
def i_heat_flux_norm(cut,Flux_ip,VZ,N,Ti,Ei,Vpar):
    '''Calculates the ion ExB and parallel heat flux. Returns Z-components.'''
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start#add temperatures
    #Inputs:
    option=1
    
    VZ_all = VZ[:,:,Flux_ip]
    ne_all = N[:,:,Flux_ip]
    Ti_perp_all = Ti[:,:,Flux_ip]
    Ei_para_all = Ei[:,:,Flux_ip]
    Vpar_all = Vpar[:,:,Flux_ip]
    
    Temperature = 2/3.*(Ti_perp_all + Ei_para_all) 
    Pi_all = np.multiply(ne_all,Temperature)
    
    br = getcutvalue_hor(bfield[:,0],cut,option)
    bz = getcutvalue_hor(bfield[:,1],cut,option)
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    B = np.array(np.sqrt(br[:]*br[:] + bz[:]*bz[:] + bzeta[:]*bzeta[:]))
    
    #Calculation of fluxes
    E_flux_Z = Pi_all*VZ_all
    Par_flux = Pi_all*Vpar_all
    
    E_Z_tor = E_flux_Z.mean(axis=1)
    E_Z = E_Z_tor.mean()
    
    Par_tor = Par_flux.mean(axis=1)
    Par = Par_tor.mean()
    Par_Z = Par * (bz[Flux_ip]/B[Flux_ip])
    
    return E_Z,Par_Z
           
def mag_i_heat_flux_norm(cut,Flux_ip,N,Ti,Ei):
    '''Calculates the ion magnetic heat fluxes. Returns Z-component.'''
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start#add temperatures
    #Inputs:
    option=1
    
    gradBR,gradBZ,_,curvR,curvZ,_ = magnetic_drifts(cut,Flux_ip,'ion')
    ne_all = N[:,:,Flux_ip]
    Ti_perp_all = Ti[:,:,Flux_ip]
    Ei_para_all = Ei[:,:,Flux_ip]
    Temperature = 2/3.*(Ti_perp_all + Ei_para_all) 
    Pi_all = np.multiply(ne_all,Temperature)    
    
    #Calculation of flux
    Mag_flux_Z = Pi_all*(gradBZ+curvZ)
   
    Mag_Z_tor = Mag_flux_Z.mean(axis=1)
    Mag_Z = Mag_Z_tor.mean()
    
    return Mag_Z
        
def mag_e_heat_flux_norm(cut,Flux_ip,N,Te,Ee):
    '''Calculates the ion magnetic heat fluxes. Returns Z-component.'''
    global VR_all, ne_all, VZ_all, br, bz, triObj,tci,Sep_ip,start#add temperatures
    #Inputs:
    option=1
    
    n_e=loader.calcNeTotal()
    gradBR,gradBZ,_,curvR,curvZ,_ = magnetic_drifts(cut,Flux_ip,'electron')
    ne_all = N[:,:,Flux_ip]
    Te_perp_all = Te[:,:,Flux_ip]
    Ee_para_all = Ee[:,:,Flux_ip]
    Temperature = 2/3.*(Te_perp_all + Ee_para_all) 
    Pe_all = np.multiply(ne_all,Temperature)    
    
    #Calculation of flux
    Mag_flux_Z = Pe_all*(gradBZ+curvZ)
   
    Mag_Z_tor = Mag_flux_Z.mean(axis=1)
    Mag_Z = Mag_Z_tor.mean()
    
    return Mag_Z

def speed_of_sound():
    m_i = 1.6723e-27
    k = 1.38e-23
    joule_conv = 1.602e-19
    Flux_ip = getiparrs_hor(0.95,1.,jcut)[-2]
    psi_val = getcutvalue_hor(psin[:],jcut,1)[Flux_ip]
    oned_location = (np.where(psin1d>psi_val))[0][0]
    Te_avg = Te1d.mean(axis=0)
    T_e = Te_avg[oned_location-1:]
    T_e_mean = T_e.mean()
    speed = np.sqrt((joule_conv*T_e_mean)/(2*m_i))
    return speed   




def Par_quiver(opt):
    option = 1
    minR = 1.001
    maxR = 2.377
    minZ = -1.363
    maxZ = 1.348
    
    
    
    unitR = (maxR-minR)/100
    unitZ = (maxZ-minZ)/99
    
    getMeshAndCuts(fileDir,minR,maxR,minZ,maxZ)
    
    Rpol = RZ[sepInds[:],0]
    Zpol = RZ[sepInds[:],1]
    
    R_points = np.array([i*unitR + minR for i in range(0,100)])
    Z_points = np.array([i*unitZ + minZ for i in range(0,99)])
    
    #V = np.load('/global/cscratch1/sd/giannos/ti255_par_vel.npy')
    V = np.load('/global/homes/g/giannos/xgc_python_dir/ti255_e_par_vel.npy')
    V_tor = V.mean(axis=2)
    Vpar = V_tor.mean(axis = 1)
    R, Z = np.meshgrid(np.linspace(minR, maxR, num=100), np.linspace(minZ, maxZ, num=99))
    br = []
    bz = []
    bzeta = []
    for i in range(0,99):
        br.append(getcutvalue_hor(bfield[:,0],i,option))
        bz.append(getcutvalue_hor(bfield[:,1],i,option))
        bzeta.append(getcutvalue_hor(bfield[:,2],i,option))
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
    if opt == 'c':
        fig, ax =plt.subplots()
        contour_levels = np.arange(-1.6,1.6,0.1)
        CS = plt.contourf(R_points, Z_points, Vpar/speed_of_sound())
        CB = plt.colorbar(CS, extend='both')
        
        CS.cmap.set_under('white')
        plt.clabel(CS,inline=1,fmt='%.2E' ,fontsize=10)
        plt.xlabel("R(m)")
        plt.ylabel("Z(m)")
        ax.plot(Rpol,Zpol,'k--')
        plt.title(r"$V_{\parallel}$ poloidal section")
        plt.show()
    if opt == 'q':
        fig,ax = plt.subplots()
        plt.style.use('ggplot')
        plt.xlabel('R(m)')
        plt.ylabel('Z(m)')
        plt.title("Projection of electron $V_{\parallel}$ on the poloidal plane")
        Q = ax.quiver(R[::2],Z[::2],VparR[::2],VparZ[::2],Vpol[::2], scale = 1e6)
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
    global  triObj,tci,Sep_ip,start,psi_list
    #Inputs:
    option=1
    cut=50
    time_range=1 #ti255_shorter:149,ti255_short:194, ti255_long:259, ti253:125
    start=0 #start where linear phase ends:ti255_shorter:490, ti255_short:445, ti255_long:380, ti253:200
    br = getcutvalue_hor(bfield[:,0],cut,1)
    temp_list = []
    psi_list = []
    R_list = []
    
    for Flux_ip in range(0,100):
        if math.isnan(br[Flux_ip]) == True:
            pass
        else: 
            Te_perp_all = np.array([[getcutvalue_hor(e_T_perp[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes-30)] for it in range(start,start+time_range)])
            Ee_para_all = np.array([[getcutvalue_hor(e_E_para[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes-30)] for it in range(start,start+time_range)])
            psi = getcutvalue_hor(psin[:],cut,option)[Flux_ip]
            Temperature = 2/3.*(Te_perp_all + Ee_para_all) 
            Temp_tor = Temperature.mean(axis=1)
            Temp = Temp_tor.mean()
            temp_list.append(Temp)
            psi_list.append(psi)
    
    fig,ax = plt.subplots()
    plt.xlabel('$\psi$')
    plt.ylabel('$T_e$')
    ax.plot(psi_list[:32],temp_list[:32],'r--',label='HFS')
    ax.plot(psi_list[32:],temp_list[32:],'b--',label='LFS')
    ax.grid(True)
    ax.legend()
    plt.show()

#parallel velocity on a line on a cut    
def par_line():
    Vpar = getcutvalue_hor(i_u_para[:,0,0],50,1)
    unitR = (loader.Rmax-loader.Rmin)/100
    R = [ip*unitR+Rmin for ip in range(100)]
    fig=plt.figure()
    plt.xlabel("R(m)")
    h = plt.ylabel(r"$\frac{V_{\parallel}}{c_s}$")
    h.set_rotation(0)
    plt.plot(R,Vpar/speed_of_sound())
    plt.axvline(x=R[43],color='k')
    plt.show()
    
#gives time history of parallel velocity    
def Vpar_time_hist(ip):
    Vpar = np.array([getcutvalue_hor(i_u_para[:,0,it],50,1)[ip] for it in range(0,100)])
    time = [t for t in range(100)]
    fig=plt.figure()
    plt.xlabel("t")
    plt.ylabel(r"$V_{\parallel}$")
    plt.plot(time,Vpar)
    plt.show()

#function that gives local value of safety factor q. Also returns the value of epsilon.
def q_loc(ip,cut):
    option=1
    R = Rmin + ip*unitR
    Z = Zmin + cut*unitZ
    r = np.sqrt(Z**2 + (R-Rmaj)**2)
    epsilon = r/Rmaj
    br = getcutvalue_hor(bfield[:,0],cut,option)
    bz = getcutvalue_hor(bfield[:,1],cut,option)
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)
    Bpol = np.sqrt(br**2 + bz**2)
    q = (r/R)*(bzeta[ip]/Bpol[ip]) 
    tf = np.sqrt(epsilon)*(1.46-0.46*epsilon)#fraction of trapped particles
    return q, epsilon,tf

#calculates local values of gyrofrequency,thermal velocity and Larmor radius for ions.
def gf(ip,iz,it,cut):
    global Zi,Ri,pot,bfield,B2
    option = 1

    m_i = 2*1.6723e-27
    m_e = 9.1e-31
    k = 1.38e-23
    joule_conv = 1.6022e-19
    charge = 1.6022e-19 
    
    Ti_perp = map_array(i_T_perp,(iz-1)%Nplanes,it,cut)[ip] 
    Ei_para = map_array(i_E_para,(iz-1)%Nplanes,it,cut)[ip] 
    Ti = (2./3.)*(Ti_perp + Ei_para) 
    Te_perp = map_array(e_T_perp,(iz-1)%Nplanes,it,cut)[ip] 
    Ee_para = map_array(e_E_para,(iz-1)%Nplanes,it,cut)[ip] 
    Te = (2./3.)*(Te_perp + Ee_para) 
    V_th_i = np.sqrt((2.*joule_conv*Ti)/m_i)
    V_th_e = np.sqrt((2.*joule_conv*Te)/m_e)
    print('V_th_i:', '%.2E' % Decimal(V_th_i))
    print('V_th_e:', '%.2E' % Decimal(V_th_e))
    
    #calculate values at the cut and find the gyrofrequency
    br = getcutvalue_hor(bfield[:,0],cut,option)[ip]
    bz = getcutvalue_hor(bfield[:,1],cut,option)[ip]
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)[ip]
    B = np.sqrt(br**2+bz**2+bzeta**2)    
    print('B:','%.2E' % Decimal(B))
    Omega = ((charge*B)/m_i)
    print('Omega:','%.2E' % Decimal(Omega))
    print('rho:','%.2E' % Decimal(V_th_i/Omega))
    return V_th_i,Omega

#evaluates an approximation for the local ion banana width.
def bw(ip,iz,it,cut):
    q, eps, _ = q_loc(ip,cut)
    V_th, Omega = gf(ip,iz,it,cut)
    bw = (math.pi/np.sqrt(2)) *((V_th*math.fabs(q))/(np.sqrt(eps)*Omega))
    return bw

#calculates local values of gyrofrequency,thermal velocity and Larmor radius for electrons.
def gf_e(ip,iz,it,cut):
    global Zi,Ri,pot,bfield,B2
    option = 1

    m_e = 9.1e-31
    k = 1.38e-23
    joule_conv = 1.602e-19
    charge = 1.602176e-19 
    permeability = 1.2566e-6
    #c = 3e8 
    Te_perp = angle.map_array_rect(e_T_perp,(iz-1)%Nplanes,it,cut)[ip] 
    Ee_para = angle.map_array_rect(e_E_para,(iz-1)%Nplanes,it,cut)[ip] 
    Te = 2/3.*(Te_perp + Ee_para) 
    V_th = np.sqrt((joule_conv*Te)/m_e)
    #print('V_th:', '%.2E' % Decimal(V_th))
    
    #calculate values at the cut and find the gyrofrequency
    br = getcutvalue_hor(bfield[:,0],cut,option)[ip]
    bz = getcutvalue_hor(bfield[:,1],cut,option)[ip]
    bzeta = getcutvalue_hor(bfield[:,2],cut,option)[ip]
    B = np.sqrt(br**2+bz**2+bzeta**2)    
    #print('B:','%.2E' % Decimal(B))
    Omega = (charge*B/(m_e))
    #print('Omega:','%.2E' % Decimal(Omega))
    #print('rho:','%.2E' % Decimal(V_th/Omega))
    return V_th,Omega

#evaluates an approximation for the local electron banana width.
def bw_e(ip,iz,it,cut):
    q, eps, _ = q_loc(ip,cut)
    V_th, Omega = gf_e(ip,iz,it,cut)
    bw = (math.pi/np.sqrt(2)) *((V_th*math.fabs(q))/(np.sqrt(eps)*Omega))
    return bw

'''Weighted line integral of function.'''
def line_int(Rlist, weight, funct):
    R_start = 0
    R_end = len(Rlist)-1
    prefactor = (Rlist[R_end]-Rlist[R_start])/(len(Rlist[R_start:R_end])-1)
    temp = weight*funct
    Summation = (1./2.)*weight[R_start]*funct[R_start] + (1./2.)*weight[R_end]*funct[R_end] + np.sum(temp[R_start+1:R_end-2])
    Integral = prefactor*Summation
    return Integral

def mode_strength(iz,cut):
    option=1
    if cut==None:
        cut = jcut
        
    dR = Ri[1] - Ri[0]
    #density
    n_all = np.array([[getcutvalue_hor(ne[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(ne.shape[2])])
    n_tor_avg = n_all.mean(axis=1)
    n_t = n_tor_avg.mean(axis=0)
    #delta n^2
    dn = np.array(n_all[:,:] - n_t)
    dn2 = np.square(dn)
    dn2_tor = dn2.mean(axis=1)
    dn2_av = dn2_tor.mean(axis=0)
    return dn2_av

def plot_ms_vpol(iz,it,cut):
    dn2_av = mode_strength(iz,cut)
    v_pol = av_flows(cut,it)
    v_pol = v_pol.mean(axis=0)
    
    R = [ip*unitR + Rmin for ip in range(100)]
    #file = open("/global/cscratch1/sd/giannos/d3d_mode(fig.5.a).txt",'a')
    #file.write("psi"+"\t"+"v_pol"+"\t"+"mode"+"\n")
    #for i in range(len(psii[3:-3])):
    #    file.write(str(psii[3+i])+"\t"+str(v_pol[3+1])+"\t"+str(np.sqrt(dn2_av[3+1]))+"\n")
    #file.close()
    #plots
    fig, ax1 = plt.subplots()
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax1.plot(psii[3:-3],v_pol[3:-3],'b-')
    ax1.axvline(x=1.,color = 'k', linestyle = '--')
    ax1.axvline(x=1.01,color='y',linestyle = '-.')
    ax1.set_xlabel(r'$\Psi_N$')
    ax1.set_ylabel(r'$V_{\theta}\, \left(\frac{m}{s}\right)$')
    ax1.yaxis.set_major_locator(plt.MaxNLocator(10))
    ax1.xaxis.set_major_locator(plt.MaxNLocator(11))
    ax1.yaxis.set_minor_locator(plt.MaxNLocator(50))
    ax1.xaxis.set_minor_locator(plt.MultipleLocator(base=0.01))
    ax1.grid(which = 'major',linestyle='-')
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax2.plot(psii[3:-3],np.sqrt(dn2_av)[3:-3],'g-')
    ax2.set_ylabel(r'$(\delta n)_{rms}$')
    fig.tight_layout()
    plt.savefig('pol_flow.png')
    plt.show()

def plot_sh_vpol(arr,iz,it,cut):
    #sh = shear(it,cut)
    sh = arr[50,:]
    v_pol = av_flows(cut,it)
    v_pol = v_pol.mean(axis=0)
    
    R = [ip*unitR + Rmin for ip in range(100)]
    #file = open("/global/cscratch1/sd/giannos/d3d_mode(fig.5.a).txt",'a')
    #file.write("psi"+"\t"+"v_pol"+"\t"+"mode"+"\n")
    #for i in range(len(psii[3:-3])):
    #    file.write(str(psii[3+i])+"\t"+str(v_pol[3+1])+"\t"+str(np.sqrt(dn2_av[3+1]))+"\n")
    #file.close()
    #plots
    fig, ax1 = plt.subplots()
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    lns1 = ax1.plot(psii[3:-3],v_pol[3:-3],'b-',label=r'$V_{\theta}$')
    ax1.axvline(x=1.,color = 'k', linestyle = '--')
    ax1.set_xlabel(r'$\Psi_N$')
    ax1.set_ylabel(r'$V_{\theta}\, \left(\frac{m}{s}\right)$')
    #ax1.yaxis.set_major_locator(plt.MaxNLocator(10))
    #ax1.xaxis.set_major_locator(plt.MaxNLocator(8))
    #ax1.yaxis.set_minor_locator(plt.MaxNLocator(35))
    #ax1.xaxis.set_minor_locator(plt.MultipleLocator(base=0.01))
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    ax1.grid(which = 'major',linestyle='-')
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    lns2 = ax2.plot(psii[3:-3],sh[3:-3]/(2*np.pi*1000),'g-',label=r'$\Omega_{E\times B}$')
    ax2.set_ylabel(r'$\Omega_{E\times B} (kHZ)$')
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=1)
    fig.tight_layout()
    plt.savefig('d3d_shear_flow.png',bbox_inches='tight')
    plt.show()

def plot_sh_mode(arr,iz,it,cut):
    #sh = shear(it,cut)
    dn2_av = mode_strength(iz,cut)
    sh = arr[50,:]
    
    R = [ip*unitR + Rmin for ip in range(100)]
    #file = open("/global/u1/g/giannos/xgc_python_dir/cmod_mode_shear(fig.6.b).txt",'a')
    #file.write("psi"+"\t"+"shear"+"\t"+"mode"+"\n")
    #for i in range(len(psii[3:-3])):
    #    file.write(str(psii[3+i])+"\t"+str(sh[3+i])+"\t"+str(np.sqrt(dn2_av[3+i]))+"\n")
    #file.close()
    #plots
    fig, ax1 = plt.subplots()
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    lns1 = ax1.plot(psii[3:-3],np.sqrt(dn2_av)[3:-3],'b-',label=r'$(\delta n)_{rms}$')
    ax1.axvline(x=1.,color = 'k', linestyle = '--')
    ax1.set_xlabel(r'$\Psi_N$')
    ax1.set_ylabel(r'$(\delta n)_{rms}$')
    #ax1.yaxis.set_major_locator(plt.MaxNLocator(10))
    #ax1.xaxis.set_major_locator(plt.MaxNLocator(11))
    #ax1.yaxis.set_minor_locator(plt.MaxNLocator(50))
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax1.xaxis.set_minor_locator(plt.MultipleLocator(base=0.01))
    ax1.grid(which = 'major',linestyle='-')
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    lns2 = ax2.plot(psii[3:-3],sh[3:-3]/(2*np.pi*1000),'g-',label=r'$\Omega_{E\times B}$')
    ax2.set_ylabel(r'$\Omega_{E\times B} (kHZ)$')
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=1)
    fig.tight_layout()
    #plt.savefig('d3d_shear_mode.png',bbox_inches='tight')
    plt.show()
    
    
    
def profiles1d(cut):
    if cut==None:
        cut = int(math.floor(Te1d.shape[0]/2.))
    Te = Te1d[cut,:]
    Ti = Ti1d[cut,:]
    ne = ne1d[cut,:]
    fig, ax1 = plt.subplots()
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    lns1=ax1.plot(psin1d[62:],Ti[62:],'b-', label='$T_i$')
    lns2=ax1.plot(psin1d[62:],Te[62:],'r-', label='$T_e$')
    ax1.set_xlabel(r'$\Psi_N$')
    ax1.set_ylabel(r'$T(eV)$')
    #ax1.set_xlim([0.9058,1.0714])
    #ax1.axvline(x=1.,color = 'k', linestyle = '--')
    #ax1.axvline(x=1.01,color='y',linestyle = '-.')
    ax1.yaxis.set_major_locator(plt.MaxNLocator(10))
    ax1.xaxis.set_major_locator(plt.MaxNLocator(8))
    ax1.yaxis.set_minor_locator(plt.MaxNLocator(50))
    ax1.xaxis.set_minor_locator(plt.MultipleLocator(base=0.01))
    ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_formatter(mtick.ScalarFormatter())
    #ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    lns3=ax2.plot(psin1d[62:],ne[62:], 'm-', label='$n$')
    #ax2.set_xlim([0.9058,1.0714])
    ax2.set_ylabel(r'$n (m^{-3})$')
    fig.tight_layout()
    ax1.grid()
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=1)

    #ax1.legend()
    #ax2.legend()
    plt.show()
    plt.savefig("d3d_prof.png", bbox_inches='tight')
    
    
    
    
def profiles(iz,cut):
    option=1
    if cut==None:
        cut = jcut
        
    dR = Ri[1] - Ri[0]
    #density
    n_all = np.array([[getcutvalue_hor(ne[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(ne.shape[2])])
    n_tor_avg = n_all.mean(axis=1)
    n_t = n_tor_avg.mean(axis=0)
    #dn/dR
    #nplus = np.roll(n_t,-1)
    #nminus = np.roll(n_t,1)
    #dndR = (nplus-nminus)/(2.*dR)
    #dndR[0] = (n_t[1]-n_t[0])/dR
    #dndR[-1] = (n_t[-1]-n_t[-2])/dR
    #delta n^2
    #dn = np.array(n_all[:,:] - n_t)
    #dn2 = np.square(dn)
    #dn2_tor = dn2.mean(axis=1)
    #dn2_av = dn2_tor.mean(axis=0)
    #Ti
    Ti_perp = getcutvalue_hor(i_T_perp[:,iz,ne.shape[2]-1],cut,option) 
    Ei_para = getcutvalue_hor(i_E_para[:,iz,ne.shape[2]-1],cut,option) 
    Ti = (2./3.)*(Ti_perp + Ei_para)
    #dTi/dR
    #Tiplus = np.roll(Ti,-1)
    #Timinus = np.roll(Ti,1)
    #dTidR = (Tiplus-Timinus)/(2.*dR)
    #dTidR[0] = (Ti[1]-Ti[0])/dR
    #dTidR[-1] = (Ti[-1]-Ti[-2])/dR
    #dTe/dR
    Te_perp = getcutvalue_hor(e_T_perp[:,iz,ne.shape[2]-1],cut,option) 
    Ee_para = getcutvalue_hor(e_E_para[:,iz,ne.shape[2]-1],cut,option) 
    Te = (2./3.)*(Te_perp + Ee_para)
    #Teplus = np.roll(Te,-1)
    #Teminus = np.roll(Te,1)
    #dTedR = (Teplus-Teminus)/(2.*dR)
    #dTedR[0] = (Te[1]-Te[0])/dR
    #dTedR[-1] = (Te[-1]-Te[-2])/dR
    '''file = open("/global/cscratch1/sd/giannos/cmod_profiles(fig.1.d).txt",'a')
    file.write("psi"+"\t"+"n"+"\t"+"T_i"+"T_e"+"\n")
    for i in range(len(psii)):
        file.write(str(psii[i])+"\t"+str(n_t[i])+"\t"+str(Ti[i])+"\t"+str(Te[i])+"\n")
    file.close()'''

    
    R = [ip*unitR + Rmin for ip in range(100)]
    
    #Integrals
    #Ln_inv = line_int(list(R[50:72]), dn2_av[50:72], -dndR[50:72])/line_int(list(R[50:72]),dn2_av[50:72],n_t[50:72])
    #Lti_inv = line_int(list(R[50:72]), dn2_av[50:72], -dTidR[50:72])/line_int(list(R[50:72]),dn2_av[50:72],Ti[50:72])
    #Lte_inv = line_int(list(R[50:72]), dn2_av[50:72], -dTedR[50:72])/line_int(list(R[50:72]),dn2_av[50:72],Te[50:72])
    
    #plots
    fig, ax1 = plt.subplots()
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    lns1=ax1.plot(psii[3:-3],Ti[3:-3],'b-', label='$T_i$')
    lns2=ax1.plot(psii[3:-3],Te[3:-3],'r-', label='$T_e$')
    ax1.set_xlabel(r'$\Psi_N$')
    ax1.set_ylabel(r'$T(eV)$')
    ax1.set_xlim([0.9058,1.0714])
    #ax1.axvline(x=1.,color = 'k', linestyle = '--')
    #ax1.axvline(x=1.01,color='y',linestyle = '-.')
    ax1.yaxis.set_major_locator(plt.MaxNLocator(10))
    ax1.xaxis.set_major_locator(plt.MaxNLocator(11))
    ax1.yaxis.set_minor_locator(plt.MaxNLocator(50))
    ax1.xaxis.set_minor_locator(plt.MultipleLocator(base=0.01))
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_formatter(mtick.ScalarFormatter())
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    #ax2.plot(psii[3:-3],np.sqrt(dn2_av[3:-3]),'g-')
    lns3=ax2.plot(psii[3:-3],n_t[3:-3], 'm-', label='$n$')
    #ax2.set_xlim([0.9058,1.0714])
    #ax2.set_ylabel(r'$(\delta n)_{rms}\,, n\times 10^{-1}$')
    ax2.set_ylabel(r'$n (m^{-3})$')
    fig.tight_layout()
    ax1.grid()
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=1)

    #ax1.legend()
    #ax2.legend()
    #plt.show()
    plt.savefig("cmod_prof.png", bbox_inches='tight')
    return 0
    #return Ln_inv, Lti_inv, Lte_inv 
    
#R-derivative of density on a cut. Returns also the scale length.
def dndR(iz,it,cut):
    option=1
    if cut==None:
        cut = jcut
       
    dR = Ri[1] - Ri[0]
    n = getcutvalue_hor(ne[:,iz,it],cut,option)
    nplus = np.roll(n,-1)
    nminus = np.roll(n,1)
    dndR = (nplus-nminus)/(2.*dR)
    dndR[0] = (n[1]-n[0])/dR
    dndR[-1] = (n[-1]-n[-2])/dR
    Ln = n/np.absolute(dndR)
    return dndR, Ln

def Lninv(cut):
    option=1
    if cut==None:
        cut = jcut
       
    dR = Ri[1] - Ri[0]
    n_all = np.array([[getcutvalue_hor(ne[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(2)])
    n_tavg = n_all.mean(axis=0)
    n = n_tavg.mean(axis=0)
    n0 = 3e20 
    nplus = np.roll(n,-1)
    nminus = np.roll(n,1)
    dndR = (nplus-nminus)/(2.*dR)
    dndR[0] = (n[1]-n[0])/dR
    dndR[-1] = (n[-1]-n[-2])/dR
    Lninv = -dndR/n
    Ln = n0/dndR
    R = [ip*unitR + Rmin for ip in range(100)]
    #fig = plt.figure()
    #plt.xlabel("$R(m)$")
    #plt.ylabel("$L^{-1}_n(m^{-1})$")
    #plt.ylabel('$n$')
    #plt.plot(R[0:100],Lninv[0:100],'b')
    #plt.plot(R[:],n[:],'r')
    #plt.show()
    return Lninv


#R-derivative of Ti on a cut. Returns also the scale length.
def dTidR(iz,it,cut):
    option=1
    if cut==None:
        cut = jcut
       
    dR = Ri[1] - Ri[0]
    Ti_perp = angle.map_array_rect(i_T_perp,(iz-1)%Nplanes,it,cut) 
    Ei_para = angle.map_array_rect(i_E_para,(iz-1)%Nplanes,it,cut) 
    Ti = 2/3.*(Ti_perp + Ei_para) 
    Tiplus = np.roll(Ti,-1)
    Timinus = np.roll(Ti,1)
    dTidR = (Tiplus-Timinus)/(2.*dR)
    dTidR[0] = (Ti[1]-Ti[0])/dR
    dTidR[-1] = (Ti[-1]-Ti[-2])/dR
    LTi = Ti/np.absolute(dTidR)
    return dTidR, LTi


#returns toroidal and time averaged inverse LTi.
def LTiinv(cut):
    option=1
    if cut==None:
        cut = jcut
       
    dR = Ri[1] - Ri[0]
    Ti_perp_all = np.array([[getcutvalue_hor(i_T_perp[:,iz,it],cut,option) for iz in range(Nplanes)]for it in range(2)]) 
    Ei_para_all = np.array([[getcutvalue_hor(i_E_para[:,iz,it],cut,option) for iz in range(Nplanes)]for it in range(2)]) 
    Ti_perp_tavg = Ti_perp_all.mean(axis=0)
    Ei_para_tavg = Ei_para_all.mean(axis=0)
    Ti_perp = Ti_perp_tavg.mean(axis=0)
    Ei_para = Ei_para_tavg.mean(axis=0)
    Ti = 2/3.*(Ti_perp + Ei_para) 
    Ti0 = 500
    Tiplus = np.roll(Ti,-1)
    Timinus = np.roll(Ti,1)
    dTidR = (Tiplus-Timinus)/(2.*dR)
    dTidR[0] = (Ti[1]-Ti[0])/dR
    dTidR[-1] = (Ti[-1]-Ti[-2])/dR
    LTiinv = -dTidR/Ti
    R = [ip*unitR + Rmin for ip in range(100)]
    #fig = plt.figure()
    #plt.xlabel("$R(m)$")
    #plt.ylabel("$L^{-1}_{Ti}(m^{-1})$")
    #plt.ylabel("$T_i$")
    #plt.plot(R[0:100],LTiinv[0:100],'b')
    #plt.show()
    return LTiinv



#R-derivative of Te on a cut. Returns also the scale length.
def dTedR(iz,it,cut):
    option=1
    if cut==None:
        cut = jcut
       
    dR = Ri[1] - Ri[0]
    Te_perp = angle.map_array_rect(e_T_perp,(iz-1)%Nplanes,it,cut) 
    Ee_para = angle.map_array_rect(e_E_para,(iz-1)%Nplanes,it,cut) 
    Te = 2/3.*(Te_perp + Ee_para) 
    Teplus = np.roll(Te,-1)
    Teminus = np.roll(Te,1)
    dTedR = (Teplus-Teminus)/(2.*dR)
    dTedR[0] = (Te[1]-Te[0])/dR
    dTedR[-1] = (Te[-1]-Te[-2])/dR
    LTe = Te/np.absolute(dTedR)
    return dTedR, LTe

#returns toroidal and time averaged inverse LTe.
def LTeinv(cut):
    option=1
    if cut==None:
        cut = jcut
       
    dR = Ri[1] - Ri[0]
    Te_perp_all = np.array([[getcutvalue_hor(e_T_perp[:,iz,it],cut,option) for iz in range(Nplanes)]for it in range(2)]) 
    Ee_para_all = np.array([[getcutvalue_hor(e_E_para[:,iz,it],cut,option) for iz in range(Nplanes)]for it in range(2)]) 
    Te_perp_tavg = Te_perp_all.mean(axis=0)
    Ee_para_tavg = Ee_para_all.mean(axis=0)
    Te_perp = Te_perp_tavg.mean(axis=0)
    Ee_para = Ee_para_tavg.mean(axis=0)
    Te = 2/3.*(Te_perp + Ee_para)
    Te0 = 400
    Teplus = np.roll(Te,-1)
    Teminus = np.roll(Te,1)
    dTedR = (Teplus-Teminus)/(2.*dR)
    dTedR[0] = (Te[1]-Te[0])/dR
    dTedR[-1] = (Te[-1]-Te[-2])/dR
    LTeinv = -dTedR/Te   
    R = [ip*unitR + Rmin for ip in range(100)]
    #fig = plt.figure()
    #plt.xlabel("$R(m)$")
    #plt.ylabel("$L^{-1}_{Te}(m^{-1})$")
    #plt.ylabel("$T_e$")
    #plt.plot(R[:],LTeinv[:],'b')
    #plt.show()
    return LTeinv



#calculates the local banana current.
def banana_current(ip,iz,it,cut):
    joule_conv = 1.602e-19
    charge = 1.602176e-19 
    q, eps, _ = q_loc(ip,cut)
    V_th, Omega = gf(ip,iz,it,cut)
    dndr, _ = dndR(iz,it,cut)
    dTdr, _ = dTidR(iz,it,cut)
    n = getcutvalue_hor(ne[:,iz,it],cut,1)
    bc = charge*(math.pi/np.sqrt(2))*(math.fabs(q)/Omega)*(math.pow(eps,1/2)*(n[ip]/2)*joule_conv*dTdr[ip] + math.pow(eps,1/2)*math.pow(V_th,2) * dndr[ip])
    bc_m = charge*(math.pi/np.sqrt(2))*(math.fabs(q)/Omega)*(math.pow(eps,3/2)*math.pow(V_th,2) * dndr[ip])
    return -bc,-bc_m, n

def Wesson_bootstrap(iz,it,cut):
    joule_conv = 1.602e-19
    temp = [q_loc(ip,cut) for ip in range(100)]
    eps = [x[1] for x in temp]
    dndr, _ = dndR(iz,it,cut)
    dTidr, _ = dTidR(iz,it,cut)
    dTedr, _ = dTedR(iz,it,cut)
    n = getcutvalue_hor(ne[:,iz,it],cut,1)
    Te_perp = angle.map_array_rect(e_T_perp,(iz-1)%Nplanes,it,cut)  
    Ee_para = angle.map_array_rect(e_E_para,(iz-1)%Nplanes,it,cut) 
    Te = 2/3.*(Te_perp + Ee_para)
    Ti_perp = angle.map_array_rect(i_T_perp,(iz-1)%Nplanes,it,cut)  
    Ei_para = angle.map_array_rect(i_E_para,(iz-1)%Nplanes,it,cut) 
    Ti = 2/3.*(Ti_perp + Ei_para)
    br = getcutvalue_hor(bfield[:,0],cut,1)
    bz = getcutvalue_hor(bfield[:,1],cut,1)
    bzeta = getcutvalue_hor(bfield[:,2],cut,1)
    Bpol = np.sqrt(br**2 + bz**2)
    Jb = -(np.sqrt(eps)/Bpol)*(2.44*joule_conv*(Te+Ti)*dndr+0.69*n*joule_conv*dTedr-0.42*n*joule_conv*dTidr)
    Jbb = -(np.sqrt(eps)/Bpol)*joule_conv*(Te+Ti)*dndr
    return Jbb

#compares the parallel flow from the data from the one that the banana current implies.
def comparison(ip,iz,it,cut):
    charge = 1.602176e-19
    V_para = angle.map_array_rect(i_u_para,(iz-1)%Nplanes,it,cut)
    n = getcutvalue_hor(ne[:,iz,it],cut,1)
    bwdth = bw(ip,iz,it,cut)
    dr = math.ceil(bwdth/unitR)
    bc,bc_m,_ = banana_current(ip-dr,iz,it,cut)
    Jbb = Wesson_bootstrap(iz,it,cut)[ip]
    #V_boot = Jb[ip]/(charge*n[ip])
    V_bboot = Jbb/(charge*n[ip])
    V_bn = bc/(charge*n[ip])
    comp = (V_para[ip]-math.fabs(V_bn))/V_para[ip]
    #comp_boot = (V_para[ip]-math.fabs(V_boot))/V_para[ip]
    #print("V_par_data = ", V_para[ip])
    #print("V_par_bn = ", V_bn)
    #print("relative error = ", comp)
    return comp, V_bn, V_para[ip],V_bboot

#Plots the comparison between banana current flow and data parallel flow.
def flow_comp():
    iz=0
    it=0
    cut=50
    vbn = []
    vpar = []
    vbboot = []
    err = []
    R = []
    for ip in range(14):
        error, V_bn, V_par,V_bboot = comparison(30+ip,iz,it,cut)
        err.append(math.fabs(error))
        vbn.append(V_bn)
        vpar.append(V_par)
        vbboot.append(V_bboot)
        R.append(Rmin+(30+ip)*unitR)
    
    fig, ax1 = plt.subplots()
    plt.xlabel("$R(m)$")
    #plt.style.use('seaborn-paper')
    h = plt.ylabel(r"$V_{\parallel} \,\, \left(\frac{m}{s}\right)$")
    h.set_rotation(0)
    ax1.yaxis.set_label_coords(-0.1, 0.5)
    ax1.plot(R,vpar,'b',marker='o',label='data flow')
    ax1.plot(R,vbn,'y',marker='o',label='banana flow')
    #ax1.plot(R,vboot,'m',marker='o',label='bootstrap flow')
    ax1.plot(R,vbboot,'m',marker='o',label='bootstrap ban')
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #ax2 = ax1.twinx()
    #ax2.plot(R,err,'-.r',marker='o',label=r'$\%$ error')
    #ax2.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: '{:1.0%}'.format(y))) 
    ax1.axvline(x=2.245,color='k',linestyle = '--')
    #fig.tight_layout()
    ax1.grid()
    ax1.legend(loc=2)
    #ax2.legend(loc=1)
    plt.show()
    
#checks for the difference between banana current flow and parallel flow
def rel_error():
    err = []
    R = []
    iz=0
    it=0
    cut=50
    V_para = angle.map_array_rect(i_u_para,(iz-1)%Nplanes,it,cut)
    psi = getcutvalue_hor(psin[:],50,1)
    sep = np.where(psi>=0.999)[0][0]
    for ip in range(20):
        err.append(math.fabs(comparison(30+ip,0,0,50)[0]))
        R.append(Rmin+(30+ip)*unitR)
        
    fig=plt.figure()
    plt.xlabel("R")
    plt.ylabel(r"relative error")
    plt.plot(R,err)
    plt.ylim((0,2))
    plt.axvline(x=2.245,color='k')
    plt.legend()
    plt.show()

#Plots inverse scale lengths        
def slp(iz,it,cut):
    _, Lti = dTidR(iz,it,cut)
    _, Lte = dTedR(iz,it,cut)
    _, Ln = dndR(iz,it,cut)
    R = [ip*unitR + Rmin for ip in range(100)]
    fig = plt.figure()
    plt.xlabel("R(m)")
    plt.ylabel("L(m)")
    plt.plot(R[0:40],Lti[0:40],'b')
    plt.plot(R[0:40],Lte[0:40],'r')
    plt.plot(R[0:40],Ln[0:40],'g')
    plt.show()

    
#Plots fluxes wrt eta_e for ETG consideration    
def etg():
    eta_file = open("/global/homes/g/giannos/xgc_python_dir/etas.txt","r")
    next(eta_file)
    angle = []
    etae = []
    etai = []
    
    
    for columns in ( raw.strip().split() for raw in eta_file):  
        angle.append(float(columns[0]))
        etae.append(float(columns[1]))
        etai.append(float(columns[2]))  
            
    
    particle_file = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/sep_flux/shorter_run/ti255_shorter_run_total.txt","r")
    next(particle_file)
    par_tur = []
    angle_par = []
    for columns in ( raw.strip().split() for raw in particle_file):
        angle_par.append(float(columns[0]))
        par_tur.append(float(columns[2]))
        
    heat_file = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/heat_flux/ion/Ion_heat_fluxes_total.txt","r")
    next(heat_file)
    heat_tur = []
    angle_heat = []
    for columns in ( raw.strip().split() for raw in heat_file):
        angle_heat.append(float(columns[0]))
        heat_tur.append(float(columns[2])/200)
        
    plt.style.use('ggplot')
    fig, ax1 = plt.subplots()
    plt.title("Turbulent Fluxes")
    plt.xlabel("angle")
    plt.ylabel("Fluxes")                 
    #plt.ylabel(r"$\eta$")
    ax1.set_xlim([0.0,360.0])
    #plt.ylim([0.0,3.0])
    ax1.plot(angle_par[:], par_tur[:],'b',marker='o',label='Particle')
    ax1.plot(angle_heat[:], heat_tur[:],'r',marker='o',label=r'Heat') 
    ax2 = ax1.twinx()
    ax2.set_xlim([0.0,360.0])
    ax2.set_ylim([0.0,3.0])
    ax2.plot(angle[:], etae[:],'g',marker='o',label=r'$\eta_e$')
    #ax2.plot(angle[:], etai[:],'g',marker='o',label=r'$\eta_i$')
    #ax2.axhline(y=2./3.,color='k')
    #ax2.axhline(y=3./2.,color='k')
    fig.tight_layout()
    ax1.grid()
    ax1.legend(loc=2)
    ax2.legend(loc=1)
    plt.show()

#Plots shear vs inverse scale lengths    
def ssl():
    shear_file = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/shear/Shear_total.txt","r")
    next(shear_file)
    angle = []
    shear = [] 
    
    
    for columns in ( raw.strip().split() for raw in shear_file):  
        angle.append(float(columns[0]))
        shear.append(float(columns[1])/(1*np.pi))
        
    sl_file = open("/global/homes/g/giannos/xgc_python_dir/Scale_lengths_2.txt")
    next(sl_file)
    angle_l = []
    lte = []
    lti = []
    ln = []
    for columns in ( raw.strip().split() for raw in sl_file):
        angle_l.append(float(columns[0]))
        lte.append(float(columns[1]))
        lti.append(float(columns[2]))
        ln.append(float(columns[3]))
    
    lte = np.asarray(lte)
    lti = np.asarray(lti)
    ln = np.asarray(ln)
    
    lpi = ln + lti
    lpe = ln + lte
    
    
    heat_file = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/heat_flux/electron/Electron_heat_fluxes_total.txt","r")
    next(heat_file)
    heat_tur = []
    angle_heat = []
    for columns in ( raw.strip().split() for raw in heat_file):
        angle_heat.append(float(columns[0]))
        heat_tur.append(float(columns[2]))
    
    heat_e = np.asarray(heat_tur)
    
    particle_file = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/sep_flux/shorter_run/ti255_shorter_run_total.txt","r")
    next(particle_file)
    par_tur = []
    angle_par = []
    for columns in ( raw.strip().split() for raw in particle_file):
        angle_par.append(float(columns[0]))
        par_tur.append(float(columns[2]))
    
    par = np.asarray(par_tur)
    
    heat_file_i = open("/global/homes/g/giannos/xgc_python_dir/ti255_analysis/heat_flux/ion/Ion_heat_fluxes_total.txt","r")
    next(heat_file_i)
    heat_tur_i = []
    angle_heat_i = []
    for columns in ( raw.strip().split() for raw in heat_file_i):
        angle_heat_i.append(float(columns[0]))
        heat_tur_i.append(float(columns[2]))
        
    heat_i = np.asarray(heat_tur_i)
    
    temp_file = open("/global/homes/g/giannos/xgc_python_dir/Sep_temperatures.txt","r")
    next(temp_file)
    ti = []
    te = []
    for columns in ( raw.strip().split() for raw in temp_file):
        ti.append(float(columns[1]))
        te.append(float(columns[2]))
    
    T_e = np.asarray(te)
    T_i = np.asarray(ti)
    
    Q_i_ov_G = (lpi*T_i)/ln
    Q_e_ov_G = (lpe*T_e)/ln
        
    plt.style.use('ggplot')
    fig, ax1 = plt.subplots()
    plt.xlabel("angle")
    plt.ylabel(r"$\frac{Q_i}{\Gamma}$")                 
    #plt.ylabel(r"$\eta$")
    ax1.set_xlim([0.0,360.0])
    plt.ylim([-500.0,1000.0])
    #ax1.plot(angle_l[:], Q_i_ov_G[:],'g',marker='o',label=r'convective')
    ax1.plot(angle_l[:], heat_i[:]/par[:],'b',marker='o',label=r'actual') 
    #ax1.plot(angle_l[:], ln[:],'b',marker ='o',label=r'$L_n$')
    #ax2 = ax1.twinx()
    #ax2.set_xlim([0.0,360.0])
    #ax2.set_ylim([0.0,3.0])
    #ax2.plot(angle[:], shear[:],'k',marker='o',label=r'$\Omega_{E\times B}\,(Hz)$')
    #ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0) ,useMathText=True)
    #ax2.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: '{:%.2e}'.format(y)))
    #ax2.plot(angle[:], etai[:],'g',marker='o',label=r'$\eta_i$')
    #ax2.axhline(y=2./3.,color='k')
    #ax2.axhline(y=3./2.,color='k')
    #fig.tight_layout()
    ax1.grid()
    ax1.legend(loc=3)
    #ax2.legend(loc=1)
    plt.show()
    
    
    
# finds average values of velocities at the midplane for feeding into the p.d.f    
def vel_distr(Flux_ip,verbose=True):
    #for the particular set of input values the point of the separatrix at jcut is between 40 and 41.
    cut = jcut
    option = 1
    m_i = 1.6723e-27
    k = 1.38e-23
    joule_conv = 1.602e-19
    start = 0
    time_range = 200
    #data loading
    temp = np.array([[getcutvRvZ(iz,it,cut) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    VR_all = temp[:,:,0,Flux_ip]
    VZ_all = temp[:,:,1,Flux_ip]
    Ti_perp_all = np.array([[getcutvalue_hor(i_T_perp[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ei_para_all = np.array([[getcutvalue_hor(i_E_para[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])
    V_para_all = np.array([[getcutvalue_hor(i_u_para[:,iz,it],cut,option)[Flux_ip] for iz in range(Nplanes)] for it in range(start,start+time_range)])   
    gradBR,gradBZ,_,curvR,curvZ,_ = magnetic_drifts(cut,Flux_ip,start)
    Temperature = 2/3.*(Ti_perp_all + Ei_para_all) 
    
    #toroidal averages
    VR = VR_all.mean(axis=1)
    VZ = VZ_all.mean(axis=1)
    Temp_tor = Temperature.mean(axis=1)
    V_par = V_para_all.mean(axis=1)
    V_th = np.sqrt((joule_conv*Temp_tor)/(2*m_i))
    
    #time averages
    VR_avg = VR.mean()
    VZ_avg = VZ.mean()
    V_par_avg = V_par.mean()
    V_th_avg = V_th.mean()
    
    #V_perp
    V_perp = np.sqrt(np.power(VR_avg+gradBR+curvR,2)+np.power(VZ_avg+gradBZ+curvZ,2))
    V_exb = np.sqrt(np.power(VR_avg,2)+np.power(VZ_avg,2))
    V_mag = np.sqrt(np.power(gradBR+curvR,2) + np.power(gradBZ+curvZ,2))
    
    
    
    print("V_th",V_th_avg)
    print("V_par",V_par_avg)
    print("V_perp", V_perp)
    print("V_exb",V_exb)
    print("V_mag", V_mag)
    if verbose == True:
        time = [x for x in range(time_range)]
        plt.style.use('ggplot')
        plt.figure()
        plt.title("Velocity distribution")
        plt.xlabel("time step")
        plt.ylabel("velocity")                 
        plt.plot(time[:], V_par[:],'b',marker='o',label='V_par')
        plt.plot(time[:], V_th[:],'r',marker='o',label='V_th')
        plt.plot(time[:], VR[:],'g',marker='o',label='VR')
        plt.plot(time[:], VR[:],'y',marker='o',label='VZ')
        plt.grid()
        plt.legend(loc=2)
        plt.show()
    else:
        pass

def frequencies():
    
    #for the particular set of input values the point of the separatrix at jcut is between 40 and 41.
    cut = jcut
    option = 1
    Z = 1
    m_i = 1.6723e-27
    m_e = 9.109e-31
    k = 1.38e-23
    joule_conv = 1.602e-19
    start = 0
    time_range = 2
    #data loading
    #Temperatures are in ev.
    n_e=loader.calcNeTotal()
    ne_all = np.array([[getcutvalue_hor(n_e[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ti_perp_all = np.array([[getcutvalue_hor(i_T_perp[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ei_para_all = np.array([[getcutvalue_hor(i_E_para[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Temperature_i = 2/3.*(Ti_perp_all + Ei_para_all) 
    V_para_all = np.array([[getcutvalue_hor(i_u_para[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(start,start+time_range)]) 
    Te_perp_all = np.array([[getcutvalue_hor(e_T_perp[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Ee_para_all = np.array([[getcutvalue_hor(i_E_para[:,iz,it],cut,option) for iz in range(Nplanes)] for it in range(start,start+time_range)])
    Temperature_e = 2/3.*(Te_perp_all + Ee_para_all) 
    #Averaging
    n_e_tor = ne_all.mean(axis=1)
    Temp_i_tor = Temperature_i.mean(axis=1)
    Temp_e_tor = Temperature_e.mean(axis=1)
    V_para_tor = V_para_all.mean(axis=1)
    ne = n_e_tor.mean(axis=0)
    Ti = Temp_i_tor.mean(axis=0)
    Te = Temp_e_tor.mean(axis=0)
    V_par = V_para_tor.mean(axis=0)
    V_th = np.sqrt((joule_conv*Ti)/(2*m_i))
    V_th_e = np.sqrt((joule_conv*Te)/(m_e))
    R = [Rmin + ip*unitR for ip in range(100)]
    q = np.array([q_loc(ip,jcut)[0] for ip in range(100)])
    eps = np.array([q_loc(ip,jcut)[1] for ip in range(100)])
    C_s = np.sqrt(((5./3.)*joule_conv*Te)/(2*m_i))
    
    #Coulomb logarithms
    LogLi = 17.3 - (1/2)*np.log(ne/10e20) + 3/2*np.log(Ti/1000) 
    LogLe = 15.2 - (1/2)*np.log(ne/10e20) + 3/2*np.log(Te/1000) 
    
    
    #Collision times
    tau_e = 0.343e12*(np.power(Te,3./2.)/(np.power(Z,2)*ne*LogLe))
    tau_i = 0.208e14*((np.power(Ti,3./2.)*(np.sqrt(2.)))/(np.power(Z,4)*ne*LogLi))
    #Bounce time
    t_b = (2*np.pi*np.absolute(q)*R)/(np.sqrt(eps)*V_th)
    t_eff = (2*np.pi*np.absolute(q)*R)/(np.power(eps,3./2.)*V_th)
    t_transit = (np.absolute(q)*R)/V_th
    t_transit_e = (np.absolute(q)*R)/V_th_e
    
    #check
    #t_e = 1.09e16*(np.power(1,3./2.)/(10e19*17))
    #t_i = 6.6e17*((np.power(1,3./2.)*(np.sqrt(2.)))/(10e19*17))
    #ten = 0.343e12*(np.power(1000,3./2.)/(10e19*17))
    #tin = 0.208e14*((np.power(1000,3./2.)*(np.sqrt(2.)))/(10e19*17))
    #print(t_i,tin,t_e,ten)
    print("ti","te","nu_i/omega_t","nu_e/omega_te")
    print(tau_i[55],tau_e[55],(1/tau_i[55])/(1/t_transit[55]),(1/tau_e[55])/(1/t_transit_e[55]),np.power(eps[55],3./2.))
    print(tau_i[56],tau_e[56],(1/tau_i[56])/(1/t_transit[56]),(1/tau_e[56])/(1/t_transit_e[56]),np.power(eps[56],3./2.))
    print(tau_i[57],tau_e[57],(1/tau_i[57])/(1/t_transit[57]),(1/tau_e[57])/(1/t_transit_e[57]),np.power(eps[57],3./2.))
    print(tau_i[58],tau_e[58],(1/tau_i[58])/(1/t_transit[58]),(1/tau_e[58])/(1/t_transit_e[58]),np.power(eps[58],3./2.))
    print("C_s",C_s[55],C_s[56],C_s[57],C_s[58])
    print("V_th_i", V_th[55],V_th[56],V_th[57],V_th[58])
    print("V_th_e", V_th_e[55],V_th_e[56],V_th_e[57],V_th_e[58])
    
    #x = [x for x in range(100)]
    #plt.plot(x,tau_e,'b-')
    #plt.plot(x,tau_i,'r-')
    #plt.grid()
    #plt.show()
    #toroidal averages
    #return ne,Ti,Te
    
#calculates the difference in temperature accross a banana width.      
def ban_P_diff():    
#banana width calculated to be 6mm and under the loading values, unitR = 1mm. Separatrix at point 39-40.
    option = 1
    Ti_perp_L = getcutvalue_hor(i_T_perp[:,0,0],jcut,option)[37]
    Ei_para_L = getcutvalue_hor(i_E_para[:,0,0],jcut,option)[37]
    Temperature_i_L = 2/3.*(Ti_perp_L + Ei_para_L) 
    Ti_perp_R = getcutvalue_hor(i_T_perp[:,0,0],jcut,option)[43]
    Ei_para_R = getcutvalue_hor(i_E_para[:,0,0],jcut,option)[43]
    Temperature_i_R = 2/3.*(Ti_perp_R + Ei_para_R) 
    
    deltaTi = Temperature_i_L - Temperature_i_R
    Ti_perp = getcutvalue_hor(i_T_perp[:,0,1],jcut,option)
    Ei_para = getcutvalue_hor(i_E_para[:,0,1],jcut,option)
    Temperature_i = 2/3.*(Ti_perp + Ei_para) 

    Te_perp = getcutvalue_hor(e_T_perp[:,0,1],jcut,option)
    Ee_para = getcutvalue_hor(e_E_para[:,0,1],jcut,option)
    Temperature_e = 2/3.*(Te_perp + Ee_para) 

    Te_perp_L = getcutvalue_hor(e_T_perp[:,0,0],jcut,option)[37]
    Ee_para_L = getcutvalue_hor(e_E_para[:,0,0],jcut,option)[37]
    Temperature_e_L = 2/3.*(Te_perp_L + Ee_para_L) 
    Te_perp_R = getcutvalue_hor(e_T_perp[:,0,0],jcut,option)[43]
    Ee_para_R = getcutvalue_hor(e_E_para[:,0,0],jcut,option)[43]
    Temperature_e_R = 2/3.*(Te_perp_R + Ee_para_R) 
    
    deltaTe = Temperature_e_L - Temperature_e_R
    print("Ti_43",Temperature_i_R)
    print("Ti_37",Temperature_i_L)
    print("Te_43",Temperature_e_R)
    print("Te_37",Temperature_e_L)
    print("DTi", deltaTi)
    print("DTe", deltaTe)
    x=[x for x in range(100)]
    plt.plot(x,Temperature_i,'b')
    plt.plot(x,Temperature_e,'r')
    plt.show()
    
def file_del():
    '''Test function for file deletion.'''
    folder = '/global/homes/g/giannos/xgc_python_dir/phis_short'
    the_file = 'ti255_phis_vs.angle_0.txt'
    file_path = os.path.join(folder, the_file)
    if os.path.isfile(file_path):
            os.unlink(file_path)
    else:
        print("File Not Found")
    
def line_read():
    '''Test function. Reads lines and discards nan's and '--' string.'''
    f = open("ti255_particle_fluxes_vs.angle_0.txt","r")
    f.seek(0)
    lines = f.readlines()
    num_of_lines = len(lines)
    f.close()
    f = open("ti255_particle_fluxes_vs.angle_0.txt","w")
    for i in range(num_of_lines):
        if "--" and "nan" not in lines[i]:
            f.write(lines[i])            
    f.close()
    
'''Functions that perform theta vs. psi contour plots. Currently used for equilibrium potential.'''

def contour_phi_plot(x,y,z,t):
    '''Does the contour plot of z.'''
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.figure()
    #lev = np.arange(-160,-20,10)
    CS = plt.contour(x,y,z,40,colors='k',linewidths=0.7)
    plt.clabel(CS, CS.levels[-20:], inline=True, fmt = '%.1f', fontsize=8)
    plt.plot(1.0,-1.8,'x')
    #plt.colorbar(CS)
    plt.xlabel(r'$\Psi_N$')
    plt.ylabel(r'$\theta$ (rad)')
    plt.title(r'$\Phi$')
    plt.savefig("cmod_isopot.png")
    plt.close()

def contour_phi_plot2(x,y,z,t,s):
    '''Does the contour plot of z.'''
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.figure()
    #lev = np.arange(-160,-20,10)
    CS = plt.contour(x,y,z,40,colors='k',linewidths=0.7)
    plt.clabel(CS, CS.levels[-20:], inline=True, fmt = '%.1f', fontsize=8)
    plt.plot(1.0,-1.8,'x')
    #plt.colorbar(CS)
    plt.xlabel(r'$\Psi_N$')
    plt.ylabel(r'$\theta$ (rad)')
    plt.title(r'%s' %(s))
    plt.savefig("%s_isopot.png" %(s))
    plt.close()

def contour_phi_plot3(x,y1,y2,t):
    plt.figure()
    plt.plot(x,y1,'r-')
    plt.plot(x[10],y1[10],'ro',label = r'$\frac{\phi}{T_e}$')
    plt.plot(x,y2,'b-')
    plt.plot(x[10],y2[10],'bo',label = r'$log\left(\frac{n_e}{n_o}\right)$')
    plt.xlabel(r'$\theta$ (rad)')
    plt.legend()
    plt.title(r'M-B $\psi = 1.0$')
    plt.savefig("MB3.png")
    plt.close()
    
    
def prepare_data2(t):
    no = 6.5e20
    Te_perp = e_T_perp[:,:,t]
    Te_para = e_E_para[:,:,t]
    Te = 2/3.*(Te_perp + Te_para)
    Te = Te.mean(axis=1)
    density = ne[:,:,t]
    density = density.mean(axis=1)
    potential = pot[:,:,t]
    potential = potential.mean(axis=1)
    psi = psin[:]
    thetas = np.arctan2(RZ[:,1],RZ[:,0]-Rmaj)
    keep = keep_nodes(psi)
    potential = keep_values(potential,keep)
    psi = keep_values(psi,keep)
    thetas = keep_values(thetas,keep)
    density = keep_values(density,keep)
    Te = keep_values(Te,keep)
    npot = [i/j for i,j in zip(potential,Te)]
    logn = [np.log(i/no) for i in density]
    return psi, thetas, npot , logn
    
    
def prepare_data(t):
    '''Loads the node data of the relevant quantities and keeps only the nodes within a psi range.'''
    potential = pot[:,:,t]
    potential = potential.mean(axis=1)
    psi = psin[:]
    thetas = np.arctan2(RZ[:,1],RZ[:,0]-Rmaj)
    keep = keep_nodes(psi)
    potential = keep_values(potential,keep)
    psi = keep_values(psi,keep)
    thetas = keep_values(thetas,keep)
    return psi, thetas, potential

def write_contour_data():
    psi, thetas, potential = prepare_data(2)
    file = open("/global/cscratch1/sd/giannos/cmod_phi_contour(fig.3.d).txt",'a')
    file.write("psi"+"\t"+"thetas"+"\t"+"potential"+"\n")
    for i in range(len(psi)):
        file.write(str(psi[i])+"\t"+str(thetas[i])+"\t"+str(potential[i])+"\n")
    file.close()

def triang_and_smooth(t):
    '''Does a triangulation on the psi-theta space, performs an interpolation there, a smoothing operation and then outputs the plot.'''
    from scipy.ndimage import gaussian_filter
    t_range = pot.shape[2]
    #for t in range(1,t_range):
    psi, thetas, potential = prepare_data(t)
    Pi = np.linspace(min(psi),max(psi),100)
    Ti = np.linspace(min(thetas),max(thetas),100)
    (PI,TI) = np.meshgrid(Pi,Ti)
    triangObject = Triangulation(psi,thetas)
    tcp = LinearTriInterpolator(triangObject,potential)
    pot_int = tcp(PI,TI)
    filtered = gaussian_filter(pot_int,1.8)
    contour_phi_plot(PI,TI,filtered,t)
    print("t is done!", t)

def triang_and_smooth2(t):
    '''Does a triangulation on the psi-theta space, performs an interpolation there, a smoothing operation and then outputs the plot.'''
    from scipy.ndimage import gaussian_filter
    t_range = pot.shape[2]
    #for t in range(1,t_range):
    psi, thetas, potential, density = prepare_data2(t)
    Pi = np.linspace(min(psi),max(psi),100)
    Ti = np.linspace(min(thetas),max(thetas),100)
    (PI,TI) = np.meshgrid(Pi,Ti)
    triangObject = Triangulation(psi,thetas)
    tcp1 = LinearTriInterpolator(triangObject,potential)
    tcp2 = LinearTriInterpolator(triangObject,density)
    pot_int = tcp1(PI,TI)
    den_int = tcp2(PI,TI)
    filtered_p = gaussian_filter(pot_int,1.8)
    filtered_d = gaussian_filter(den_int,1.8)
    #contour_phi_plot2(PI,TI,filtered_p,t,'phiovTe')
    #contour_phi_plot2(PI,TI,filtered_d,t,'lognovno)')
    contour_phi_plot3(TI,filtered_p,filtered_d,t)
    print("t is done!", t)
    
    
def keep_nodes(arr):
    '''Takes an array of nodes and keeps only those within a psi range.'''
    threshold_l = 0.9999#0.95
    threshold_u = 1.00001#1.05
    
    store = []
    for i in range(arr.shape[0]):
        if arr[i] > threshold_l:
            if arr[i] < threshold_u:
                store.append(i)
            else:
                pass
        else:
            pass
    return store

def keep_values(arr, keep_list):
    '''Takes a list of indices that need to be kept in another array and returns the second array with only those values.'''
    new_list = []        
    for i in range(len(keep_list)):
        new_list.append(arr[keep_list[i]])
        
    return new_list    
        
def smoothing_alg(L):
    LF = np.roll(L,1)
    LB = np.roll(L,-1)
    LFF = np.roll(L,2)
    LBB = np.roll(L,-2)
    newL = (2.0*(LF + LB)+LFF+LBB)/6.0
    return newL

def smooth(L,n):
    sL = L
    for i in range(0,int(n)):
        sL = smoothing_alg(sL)
    return sL


def write_mat():
    #B_R = np.array(getframe(bfield[:,0]))
    #np.save('/global/cscratch1/sd/giannos/ti344_arrays/B_R_mid', B_R)
    #print('BR done')
    #B_Z = np.array(getframe(bfield[:,1]))
    #np.save('/global/cscratch1/sd/giannos/ti344_arrays/B_Z_mid', B_Z)
    #print('BZ done')
    #B_zeta = np.array(getframe(bfield[:,2]))
    #np.save('/global/cscratch1/sd/giannos/ti344_arrays/B_zeta_mid', B_zeta)
    #print('Bzeta done')
    ne_all = np.array([[getframe(ne[:,iz,it]) for iz in range(0,Nplanes)] for it in range(0,200)])
    #np.save('/global/cscratch1/sd/giannos/ti344_arrays/ne_mid', ne_all)
    np.save('ti255_ne_short.npy',ne_all)
    print('ne done')
    