import numpy as np
import xgc
import math
import pylab as py
from matplotlib.tri import Triangulation, LinearTriInterpolator
import matplotlib.pyplot as plt
import glob
import angle
import numpy.ma as ma
from decimal import Decimal
import core
import Reyn

#set up parameters for plotting:
plt.style.use('default')
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.labelweight'] = 'heavy'
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12
plt.rcParams['text.usetex']=True


#core.getMeshAndCuts(core.fileDir,None,None,None,None)


def loading_fun():
    '''Function that reads the 2D values of psin,R,B,Bp,b_zeta and phi.'''
    global psin,unitR, unitZ, phi, B, Bp, Bzetai, R, Z, Rmin,Rmax,Zmin,Zmax
    option = 1
    Rmax = np.amax(core.RZ[:,0])
    Rmin = np.amin(core.RZ[:,0])
    Zmax = np.amax(core.RZ[:,1])
    Zmin = np.amin(core.RZ[:,1])
    unitR = (Rmax-Rmin)/len(core.Ri)
    unitZ = (Zmax-Zmin)/len(core.Zi)
    
    phi = np.asarray([[core.getcutvalue_hor(core.pot[:,iz,0],cut,option) for iz in range(0,core.Nplanes)]for cut in range(0,len(core.Zi))]).mean(axis=1)
    psin = np.asarray([core.getcutvalue_hor(core.psin[:],cut,option) for cut in range(0,len(core.Zi))])
    BRi = np.asarray([core.getcutvalue_hor(core.bfield[:,0],cut,option) for cut in range(0,len(core.Zi))])
    BZi = np.asarray([core.getcutvalue_hor(core.bfield[:,1],cut,option)for cut in range(0,len(core.Zi))])
    Bzetai = np.asarray([core.getcutvalue_hor(core.bfield[:,2],cut,option)for cut in range(0,len(core.Zi))])
    Bp = np.sqrt(np.square(BRi)+np.square(BZi))
    B = np.sqrt(np.square(BRi)+np.square(BZi)+np.square(Bzetai))
    R = np.array([ip*unitR+Rmin for ip in range(0,len(core.Zi))])
    Z = np.array([ip*unitZ+Zmin for ip in range(0,len(core.Zi))])
    
    return psin, phi, B, Bp, Bzetai, R, Z


def magn_axis():
    '''Function that locates coordinates of magnetic axis'''
    temp = map(np.nanmin,psin)#creates a 1D array with minimums of each row. The row number is the Z-number.
    Z_loc = temp.index(min(temp))#finds the cut where the zero of psin is located
    temp2 = psin[Z_loc,:].tolist()
    R_loc = temp2.index(np.nanmin(temp2))#finds the location of the R-minimum
    Z_pos = Z_loc*unitZ+Zmin
    R_pos = R_loc*unitR+Rmin
    ntemp2 = np.asarray(temp2)
    ntemp2 = ntemp2[~np.isnan(ntemp2)]
    #print(ntemp2)
    R_sep_loc = np.where(np.abs(ntemp2-1.0)<1e-1)[0]
    #diff_arr = np.abs(ntemp2-1.0).tolist() #for when we run with more accuracy
    #R_sep_loc = diff_arr.index(min(diff_arr))
    #print('Rsep=',R_sep_loc[-1])
    return R_pos,Z_pos, R_loc, Z_loc, R_sep_loc[-1]


def evaluation_at_tip(R_star,Z_star):
    '''Function that evaluates quantities at the banana tip.'''
    #_, _, _ , Z_star,_ = magn_axis()
    diff_mat_R = np.abs(R-R_star)
    diff_mat_Z = np.abs(Z-Z_star)
    diff_list_R = diff_mat_R.tolist()
    diff_list_Z = diff_mat_Z.tolist()
    R_star_loc = diff_list_R.index(min(diff_list_R))
    Z_star_loc = diff_list_Z.index(min(diff_list_Z))
    B_star = B[Z_star_loc,R_star_loc]
    phi_star = phi[Z_star_loc,R_star_loc]
    
    return B_star,phi_star


def psi_der(cut):
    
    dZ = core.Zi[1] - core.Zi[0]
    dR = core.Ri[1] - core.Ri[0]
    
    #d/dZ of psi
    arrim1 = psin[(cut-1)%len(core.Zi),:]
    arri = psin[cut,:]
    arrip1 = psin[(cut+1)%len(core.Zi),:]   
    darrdZ = (arrip1-arrim1)/(2.*dZ)
    
     #d/dR of psi
    arriplus = np.roll(arri,-1,axis=0)
    arriminus = np.roll(arri,1,axis=0)
    darrdR = (arriplus-arriminus)/(2.*dR)
    darrdR[0] = (arri[1]-arri[0])/dR
    darrdR[-1] = (arri[-1]-arri[-2])/dR

    return darrdZ,darrdR


def two_d_psi_der():
    '''Function that calculates dpsin/dR, dpsin/dZ, abs(grad psin).'''
    dpsidZ,dpsidR = zip(*[psi_der(cut) for cut in range(0,len(core.Zi))])
    dPsdZ = np.asarray(dpsidZ)
    dPsdR = np.asarray(dpsidR)
    
    #abs(grad psin)
    abs_grad_psin = abs(np.sqrt(np.square(dPsdZ) + np.square(dPsdR)))
    
    return dPsdZ,dPsdR,abs_grad_psin


def Dpsin_Int():
    '''Function that performs the integral for Dpsin. The integral is a dR one and starts where psin=0 and ends at the separatrix where
    psin=1. It is performed following the line on the cut that corresponds to psin=0 (Z_st_loc)'''
    _,_,R_st_loc,Z_st_loc,R_sep_loc = magn_axis()
    R_end_loc = R_sep_loc
    prefactor = (R[R_end_loc]-R[R_st_loc])/(len(R[R_st_loc:R_end_loc])-1)
    _,dPsdR,abs_grad_psin = two_d_psi_der()
    weight = ((R[R_st_loc:R_end_loc]*Bp[Z_st_loc,R_st_loc:R_end_loc])/abs_grad_psin[Z_st_loc,R_st_loc:R_end_loc])*dPsdR[Z_st_loc,R_st_loc:R_end_loc]
  
    Summation = (1/2)*weight[0]+(1/2)*weight[-1]+np.sum(weight[1:-2])
    Integral = prefactor*Summation
    return Integral


def banana_width(R_star,W_star):
    '''Function that evaluates the babana width. Takes banana tip, R_star in (m) and ion energy, W_star in (eV).'''
    global A_i, Z_i
    A_i = 2.0 #ion mass number
    Z_i = 1.0 #ion charge 
    
    ban_width = (1.44*1.0e-4)*(R_star*np.sqrt(W_star)*np.sqrt(A_i))/(Dpsin_Int()*Z_i)
    
    return ban_width


def psi_star(R_star,Z_star,W_star):
    '''Function that evaluates normalized toroidal momentum and returns the positive and negative branch.'''
    e= 0.0
    B_star,phi_star = evaluation_at_tip(R_star,Z_star)
    psin_star_p = psin + banana_width(R_star,W_star)*((R*(Bzetai/B))/R_star)*np.sqrt(1-(B/B_star)+(Z_i*e*(phi_star-phi))/W_star)
    psin_star_m = psin - banana_width(R_star,W_star)*((R*(Bzetai/B))/R_star)*np.sqrt(1-(B/B_star)+(Z_i*e*(phi_star-phi))/W_star)
    
    return psin_star_p,psin_star_m
    
def cont_test(R_star,Z_star,W_star):
    e=0.2
    B_star,phi_star = evaluation_at_tip(R_star,Z_star)
    p = 1-(B/B_star)+(Z_i*e*(phi_star-phi))/W_star
    fig, ax =plt.subplots(figsize=(10,10))
    CS = plt.contour(R, Z, p,colors='b',linewidths=2)
    plt.clabel(CS)
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.show()

    
def contour_banana_plot(R_star,Z_star,W_star):
    '''Contour plot of the normalized toroidal momentum branches, i.e., the banana orbits.'''
    psin_star_p,psin_star_m = psi_star(R_star,Z_star,W_star)
    psi_list = []
    for cut in range(len(core.Zi)):
        psi_list.append(core.getcutvalue_hor(core.psin[:],cut,1))
    psi_mat = np.reshape(psi_list,(len(core.Ri),len(core.Zi)))
    fig, ax =plt.subplots(figsize=(10,10))
    contour_levels = np.arange(0.9,1.2,0.05)
    CS = plt.contour(R, Z, psi_mat,contour_levels,colors='y')
    #CB = plt.colorbar(CS, extend='both')
    p_levels = np.arange(-0.06468,1.6884,0.5)
    m_levels = np.arange(-0.06468,1.9245,0.5)
    CP = plt.contour(R, Z, psin_star_p,p_levels,colors='b',linewidths=2)
    CM = plt.contour(R, Z, psin_star_m,m_levels,colors='r',linewidths=2)
    #print("pmax",np.nanmax(psin_star_p))
    #print("pmin",np.nanmin(psin_star_p))
    #print("mmax",np.nanmax(psin_star_m))
    #print("mmin",np.nanmin(psin_star_p))
    #CB = plt.colorbar(CP, extend='both')
    #CD = plt.colorbar(CM, extend='both')
    #plt.clabel(CP, inline=1,fmt='%.2E', fontsize=10)
    plt.xlabel("R(m)")
    plt.ylabel("Z(m)")
    ax.plot(core.loader.RZ[core.sepInds,0],core.loader.RZ[core.sepInds,1],'k--')
    plt.title('Banana Orbits, $R_\star = %s,\,\, E_i = %s\, eV$'%(R_star,W_star))
    plt.show()




'''Functions that evaluate the banana plots from the nodes on an psi-theta plane.'''    
def prepare_data():
    '''Loads the node data of the relevant quantities and keeps only the nodes within a psi range.'''
    global psi, thetas, potential, B_mag, Bzeta, Radius
    potential = core.pot[:,:,0]
    potential = potential.mean(axis=1)
    BR = core.bfield[:,0]
    BZ = core.bfield[:,1]
    Bzeta = core.bfield[:,2]
    Radius = core.RZ[:,0]
    psi = core.psin[:]
    thetas = np.arctan2(core.RZ[:,1],core.RZ[:,0]-core.Rmaj)
    keep = keep_nodes(psi)
    potential = keep_values(potential,keep)
    psi = keep_values(psi,keep)
    thetas = keep_values(thetas,keep)
    BR = keep_values(BR,keep)
    BZ = keep_values(BZ,keep)
    Bzeta = keep_values(Bzeta,keep)
    Radius = keep_values(Radius,keep)
    B_mag = np.sqrt(np.square(BR)+np.square(BZ)+np.square(Bzeta))
    return psi, thetas, potential, B_mag, Bzeta, Radius

def psi_star_nodes(star_node,W_star):#give the core.sepInds node that the tip of the banana lies.
    '''Function that evaluates normalized toroidal momentum and returns the positive and negative branch.'''
    #evaluating values at the banana tip
    R_star = core.RZ[core.sepInds[star_node],0]
    B_star_R = core.bfield[core.sepInds[star_node],0]
    B_star_Z = core.bfield[core.sepInds[star_node],1]
    B_star_zeta = core.bfield[core.sepInds[star_node],2]
    B_star = np.sqrt(np.square(B_star_R)+np.square(B_star_Z)+np.square(B_star_zeta))
    phi_star = core.pot[core.sepInds[star_node],:,0]
    phi_star = phi_star.mean()
   
    #loading data on the nodes
    psi, thetas, potential, B_mag, Bzeta, Radius = prepare_data()
    
    #Include the electrostatic potential
    e= 1.0
    
    psin_star_p = psi + banana_width(R_star,W_star)*((Radius*(Bzeta/B_mag))/R_star)*np.sqrt(1-(B_mag/B_star)+(Z_i*e*(phi_star-potential))/W_star)
    psin_star_m = psi - banana_width(R_star,W_star)*((Radius*(Bzeta/B_mag))/R_star)*np.sqrt(1-(B_mag/B_star)+(Z_i*e*(phi_star-potential))/W_star)
    
    return psin_star_p,psin_star_m


def triang_and_interp(star_node, W_star):
    '''Does a triangulation on the psi-theta space, performs an interpolation there and then outputs the plot.'''
    from scipy.ndimage import gaussian_filter
    psi_star_p, psi_star_m = psi_star_nodes(star_node, W_star)
    Pi = np.linspace(min(psi),max(psi),400)
    Ti = np.linspace(min(thetas),max(thetas),400)
    (PI,TI) = np.meshgrid(Pi,Ti)
    triangObject = Triangulation(psi,thetas)
    tcp = LinearTriInterpolator(triangObject,psi_star_p)
    tcm = LinearTriInterpolator(triangObject,psi_star_m)
    tcp_int = tcp(PI,TI)
    tcm_int = tcm(PI,TI)
    p_filt = gaussian_filter(tcp_int,1.5)
    m_filt = gaussian_filter(tcm_int,1.5)
    psi_theta_plot(PI,TI,p_filt, m_filt, W_star)
        
def psi_theta_plot(psi, theta, P, M, W_star):
    plt.figure()
    #p_levels = [0.992]#np.arange(0.999,1.99,0.01)
    #m_levels = [0.992]#np.arange(0.999,1.99,0.01)
    CP = plt.contour(psi, theta, P, colors='b',linewidths=2)
    #plt.clabel(CP)
    CM = plt.contour(psi, theta, M, colors='r',linewidths=2)
    #plt.clabel(CM)
    plt.xlabel(r"$\Psi_N$")
    plt.ylabel(r"$\theta\, (rad)$")
    plt.title('C-Mod Banana Orbits, $E_i = %s\, eV$'%(W_star))
    plt.show()    
    
def keep_nodes(arr):
    '''Takes an array of nodes and keeps only those within a psi range.'''
    threshold_l = 0.95
    threshold_u = 1.05
    
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
