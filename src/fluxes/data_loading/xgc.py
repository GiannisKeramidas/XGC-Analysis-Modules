
"""Module of the XGC loader for general use, taken from Loic's load_XGC_local for BES.

It reads the data from the simulation and remove the points not wanted.
This file is based on an other code written by Lei Shi (:download:`code <../../../../FPSDP/Plasma/XGC_Profile/load_XGC_profile.py>`).

"""

import numpy as np
import os.path
import sys
import glob
from scipy.interpolate import splrep, splev
from scipy.interpolate import LinearNDInterpolator, CloughTocher2DInterpolator

#convenience gateway to load XGC1 or XGCa data
def load(*args,**kwargs):
    global file_path
    file_path = os.path.join(args[0],'')

    if len(glob.glob(file_path+'xgc.3d*')) > 0:
        return xgc1Load(*args,**kwargs)
    
    elif len(glob.glob(file_path+'xgc.2d*')) > 0:
        return xgcaLoad(*args,**kwargs)

    else:
        raise ValueError('XGC files not found in '+file_path)


class _load(object):
    """Loader class for general use.

    The idea of this loader is to load all data, limiting spatially and temporally as user-specified.

    INPUTS
    :param str xgc_path: Name of the directory containing the data

    OPTIONAL
    :param int t_start: Time step at which starting the diagnostics [1-based indexed (to match file)]
    :param int t_end: Time step at which stoping the diagnostics [1-based indexed (to match file)]
    :param int dt: Interval between two time step that the diagnostics should compute
    :param float Rmin: Limit minimum major radius of data
    :param float Rmax: Limit maximum major radius of data
    :param float Zmin: Limit minimum vertical coordinate of data
    :param float Zmax: Limit maximum vertical coordinate of data
    :param float psinMin: Limit minimum normalized poloidal flux ("psi-normal") of data
    :param float psinMax: Limit maximum normalized poloidal flux ("psi-normal") of data
    :param int phi_start: Toroidal plane index to start data [0-based indexed]
    :param int phi_end: Toroidal plane index to stop data [0-based indexed]


    :param str kind: Order of the interpolation method (linear or cubic))
    """

    def __init__(self,xgc_path,t_start=890,t_end=901,dt=1,
        Rmin=None,Rmax=None,Zmin=None,Zmax=None,
        psinMin=None,psinMax=None,phi_start=0, phi_end=None,
        kind='linear'):
        """Copy all the input values and call all the functions that compute the equilibrium and the first
        time step.
        """
        def readAdios(x,v):
            if '/' in v: v = '/'+v
            if type(x) is adios.file:
                return x[v][:]
            else:
                return adios.file(x+'.bp')[v][:]

        def readHDF5(x,v):
            return h5py.File(x+'.h5','r')[v][:]

        print('Loading XGC output data')
        
        self.xgc_path = os.path.join(xgc_path,'')  #get file_path, add path separator if not there
        self.mesh_file=self.xgc_path+'xgc.mesh'
        #check if files are in HDF5 or ADIOS format
        if os.path.exists(self.mesh_file+'.bp'):
            import adios
            self.readCmd=readAdios
        elif os.path.exists(self.mesh_file+'.h5'):
            import h5py
            self.readCmd=readHDF5
        else:
            raise ValueError('No xgc.mesh file found')

        print('from directory:'+ self.xgc_path)

        #read in units file
        self.unit_file = self.xgc_path+'units.m'
        self.unit_dic = self.load_m(self.unit_file)

        #read in time
        self.oneddiag_file=self.xgc_path+'xgc.oneddiag'
        self.mask1d = self.oned_mask()
        self.time = self.readCmd(self.oneddiag_file,'time')[self.mask1d]
        self.t_start=t_start
        self.t_end=t_end        
        if self.t_end is None: self.t_end=len(self.time)
        #self.time = self.time[(self.t_start-1):(self.t_end):dt]
        self.time = self.time[0:(self.t_end):dt]
        self.time_steps = np.arange(self.t_start,self.t_end+1,dt) #1-based for file names
        self.tstep = self.unit_dic['sml_dt']*self.unit_dic['diag_1d_period']
        self.dt = self.tstep * dt
        self.Ntimes = len(self.time)

        #magnetics file
        self.bfield_file=self.xgc_path+'xgc.bfield'

        # limits of the mesh in tokamak coordinates. Set to min,max of arrays in loadMesh()
        #if unspecified by user
        self.Rmin = self.unit_dic['eq_x_r'] if 'x' in str(Rmin).lower() else Rmin
        self.Rmax = self.unit_dic['eq_x_r'] if 'x' in str(Rmax).lower() else Rmax
        self.Zmin = self.unit_dic['eq_x_z'] if 'x' in str(Zmin).lower() else Zmin
        self.Zmax = self.unit_dic['eq_x_z'] if 'x' in str(Zmax).lower() else Zmax
        self.psinMin=psinMin
        self.psinMax=psinMax

        self.kind = kind
        
        
        #read in mesh, equilibrium data, and finally fluctuation data
        print('Loading mesh and psi...')
        self.loadMesh()
        print('\tmesh and psi loaded.')
        
        print('Loading magnetics...')
        self.loadBfield()
        print('\tmagnetics loaded.')
        
        print('Loading equilibrium...')
        self.loadEquil()
        print('\tequlibrium loaded.')

    
    def load_m(self,fname):
        """load the whole .m file and return a dictionary contains all the entries.
        """
        f = open(fname,'r')
        result = {}
        for line in f:
            words = line.split('=')
            key = words[0].strip()
            value = words[1].strip(' ;\n')
            result[key]= float(value)
        f.close()
        return result

    def loadMesh(self):
        """load R-Z mesh and psi values, then create map between each psi 
           value and the series of points on that surface.
        """
        # get mesh R,Z and psin
        RZ = self.readCmd(self.mesh_file,'coordinates/values')
        R=RZ[:,0]
        Z=RZ[:,1]
        psi = self.readCmd(self.mesh_file,'psi')
        psin = psi/self.unit_dic['psi_x']
        tri=self.readCmd(self.mesh_file,'cell_set[0]/node_connect_list') #already 0-based

        # set limits if not user specified
        if self.Rmin is None: self.Rmin=np.min(R)
        if self.Rmax is None: self.Rmax=np.max(R)
        if self.Zmin is None: self.Zmin=np.min(Z)
        if self.Zmax is None: self.Zmax=np.max(Z)
        if self.psinMin is None: self.psinMin=np.min(psin)
        if self.psinMax is None: self.psinMax=np.max(psin)

        #limit to the user-input ranges        
        self.rzInds = ( (R>=self.Rmin) & (R<=self.Rmax) & 
            (Z>=self.Zmin) & (Z<=self.Zmax) & 
            (psin>=self.psinMin) & (psin<=self.psinMax) )

        self.RZ = RZ[self.rzInds,:]
        self.psin = psin[self.rzInds]

        # psi interpolant
        fill_ = np.nan
        if self.kind == 'linear':
            self.psi_interp = LinearNDInterpolator(
                self.RZ, self.psin, fill_value=fill_)
        elif self.kind == 'cubic':
            self.psi_interp = CloughTocher2DInterpolator(
                self.RZ, self.psin, fill_value=fill_)
        else:
            raise NameError("The method '{}' is not defined".format(self.kind))        

        #get the triangles which are all contained within the vertices defined by
        #the indexes igrid
        #find which triangles are in the defined spatial region
        tmp=self.rzInds[tri] #rzInds T/F array, same size as R
        goodTri=np.all(tmp,axis=1) #only use triangles who have all vertices in rzInds
        self.tri=tri[goodTri,:]
        #remap indices in triangulation
        indices=np.where(self.rzInds)[0]
        for i in range(len(indices)):
            self.tri[self.tri==indices[i]]=i


    def loadEquil(self):
        """Load equilibrium profiles and compute the interpolant
        """
        #read in 1d psin data
        self.psin1d = self.readCmd(self.oneddiag_file,'psi')
        if self.psin1d.ndim > 1: self.psin1d = self.psin1d[0,:]

        #read electron temperature
        try:
          etemp_par=self.readCmd(self.oneddiag_file,'e_parallel_mean_en_avg')
          etemp_per=self.readCmd(self.oneddiag_file,'e_perp_temperature_avg')
        except:
          #etemp_par=self.readCmd(self.oneddiag_file,'e_parallel_mean_en_1d')
          #etemp_per=self.readCmd(self.oneddiag_file,'e_perp_temperature_1d')
          etemp_par=self.readCmd(self.oneddiag_file,'e_parallel_mean_en_df_1d')
          etemp_per=self.readCmd(self.oneddiag_file,'e_perp_temperature_df_1d')  
        self.Te1d=(etemp_par[self.mask1d,:]+etemp_per[self.mask1d,:])*2./3

        #read ion temperature
        try:
          itemp_par=self.readCmd(self.oneddiag_file,'i_parallel_mean_en_avg')
          itemp_per=self.readCmd(self.oneddiag_file,'i_perp_temperature_avg')
        except:
          #itemp_par=self.readCmd(self.oneddiag_file,'i_parallel_mean_en_1d')
          #itemp_per=self.readCmd(self.oneddiag_file,'i_perp_temperature_1d')
          itemp_par=self.readCmd(self.oneddiag_file,'i_parallel_mean_en_df_1d')
          itemp_per=self.readCmd(self.oneddiag_file,'i_perp_temperature_df_1d')  
        self.Ti1d=(itemp_par[self.mask1d,:]+itemp_per[self.mask1d,:])*2./3
        self.Tiperp=itemp_per[self.mask1d,:]
        self.Tipar=itemp_par[self.mask1d,:]

        #read flows for current
        '''try:
          i_pol_flow=self.readCmd(self.oneddiag_file,'i_gc_poloidal_flow_df_1d')
          e_pol_flow=self.readCmd(self.oneddiag_file,'e_gc_poloidal_flow_df_1d')
        except:
          i_pol_flow=self.readCmd(self.oneddiag_file,'i_gc_poloidal_flow_1d')
          e_pol_flow=self.readCmd(self.oneddiag_file,'e_gc_poloidal_flow_1d')
        self.i_pol_f=i_pol_flow[self.mask1d,:]
        self.e_pol_f=e_pol_flow[self.mask1d,:]'''
        
        #read electron particle fluxes
        '''try:
            e_radial_flux = self.readCmd(self.oneddiag_file, 'e_radial_flux_df_1d')
            e_radial_flux_EXB = self.readCmd(self.oneddiag_file, 'e_radial_flux_ExB_df_1d')
        except:
            e_radial_flux = self.readCmd(self.oneddiag_file, 'e_radial_flux_1d')
            e_radial_flux_EXB = self.readCmd(self.oneddiag_file, 'e_radial_flux_ExB_1d')
        self.e_rad_fl = e_radial_flux[self.mask1d,:]
        self.e_rad_exb_fl = e_radial_flux_EXB[self.mask1d,:]'''
        
        #read ion particle fluxes
        try:
            i_radial_flux = self.readCmd(self.oneddiag_file, 'i_radial_flux_df_1d')
            i_radial_flux_EXB = self.readCmd(self.oneddiag_file, 'i_radial_flux_ExB_df_1d')
        except:
            i_radial_flux = self.readCmd(self.oneddiag_file, 'i_radial_flux_1d')
            i_radial_flux_EXB = self.readCmd(self.oneddiag_file, 'i_radial_flux_ExB_1d')
        self.i_rad_fl = i_radial_flux[self.mask1d,:]
        self.i_rad_exb_fl = i_radial_flux_EXB[self.mask1d,:]
        
        
        
        #read electron heat fluxes
        '''try:
            e_radial_en_flux = self.readCmd(self.oneddiag_file, 'e_radial_en_flux_df_1d')
            e_radial_en_flux_EXB = self.readCmd(self.oneddiag_file, 'e_radial_en_flux_ExB_df_1d')
        except:
            e_radial_en_flux = self.readCmd(self.oneddiag_file, 'e_radial_en_flux_1d')
            e_radial_en_flux_EXB = self.readCmd(self.oneddiag_file, 'e_radial_en_flux_ExB_1d')
        self.e_rad_en_fl = e_radial_en_flux[self.mask1d,:]
        self.e_rad_en_exb_fl = e_radial_en_flux_EXB[self.mask1d,:]'''
        
        #read ion heat fluxes
        '''try:
            i_radial_en_flux = self.readCmd(self.oneddiag_file, 'i_radial_en_flux_df_1d')
            i_radial_en_flux_EXB = self.readCmd(self.oneddiag_file, 'i_radial_en_flux_ExB_df_1d')
        except:
            i_radial_en_flux = self.readCmd(self.oneddiag_file, 'i_radial_en_flux_1d')
            i_radial_en_flux_EXB = self.readCmd(self.oneddiag_file, 'i_radial_en_flux_ExB_1d')
        self.i_rad_en_fl = i_radial_en_flux[self.mask1d,:]
        self.i_rad_en_exb_fl = i_radial_en_flux_EXB[self.mask1d,:]'''
        
        

        #read electron density
        #self.ne1d = self.readCmd(self.oneddiag_file,'e_gc_density_1d')[self.mask1d,:]
        self.ne1d = self.readCmd(self.oneddiag_file,'e_gc_density_df_1d')[self.mask1d,:]

        #read n=0,m=0 potential
        try:
            self.psin001d = self.readCmd(self.oneddiag_file,'psi00_1d')/self.unit_dic['psi_x']
        except:
            self.psin001d = self.readCmd(self.oneddiag_file,'psi00')/self.unit_dic['psi_x']
        if self.psin001d.ndim > 1: self.psin001d = self.psin001d[0,:]
        self.pot001d = self.readCmd(self.oneddiag_file,'pot00_1d')[self.mask1d,:]
        
        #create splines for t=0 data
        self.te0_sp = splrep(self.psin1d,self.Te1d[0,:],k=1)
        self.ne0_sp = splrep(self.psin1d,self.ne1d[0,:],k=1)

    def loadBfield(self):
        """Load magnetic field
        """

        self.bfield = self.readCmd(self.bfield_file,'node_data[0]/values')[self.rzInds,:]

    def oned_mask(self):
        """Match oned data to 3d files, in cases of restart.
           Use this on oneddiag variables, e.g. n_e1d = ad.file('xgc.oneddiag.bp','e_gc_density_1d')[mask1d,:]
        """
        try:
            step = self.readCmd(self.oneddiag_file,'step')
            dstep = step[1] - step[0]
        
            idx = np.arange(step[0]/dstep,step[-1]/dstep+1)
        
            mask1d = np.zeros(idx.shape,dtype=np.int32)
            for i in idx:
                mask1d[i-1] = np.where(step == i*dstep)[0][-1] #get last occurence
        except:
            mask1d = Ellipsis #pass variables unaffected
        
        return mask1d


class xgc1Load(_load):
    def __init__(self,xgc_path,phi_start=0,phi_end=None,**kwargs):
        #call parent loading init, including mesh and equilibrium
        #super().__init__(*args,**kwargs)
        super(xgc1Load,self).__init__(xgc_path,**kwargs)

        #read in number of planes
        fluc_file0 = self.xgc_path + 'xgc.3d.' + str(self.time_steps[0]).zfill(5)
        self.Nplanes=self.readCmd(fluc_file0,'dpot').shape[1]
        self.phi_start=phi_start
        self.phi_end = phi_end
        if phi_end is None: self.phi_end=self.Nplanes-1
        self.Nplanes=self.phi_end-self.phi_start+1

        
        '''if len(glob.glob(file_path+'xgc.f3d*')) > 0:
            print('Loading f3d files...')
            self.loadf3d()
            print('f3d files loaded')
        else:
            pass'''
        
        print('Loading fluctuations...')
        self.loadFluc()
        print('fluctuations loaded')

    def loadFluc(self):
        """Load non-adiabatic electron density, electrical static 
        potential fluctuations, and n=0 potential for 3D mesh.
        The required planes are calculated and stored in sorted array.
        fluctuation data on each plane is stored in the same order.
        Note that for full-F runs, the perturbed electron density 
        includes both turbulent fluctuations and equilibrium relaxation,
        this loading method doesn't differentiate them and will read all of them.
        
        """
        self.eden = np.zeros( (len(self.RZ[:,0]), self.Nplanes, self.Ntimes-self.t_start) )
        
        self.dpot = np.zeros( (len(self.RZ[:,0]), self.Nplanes, self.Ntimes-self.t_start) )
        self.pot0 = np.zeros( (len(self.RZ[:,0]), self.Ntimes-self.t_start) )
        
        for i in range(self.t_start,self.t_end):
            flucFile = self.xgc_path + 'xgc.3d.'+str(i).zfill(5)
            sys.stdout.write('\r\tLoading file ['+str(i)+'/'+str(self.Ntimes)+']')
            self.dpot[:,:,i-self.t_start] = self.readCmd(flucFile,'dpot')[self.rzInds,self.phi_start:(self.phi_end+1)]
            self.pot0[:,i-self.t_start] = self.readCmd(flucFile,'pot0')[self.rzInds]
            self.eden[:,:,i-self.t_start] = self.readCmd(flucFile,'eden')[self.rzInds,self.phi_start:(self.phi_end+1)]
            
        
        if self.Nplanes == 1:
            self.dpot = self.dpot.squeeze()
            self.eden = self.eden.squeeze()

    def loadf3d(self):
        
        self.i_T_perp = np.zeros((len(self.RZ[:,0]),self.Nplanes, self.Ntimes-self.t_start))
        self.i_E_para = np.zeros((len(self.RZ[:,0]),self.Nplanes, self.Ntimes-self.t_start))
        self.e_T_perp = np.zeros((len(self.RZ[:,0]),self.Nplanes, self.Ntimes-self.t_start))
        self.e_E_para = np.zeros((len(self.RZ[:,0]),self.Nplanes, self.Ntimes-self.t_start))
        self.i_u_para = np.zeros((len(self.RZ[:,0]),self.Nplanes, self.Ntimes-self.t_start))
        self.e_u_para = np.zeros((len(self.RZ[:,0]),self.Nplanes, self.Ntimes-self.t_start))
        
        for i in range(self.t_start,self.t_end):
            flucFilef = self.xgc_path + 'xgc.f3d.'+str(i).zfill(5)
            sys.stdout.write('\r\tLoading file ['+str(i)+'/'+str(self.Ntimes)+']')
            self.i_T_perp[:,:,i-self.t_start] = self.readCmd(flucFilef,'i_T_perp')[self.rzInds,self.phi_start:(self.phi_end+1)]
            self.i_E_para[:,:,i-self.t_start] = self.readCmd(flucFilef,'i_E_para')[self.rzInds,self.phi_start:(self.phi_end+1)]
            self.e_T_perp[:,:,i-self.t_start] = self.readCmd(flucFilef,'e_T_perp')[self.rzInds,self.phi_start:(self.phi_end+1)]
            self.e_E_para[:,:,i-self.t_start] = self.readCmd(flucFilef,'e_E_para')[self.rzInds,self.phi_start:(self.phi_end+1)]
            self.i_u_para[:,:,i-self.t_start] = self.readCmd(flucFilef,'i_u_para')[self.rzInds,self.phi_start:(self.phi_end+1)]
            self.e_u_para[:,:,i-self.t_start] = self.readCmd(flucFilef,'e_u_para')[self.rzInds,self.phi_start:(self.phi_end+1)]
        
        if self.Nplanes == 1:
            self.i_T_perp = self.i_T_perp.squeeze()
            self.i_E_para = self.i_E_para.squeeze()
            self.e_T_perp = self.e_T_perp.squeeze()
            self.e_E_para = self.e_E_para.squeeze()
            self.i_u_para = self.i_u_para.squeeze()
            self.e_u_para = self.e_u_para.squeeze()



    def calcNeTotal(self,psin=None):
        """Calculate the total electron at the wanted points.

        :param np.array[N] psin

        :returns: Total density
        :rtype: np.array[N]
        
        """
        
        if psin is None: psin=self.psin

        # temperature and density (equilibrium) on the psi mesh
        te0 = splev(psin,self.te0_sp)
        # avoid temperature <= 0
        te0[te0<np.min(self.Te1d)/10] = np.min(self.Te1d)/10
        ne0 = splev(psin,self.ne0_sp)
        ne0[ne0<np.min(self.ne1d)/10] = np.min(self.ne1d)/10
        

        #neAdiabatic = ne0*exp(dpot/te0)
        factAdiabatic = np.exp(np.einsum('i...,i...->i...',self.dpot,1./te0))
        #factAdiabatic = 1+np.einsum('i...,i...->i...',self.dpot,1./te0)
        self.neAdiabatic = np.einsum('i...,i...->i...',ne0,factAdiabatic)

        #ne = neAdiatbatic + dneKinetic
        self.n_e = self.neAdiabatic + self.eden

        #TODO I've ignored checking whether dne<<ne0, etc. may want to add
        return self.n_e
    
    def expo_arg(self,psin=None):
        if psin is None: psin=self.psin

        # temperature and density (equilibrium) on the psi mesh
        te0 = splev(psin,self.te0_sp)
        # avoid temperature <= 0
        te0[te0<np.min(self.Te1d)/10] = np.min(self.Te1d)/10
        arg = np.einsum('i...,i...->i...',self.dpot,1./te0)
        return arg

    def calcNeAdiabatic(self,psin=None):
        """Calculate the adiabatic electron density at the wanted points.

        :param np.array[N] psin

        :returns: Adiabatic density
        :rtype: np.array[N]
        
        """
        
        if psin is None: psin=self.psin

        # temperature and density (equilibrium) on the psi mesh
        te0 = splev(psin,self.te0_sp)
        # avoid temperature <= 0
        te0[te0<np.min(self.Te1d)/10] = np.min(self.Te1d)/10
        ne0 = splev(psin,self.ne0_sp)
        ne0[ne0<np.min(self.ne1d)/10] = np.min(self.ne1d)/10
        

        #neAdiabatic = ne0*exp(dpot/te0)
        factAdiabatic = np.exp(np.einsum('i...,i...->i...',self.dpot,1./te0))
        #factAdiabatic = 1+np.einsum('i...,i...->i...',self.dpot,1./te0)
        self.neAdiabatic = np.einsum('i...,i...->i...',ne0,factAdiabatic)

        
        #TODO I've ignored checking whether dne<<ne0, etc. may want to add
        return self.neAdiabatic

    
    
    def calcPotential(self):
        self.pot = self.pot0[:,np.newaxis,:] + self.dpot
        return self.pot



class xgcaLoad(_load):
    def __init__(self,xgc_path,**kwargs):
        #call parent loading init, including mesh and equilibrium
        #super().__init__(*args,**kwargs)
        super(xgcaLoad,self).__init__(xgc_path,**kwargs)

        print('Loading f0 data...')
        self.loadf0()
        print('f0 data loaded')

    def loadf0(self):
        ##f0 mesh data
        self.f0_file = self.xgc_path+'xgc.f0.mesh'
        #load velocity grid parallel velocity
        f0_nvp = self.readCmd(self.f0_file,'f0_nvp')
        self.nvpa = 2*f0_nvp+1 #actual # of Vparallel velocity pts (-vpamax,0,vpamax)
        self.vpamax = self.readCmd(self.f0_file,'f0_vp_max')
        self.vpar = np.linspace(-self.f0_vp_max,self.vpamax,self.nvpa)
        #load velocity grid perpendicular velocity
        f0_nmu = self.readCmd(self.f0_file,'f0_nmu')
        self.nvpe = f0_nmu + 1 #actual # of Vperp velocity pts (0,vpemax)
        self.vpemax = self.readCmd(self.f0_file,'f0_smu_max')
        self.vperp = np.linspace(0,self.f0_smu_max,self.nvpe)
        #load velocity grid density
        self.f0_ne = self.readCmd(self.f0_file,'f0_den')
        #load velocity grid electron and ion temperature
        f0_t_ev = self.readCmd(self.f0_file,'f0_t_ev')
        self.f0_Te = f0_t_ev[0,:]
        self.f0_Ti = f0_t_ev[1,:]

        self.f0_grid_vol_vonly = self.readCmd(self.f0_file,'f0_grid_vol_vonly')

    def load2D():
        self.iden = np.zeros( (len(self.RZ[:,0]), self.Ntimes) )
        
        self.dpot = np.zeros( (len(self.RZ[:,0]), self.Ntimes) )
        self.pot0 = np.zeros( (len(self.RZ[:,0]), self.Ntimes) )
        self.epsi = np.zeros( (len(self.RZ[:,0]), self.Ntimes) )
        self.etheta = np.zeros( (len(self.RZ[:,0]), self.Ntimes) )
        
        
        for i in range(self.t_start,self.t_end+1):
            twodFile = self.xgc_path + 'xgc.2d.'+str(i).zfill(5)

            self.iden[:,i-1] = self.readCmd(flucFile,'iden')[self.rzInds]

            self.dpot[:,i-1] = self.readCmd(flucFile,'dpot')[self.rzInds]
            self.pot0[:,i-1] = self.readCmd(flucFile,'pot0')[self.rzInds]
            self.epsi[:,i-1] = self.readCmd(flucFile,'epsi')[self.rzInds]
            self.etheta[:,i-1] = self.readCmd(flucFile,'etheta')[self.rzInds]
            

        


class gengridLoad():
    def __init__(self,file_path):
        self.file_path = file_path

        #read in the grid data from node file
        #TODO Read in ele and poly components also
        f = open(self.file_path,'r')
        Nlines = int(f.readline().split()[0])
        self.RZ = np.empty([Nlines,2])
        for i,line in enumerate(f):
            if i >= Nlines: break
            self.RZ[i,0:2] = np.array(line.split()[1:3],dtype='float')