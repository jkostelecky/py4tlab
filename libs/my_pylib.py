#!/usr/bin/env python3
import os
import sys 
import glob
import my_pylib_tlab as mpt
import numpy as np
from   numpy import linalg as LA
import math as math
import subprocess
from scipy import integrate
from scipy.special import erf
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import netCDF4  as nc 
import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 

"""
#############################################################################
# Version 0.01 - Written by Jonathan Kostelecky
#---------------------------------------------------------------------------#
Python library for tlab.

Created  on & by: 2022/03/29 - J. Kostelecky (j.kostelecky@posteo.de)

Modified on & by: 2023/10/26 - J. Kostelecky
Modified on & by: ...        - ...
#---------------------------------------------------------------------------#
Contains:
            - classes & funtions for tlab standard data 
              and pre-/postprocessing
#############################################################################
"""

#############################################################################
#---------------------------------------------------------------------------#
# read 3d binary fields (flow, scalar, pressure, geometry)
#
# Variables:
#   - path : path to file(s)
#   - var  : u, v, w, phi1, phi2, scal, (phi1, phi2), flow (all velocity fields)
#            eps (IBM; format: bin, int, real), pre (pressure), etc. 
#            (cf. in field_name)
#   - index: index of the fields
#
#---------------------------------------------------------------------------#
#############################################################################
class Field:
    """
        Read and process binary 3D-field(s) from tlab simulation.
        Header: 5-numbers with dtype '<i4' (headerlength(52bits),imax,jmax,kmax,iterationstep )
                4-numbers with dtype '<f8' (viscosity, rtime, visc(1/reynolds), froude, rossby)
    """
    #-----------------------------------------------------------------------#
    def __init__(self, path='', var='flow', index=1, forcing='none'):
        # params
        self.path    = path
        self.var     = var
        self.index   = index
        self.forcing = forcing
        self.data    = dict()
        # data types
        self.type_i1           = np.dtype('<i1') 
        self.type_i4           = np.dtype('<i4') 
        self.type_f8           = np.dtype('<f8')
        self.sizeofdata_int1   = 1 
        self.sizeofdata_int4   = 4
        self.sizeofdata_float  = 8
        # header composition
        self.head_params_int   = 5 
        self.head_params_float = 4
        # file names
        self.fname = self.field_name(self.path, self.var, self.index)
        # get header informations
        self.head_int = self.read_binary(self.fname[0],self.type_i4,self.head_params_int,0)
        if not (var=='eps_int' or var=='eps_bit'): 
            self.head_float = self.read_binary(self.fname[0],self.type_f8,self.head_params_float,self.head_params_int*self.sizeofdata_int4)
            if var=='scal':
                self.head_float2 = self.read_binary(self.fname[1],self.type_f8,self.head_params_float,self.head_params_int*self.sizeofdata_int4)
        # print header
        print("--------------------------------------------------")  
        print("Reading 3d - fields   :  " , str(var))
        self.headersize  = self.head_int[0];    print('Headerlength          :   ' + str(self.headersize) + ' - Bytes')
        self.nx          = self.head_int[1]; 
        if var=='eps_bit': self.nx = self.nx * 8   
        self.ny          = self.head_int[2];
        self.nz          = self.head_int[3];    print('Gridsize (nx*ny*nz)   :   ' + str(self.nx)+'*'+str(self.ny)+'*'+str(self.nz))
        self.points      = np.int64(self.nx) * np.int64(self.ny)  * np.int64(self.nz) 
        self.step        = self.head_int[4];    print('Iteration step nr.    :   ' + str(self.step))
        self.rho         = 1.;                  print('Density               :   ' + str(self.rho)) 
        if not var[:3]=='eps':
            self.rtime   = self.head_float[0];  print('Physical sim. time    :   ' + str(self.rtime))
            self.visc    = self.head_float[1];  print('Viscositiy (1/Re)     :   ' + str(self.visc))
            self.re      = 1/self.visc;         print('Reynolds    number    :   ' + str(self.re))
            if var=='phi1' or var=='phi2':
                self.sc  = self.head_float[2];  print('Schmidt     number    :   ' + str(self.sc))
                self.da  = self.head_float[3];  print('Damkoehler  number    :   ' + str(self.da))
            elif var=='scal':
                self.sc2 = self.head_float2[2]; print('Schmidt2    number    :   ' + str(self.sc2))
                self.da2 = self.head_float2[3]; print('Damkoehler2 number    :   ' + str(self.da2))
            else:
                self.fr  = self.head_float[2];  print('Froude      number    :   ' + str(self.fr))
                self.ro  = self.head_float[3];  print('Rossby      number    :   ' + str(self.ro))
        if self.forcing == 'cpg':
            self.re_cl  = self.re;                    print('Centerl. Reynolds nr. :   ' + str(self.re_cl ))
            if self.forcing == 'cpg_open': self.re_tau = ((4/3) / 0.116)**-1 * self.re_cl ** 0.88
            else: self.re_tau = 0.116 * self.re_cl ** 0.88
            print('Friction Reynolds nr. :   ' + str(self.re_tau))
        if var=='eps_bit' :  byte_size = self.sizeofdata_int1/8
        elif var=='eps_int': byte_size = self.sizeofdata_int1
        else:                byte_size = self.sizeofdata_float
        print('Binary 3D-field size  :   ' + str(round((self.points*byte_size + self.headersize)/10**6,1)) + ' - Megabytes')
        return
    #-----------------------------------------------------------------------#
    def read_binary(self, fname, t, count, seek):
        f = open(fname,'rb')
        f.seek(seek,0)
        rec = np.fromfile(f, t, count)
        f.close()
        return rec
    #-----------------------------------------------------------------------#   
    def remove_mean(self, mean):  
        # input and output name
        i_fname = self.fname[0]
        o_fname = self.field_name(self.path, self.var, self.index+1)[0]
        # open input-/output-files
        i_f = open(i_fname,'rb') # input 
        o_f = open(o_fname,'wb') # output 
        # write header of ofile
        print('--------------------------------------------------')
        print('Remove mean of file   : ', i_fname)
        print('Write  to file        : ', o_fname)
        data = np.array([self.head_int],   self.type_i4); data.tofile(o_f)
        data = np.array([self.head_float], self.type_f8); data.tofile(o_f)
        # replace mean
        mean = np.array(mean)        # remove possible mask
        i_f.seek(self.headersize, 0) # jump over header
        nxy = self.nx*self.ny
        # buffered io
        for k in range(self.nz):
            # read input xy - plane
            i_dat = np.fromfile(i_f, self.type_f8, nxy)
            i_dat = i_dat.reshape((self.nx,self.ny), order='F')
            # replace mean
            i_dat = (i_dat - i_dat.mean(axis=(0))[np.newaxis, :]) + mean[np.newaxis, :]
            # write new file
            i_dat = i_dat.reshape((nxy), order='F')
            i_dat[:].astype(self.type_f8).tofile(o_f)
        # close files
        i_f.close()
        o_f.close()
        print('Write  to file        : DONE')
        return
    #-----------------------------------------------------------------------#   
    def smooth_upper(self, y, mean):  
        # input and output name
        i_fname = self.fname[0]
        o_fname = self.field_name(self.path, self.var, self.index+1)[0]
        # open input-/output-files
        i_f = open(i_fname,'rb') # input 
        o_f = open(o_fname,'wb') # output 
        # write header of ofile
        print('--------------------------------------------------')
        print('Remove mean of file   : ', i_fname)
        print('Write  to file        : ', o_fname)
        data = np.array([self.head_int],   self.type_i4); data.tofile(o_f)
        data = np.array([self.head_float], self.type_f8); data.tofile(o_f)
        # mean function and vertical axis
        mean = np.array(mean)     # remove possible mask
        y    = np.array(y); dy = np.diff(y)
        # blending parameters - change here
        delta_smooth = 9
        loc   = 15
        start = 564
        delta_smooth = 12
        loc   = 10
        start = 605
        # do not change
        bloc   = y[start+loc]
        dblend = delta_smooth*dy[-1]
        # plot blending
        plt.figure()
        plt.xlabel('y'); plt.ylabel('erf(y) - blending')
        plt.plot(y, ( - erf((y - bloc)/dblend) + 1) / 2)
        plt.axvline(y[592], color='red',  linestyle='dotted')
        plt.axvline(y[564], color='gray', linestyle='dotted')
        plt.axvline(y[620], color='gray', linestyle='dotted')
        plt.grid(True)
        plt.show
        # io blending 
        i_f.seek(self.headersize, 0) # jump over header
        nxy = self.nx*self.ny
        # buffered io
        for k in range(self.nz):
            # read input xy - plane
            i_dat = np.fromfile(i_f, self.type_f8, nxy)
            i_dat = i_dat.reshape((self.nx,self.ny), order='F')
            # blending
            for j in range(self.ny):
                i_dat[:,j] = i_dat[:,j]            * (    ( - erf((y[j] - bloc)/dblend) + 1) / 2) + \
                              mean[np.newaxis,j]   * (1 - ( - erf((y[j] - bloc)/dblend) + 1) / 2)
            # write new file
            i_dat = i_dat.reshape((nxy), order='F')
            i_dat[:].astype(self.type_f8).tofile(o_f)
        # close files
        i_f.close()
        o_f.close()
        print('Write  to file        : DONE')
        return
    #-----------------------------------------------------------------------#   
    def write_iniscal(self, mean, sc=1, da=0):
        # needs reference field for header information
        # input and output name
        i_fname = self.fname[0]
        o_fname = self.field_name(self.path, 'phi1', self.index+1)[0]
        # open input-/output-files
        i_f = open(i_fname,'rb') # input 
        o_f = open(o_fname,'wb') # output 
        # write header of ofile
        print('--------------------------------------------------')
        print('Write iniscal field   : ', o_fname)
        data = np.array([self.head_int],   self.type_i4); data.tofile(o_f)
        data = np.array([self.head_float], self.type_f8)
        data[0,-2] = sc; data[0,-1] = da; data.tofile(o_f)
        # write mean
        mean = np.array(mean)        # remove possible mask
        i_f.seek(self.headersize, 0) # jump over header
        nxy = self.nx*self.ny
        # buffered io
        i_dat = np.zeros((self.nx, self.ny))
        i_dat[:,:] = mean[np.newaxis, :]
        i_dat = i_dat.reshape((nxy), order='F')
        for k in range(self.nz):
            i_dat[:].astype(self.type_f8).tofile(o_f)
        i_f.close(); o_f.close()
        print('Write  to file        : DONE')
        return
    #-----------------------------------------------------------------------#   
    def read_3d_field(self):
        # skip header, read fields, reshape to 3d array and save
        print("Reading binary fields :   ...")
        dtype  = self.type_f8
        if self.var=='flow':   list_var = ['u', 'v', 'w']
        elif self.var=='scal': list_var = ['phi1', 'phi2']
        elif self.var=='vort': list_var = ['vortx', 'vorty', 'vortz']
        elif self.var[:3]=='eps': 
            list_var = ['eps']
            if self.var=='eps_int' or self.var=='eps_bit': dtype = self.type_i1
        else: list_var = [self.var]
        # read fields 
        for i in range(len(self.fname)):
            field = self.read_binary(self.fname[i], dtype, self.points, self.headersize)
            if self.var=='eps_bit':
                rsize = self.points; eps = np.zeros(rsize)
                field = int2bit_1(eps, field) # eps = int2bit_2(eps,data) # faster
            self.data[list_var[i]] = field.reshape((self.nx, self.ny, self.nz), order='F')
            del field
        print("Reading binary fields :   DONE")
        return
    #-----------------------------------------------------------------------#       
    def field_name(self, path, var, index):
        # ini
        num = []; fname = []
        # lists
        flow_list = ['u', 'rand_u', 'v', 'rand_v', 'w', 'rand_w', 'flow']
        scal_list = ['phi1', 'phi2', 'rand_phi1', 'rand_phi2', 'scal']
        geom_list = ['eps_bit', 'eps_int', 'eps_real']
        pres_list = ['pre']
        vort_list = ['vort', 'vortx', 'vorty', 'vortz']
        list_1    = ['u', 'rand_u', 'phi1', 'rand_phi1', 'eps_bit', 'eps_int', 'eps_real', 'pre', 'vortx']
        list_2    = ['v', 'rand_v', 'phi2', 'rand_phi2', 'vorty']
        list_3    = ['w', 'rand_w', 'vortz']
        # prefixes
        if   var in flow_list: prefix = 'flow.'
        elif var in scal_list: prefix = 'scal.'
        elif var in geom_list: prefix = 'eps'
        elif var in pres_list: prefix = 'Pressure'
        elif var in vort_list: prefix = 'VorticityVector'
        elif var[:4]=='rand':  prefix = prefix + 'rand'
        else: raise ValueError('Variable not implemented')
        # suffix
        if   var=='flow':   num = [1,2,3] 
        elif var=='scal':   num = [1,2] 
        elif var=='vort':   num = [1,2,3] 
        elif var in list_1: num = [1] 
        elif var in list_2: num = [2] 
        elif var in list_3: num = [3] 
        # build names
        for i in num:
            name = ''
            if not var in ['pre', 'vort', 'vortx', 'vorty', 'vortz']:
                name = prefix             
                if not var[:4]=='rand': name = name + str(index) 
                name = name + '.' + str(i)
            else:
                name = prefix + format(index, '0>6d') + '.' + str(i) 
                # name = prefix + format(index, '0>6d') 
                # if not var=='pre': name += '.' + str(i) 
            fname.extend([path + name])
        return fname
#############################################################################    
#---------------------------------------------------------------------------#
# read grid binary file / grid fortran record (Cedrick)
#
# Variables:
#   - path : path to file
#   - fname: 'grid'
#
#---------------------------------------------------------------------------#  
#############################################################################          
class DnsGrid: 
    def __init__(self, path, fname):
        # params
        self.path    = path
        self.name    = path + fname
        self.type_i4 = np.dtype('<i4') 
        self.type_f8 = np.dtype('<f8') 
        
        f_h = open(self.name, 'rb') 
    
        [n_isize, rec_isize] = self.read_fortran_record(f_h, self.type_i4) 
        [n_fsize, rec_fsize] = self.read_fortran_record(f_h, self.type_f8) 
        
        self.nx = rec_isize[0];  self.lx = rec_fsize[0]
        self.ny = rec_isize[1];  self.ly = rec_fsize[1]
        self.nz = rec_isize[2];  self.lz = rec_fsize[2]
    
        [nx_test, self.x] = self.read_fortran_record(f_h, self.type_f8) 
        [ny_test, self.y] = self.read_fortran_record(f_h, self.type_f8) 
        [nz_test, self.z] = self.read_fortran_record(f_h, self.type_f8)
        
        self.dx = self.x[1:] - self.x[:-1]
        self.dy = self.y[1:] - self.y[:-1]
        self.dz = self.z[1:] - self.z[:-1]
        
        print("--------------------------------------------------")       
        print("Reading grid file         ...")
        print('Domain size (Lx*Ly*Lz):   '+'('+str(round(self.lx,5))+' x '+str(round(self.ly,5))+' x '+str(round(self.lz,5))+')')
        print('Grid size   (nx*ny*nz):   '+'('+str(self.nx)+' x '+str(self.ny)+' x '+str(self.nz)+')')
        print('Grid dx     (min-max) :   '+'min='+str(round(self.dx.min(),5))+', max='+str(round(self.dx.max(),5)))
        print('Grid dy     (min-max) :   '+'min='+str(round(self.dy.min(),5))+', max='+str(round(self.dy.max(),5)))
        if self.z[-1] > 0.0:
            print('Grid dz     (min-max) :   '+'min='+str(round(self.dz.min(),5))+', max='+str(round(self.dz.max(),5)))

        if not ( nx_test == self.nx and ny_test == self.ny and nz_test == self.nz ): 
            print('Grid Dimensions do not match') 
            sys.exit()
        return
    #------------------------------------------------------------------------#        
    def read_fortran_record(self, f_h, t):      # reading step by step with checking
        t_i4  = np.dtype('<i4')                
        dum1, = np.fromfile(f_h, t_i4, count=1)  
        n     = dum1//t.itemsize                # number of entries to be read
        rec   = np.fromfile(f_h, t,    count=n)
        dum2  = np.fromfile(f_h, t_i4, count=1) # at the end again the number of read entries
        if dum1 != dum2:                        # checking 
            print('ERROR READING RECORD', dum1, dum2) 
            sys.exit()
        else:
            return [n, rec]
        
#############################################################################    
#---------------------------------------------------------------------------#
# read InputParams from dns.ini
#
# Variables:
#
#---------------------------------------------------------------------------#  
#############################################################################  
class InputParam:
    def __init__(self, path, pathout, fout, sim_type='ekman'):
        self.path  = path
        self.pathout = pathout
        self.fout  = fout
        self.name  = 'dns.ini'
        self.fname = self.path + '/' + self.name
        self.param = dict()
        self.type  = sim_type
        # read default params
        if os.path.exists(self.path):
            self.get_default_values()
            self.get_input_param(sim_type=self.type)
    #------------------------------------------------------------------------#        
    def get_default_values(self,):
        # set default values
        self.param['rho']   = 1
        self.param['f_cor'] = 1
        self.param['geo']   = 1
        self.param['Ro']    = 1
        return 
    #------------------------------------------------------------------------#            
    def get_value(self, list, topic, parameter):
        i = 0
        for p in list:
            if p == topic:
                break
            else: 
                i+=1
                continue
        for par,val in list[i+1:]:
            if par == parameter:
                params = [float(i) for i in (val.split(','))]
                break
        return params
    #------------------------------------------------------------------------#        
    # def get_entry(self,l,key):
    #     i = 0
    #     for p in l:
    #         if p == key:
    #             break
    #         else: 
    #             i+=1
    #             continue
    #     return i
    #------------------------------------------------------------------------#        
    def get_input_param(self,sim_type):
        f = open(self.fname, "r")
        l_param = []
        for line in f:
            line = line.split(" ")[0].strip()
            if len(line)==0:
                continue
            else:
                if line[0] == '[':
                    l_param.append(line[1:-1])
                else:
                    l_param.append(line.split('=', 2))
                                
        # get rotation values
        if sim_type=='ekman':
            rot_param = self.get_value(l_param, 'Rotation', 'Parameters')      
            self.param['u_geo'] =   np.cos(rot_param[0])*rot_param[1]
            self.param['w_geo'] = - np.sin(rot_param[0])*rot_param[1]

        # get dimensionless numbers
        if sim_type=='ekman':
            re = 'Re_ro'
        else:
            re = 'Re'
        self.param[re]   = self.get_value(l_param, 'Parameters', 'Reynolds')[0]        
        self.param['Pr'] = self.get_value(l_param, 'Parameters', 'Schmidt')[0]        

        # compute further values
        self.param['visc']      = self.param[re]**-1
        if sim_type=='ekman':
            self.param['lambda_ro'] = self.param['geo'] / self.param['f_cor']
            self.param['d_lam']     = np.sqrt(2*self.param['visc'] / self.param['f_cor'])
            self.param['Re_d']      = np.sqrt(2*self.param['Re_ro'])
        
        # save params to nc file
        write_dict2nc(self.param, self.pathout, self.fout)
        return
       
#############################################################################    
#---------------------------------------------------------------------------#
# read dns.out file 
#
# Variables:
#   - path  : path to file
#   - fname : 'dns.out'
#   - vlevel: verbosity level
#
#---------------------------------------------------------------------------#  
#############################################################################        
class DnsOut:
    def __init__(self, path, fname, vlevel=2): 
        # params
        self.path = path
        self.name = path + fname
        self.data = dict()
        print("--------------------------------------------------")       
        print("Reading dns.out file      ...")
        # open dns.out file and process data
        with open(self.name, 'r') as f:
            dns_out = f.readlines()
            dns_out = np.asanyarray(dns_out)
            #---------------------------------------------#
            # get part of header
            dns_header = dns_out[0][15:20]
            # clean up
            dns_out_clean = dns_out
            index_clean   = []
            for i in range(dns_out.size):
                if dns_out[i][15:20]==dns_header:
                    index_clean.append(i)
            j = 0
            index_clean = index_clean[::2]
            for i in index_clean:
                dns_out_clean = np.delete(dns_out_clean,range(i-j,i+3-j), axis=0)
                j = j + 3
            # get array size
            col_out = len(dns_out_clean[0].split())
            lin_out = len(dns_out_clean)
            # data array
            if vlevel < 2:
                out = np.zeros((lin_out,col_out-1)) # skip first column with zeros
                # split
                i = 0
                for line in dns_out_clean:
                    data_line = line.split()
                    for j in range(col_out-1):
                        out[i,j] = float(data_line[j+1])
                    i += 1
            else:
                out = np.zeros((lin_out,col_out-1)) # skip second column with zeros
                # split
                i = 0
                for line in dns_out_clean:
                    data_line = line.split()
                    out[i,0] = float(data_line[0][1:-1])
                    for j in range(1,col_out-1):
                        out[i,j] = float(data_line[j+1])
                    i += 1
        # stored variables
        list_print = []
        list_var   = []
        if vlevel >= 2: list_print.append('timestamp (date)'); list_var.append('date')        
        list_print.extend(['iteration_step (it)', 'time_total (time)','time_step (dt)', \
                           'cfl_number (clf)','dif_number (dif)','visc  (visc)',        \
                           'dil_min  (dil_min)','dil_max (dilmax)'])
        list_var.extend(['it', 'time','dt','cfl','dif','visc','dil_min','dil_max'])
        # print
        i = 0
        for var in list_print:
            print('Column ',str(i), '    Variable:  ', var )
            self.data[list_var[i]] = out[:,i]
            i += 1  
        del out, dns_out, dns_out_clean, f, line
        return
#############################################################################
#---------------------------------------------------------------------------#
# process avg.nc files
#
# Variables:
#   - path  : path to file(s)
#   - fname : type of statistic files ('avg', 'avg1s', 'avg2s')
#
#---------------------------------------------------------------------------#
#############################################################################
class Statistics:
    """
        Read and process avgxxxx.nc files for channel flow
    """
    #-----------------------------------------------------------------------#
    def __init__(self, path, fname='avg', index=0):
        self.path  = path 
        self.fname = fname # avg, avg1s, avg2s ... 
        self.index = index
        return
    #-----------------------------------------------------------------------#
    def list_var(self):
        # Lists with variables -- example
            # list_var = avg.list_var()
        # flow lists
            # list_dim           = ['t', 'y', 'it']
            # list_means         = ['rR', 'rU', 'rV', 'rW', 'rP', 'rT', 're', 'rh', 'rs', 'rB', 'fU', 'fV', 'fW', 'fT', 'fe', 'fh', 'fs']
            # list_fluctuations  = ['Tke', 'Rxx', 'Ryy', 'Rzz', 'Rxy', 'Rxz', 'Ryz', 'rP2', 'rR2', 'rT2', 'fT2', 're2', 'fe2', 'rh2', 'fh2', 'rs2', 'fs2']
            # list_vorticity     = ['Wx', 'Wy', 'Wz', 'Wx2', 'Wy2', 'Wz2']
            # list_rxxbudget     = ['Rxx_t', 'Bxx', 'Cxx', 'Pxx', 'Exx', 'PIxx', 'Fxx', 'Txxy_y', 'Txxy',  'Gxx', 'Dxx']
            # list_ryybudget     = ['Ryy_t', 'Byy', 'Cyy', 'Pyy', 'Eyy', 'PIyy', 'Fyy', 'Tyyy_y', 'Tyyy',  'Gyy', 'Dyy']
            # list_rzzbudget     = ['Rzz_t', 'Bzz', 'Czz', 'Pzz', 'Ezz', 'PIzz', 'Fzz', 'Tzzy_y', 'Tzzy',  'Gzz', 'Dzz']
            # list_rxybudget     = ['Rxy_t', 'Bxy', 'Cxy', 'Pxy', 'Exy', 'PIxy', 'Fxy', 'Txyy_y', 'Txyy',  'Gxy', 'Dxy']
            # list_rxzbudget     = ['Rxz_t', 'Bxz', 'Cxz', 'Pxz', 'Exz', 'PIxz', 'Fxz', 'Txzy_y', 'Txzy',  'Gxz', 'Dxz']
            # list_ryzbudget     = ['Ryz_t', 'Byz', 'Cyz', 'Pyz', 'Eyz', 'PIyz', 'Fyz', 'Tyzy_y', 'Tyzy',  'Gyz', 'Dyz']
            # list_tkebudget     = ['Tke_t', 'Buo', 'Con', 'Prd', 'Eps', 'Pi', 'Trp', 'Trp1', 'Trp2', 'Trp3', 'Trp1_y', 'Trp2_y', 'Trp3_y', 'G', 'D', 'Phi', 'UgradP']
            # list_higherorder   = ['rU3', 'rU4', 'rV3', 'rV4', 'rW3', 'rW4']
            # list_derivfluct    = [ ]; # list_rhobudget = [ ]; # list_stratified= [ ]
        # scalar 
            # list_sdim          = ['t', 'y', 'it']
            # list_smeans        = ['rS', 'fS', 'rS_y', 'fS_y', 'rQ', 'fQ']
            # list_sfluctuations = ['Rsu', 'Rsv', 'Rsw', 'fS2', 'fS3', 'fS4', 'rS2', 'rS3', 'rS4']
            # list_rssbudget     = ['Rss_t', 'Css', 'Pss', 'Ess', 'Tssy1', 'Tssy2', 'Tssy_y', 'Dss', 'Qss']
            # list_rsubudget     = [ ]; # list_rsvbudget    = [ ] # list_rswbudget = [ ]
            # list_sderivfluct   = [ ]; # list_crossscalars = [ ]
        return list(self.data.variables.keys())
    #-----------------------------------------------------------------------#
    def read_nc(self, path, name='avg', index=0):    
        self.file  = path + name + str(index) + '.nc'            
        self.data  = nc.Dataset(self.file, 'r', format='netCDF-4')
        print("--------------------------------------------------")       
        print('Reading file          :  ', self.file)
        return
    #-----------------------------------------------------------------------#
    def merge_nc_all(self, delete=True):
        # go to directory
        path_current = os.path.abspath(os.curdir)
        print("--------------------------------------------------")       
        print('Current path          :  ', path_current)
        print('Move to path          :  ', self.path)
        os.chdir(self.path)
        # check if file already exists
        file_merge = self.path+self.fname+str(self.index)+'.nc'
        if os.path.exists(file_merge):
            if not delete:
                print('ERROR File exists already, no merging possible')
                sys.exit()
            else:
                print('Delete file           :  ', file_merge)
                os.remove(file_merge)    
        # ordered list of nc-files for merging
        list_nc = [n[len(self.path+self.fname):-3] for n in glob.glob(self.path + self.fname + '*.nc')]
        if '_t0'  in list_nc: list_nc.remove('_t0')
        if '_in0' in list_nc: list_nc.remove('_in0')
        list_nc = list(map(int,list_nc)); list_nc.sort(); list_merge = ''
        for i in list_nc: list_merge += self.fname + str(i) + '.nc '
        # merge
        print('Merging all .nc fiels :  ', 'indices --- start:', list_nc[0], 'end:', list_nc[-1], 'step:', list_nc[1]-list_nc[0], 'num:', len(list_nc))
        command = 'ncrcat --4 -L 9 {} {}{}.nc'.format(list_merge,self.fname,self.index)
        print('Merging with command  :  ', command[:15], '... all files ...', command[-9:]) 
        os.system(command)
        # move back
        print('Move back to path     :  ', path_current)
        os.chdir(path_current)
        return
    #-----------------------------------------------------------------------
    def avg_all(self, fin='avg', fout='avg_t', skip=0, ensemble=True, twa=False, delete=True): 
        # name of new/destiny file
        fname_dst = self.path + fout + str(self.index) + '.nc'
        dtype     = 'f8'
        
        # check if file already exists
        file_dst = file_exist(self.path, fname_dst, delete)

        # read 
        self.read_nc(self.path, fin, self.index)

        # lists of all variables to be processed, remove dimensions
        list_var = self.list_var()
        for key in ['t', 'y','it']: list_var.remove(key)
                
        # read dimensions 
        t  = self.data.variables['t'][skip:]
        dt = np.diff(t)
        y  = self.data.variables['y'][:]
        it = self.data.variables['it'][skip:]

        # add dimensions
        file_dst.createDimension('y', y.size)
        file_dst.createDimension('t', t.size)
        
        # add netcdf variables
        y_dst  = file_dst.createVariable('y',  dtype, ('y',))
        t_dst  = file_dst.createVariable('t',  dtype, ('t',))
        it_dst = file_dst.createVariable('it', dtype, ('t',))
        
        # write/assign new variables
        y_dst[:]  = y[:]
        t_dst[:]  = t[:]
        it_dst[:] = it[:]

        # remove 1d-variables list
        list_var_1d = []
        for key in list_var:
            if len(self.data.variables[key].shape) < 2: 
                list_var_1d.append(key); list_var.remove(key)
            
        # averaging over iteration steps (ensemble) / time steps (twa)
        for key in list_var:
            if ensemble:
                var_dst    = file_dst.createVariable(key, dtype, ('y',))
                var_mean   = self.data.variables[key][skip:,:].mean(axis=0)
            elif twa:
                var_dst    = file_dst.createVariable(key, dtype, ('y',))
                var_all    = self.data.variables[key][skip:,:]
                area_sum   = np.zeros(y.size)
                for i in range(t.size-1):
                    area_sum[:] += dt[i] / 2 * (var_all[i,:] + var_all[i+1,:])
                var_mean   = area_sum / (t[-1] - t[0])  
            var_dst[:] = var_mean[:]
        # save temporal evolution of global parameters
        for key in list_var_1d:
            var_dst    = file_dst.createVariable(key, dtype, ('t',))
            var_tmp    = self.data.variables[key][skip:]
            var_dst[:] = var_tmp[:]

        # close new file 
        file_dst.close()
        return
#############################################################################
#---------------------------------------------------------------------------#
# AVG 3d-fields
#
# Variables:
#   - path : where the fields are located
#   - index: index of the avg
#
#---------------------------------------------------------------------------#
#############################################################################
class AVG_Fields:
    """
        Read and process 3d-fields for channel flow
    """
    #-----------------------------------------------------------------------#
    def __init__(self, path, ftype='flow', scalar=True, pressure=True, num_scal=1, skip=0):
        self.path     = path
        self.var      = []
        if len(ftype) > 0: self.var.extend([ftype])
        self.scalar   = scalar
        self.pressure = pressure
        self.num_scal = num_scal
        self.skip     = skip
        self.data     = dict()
        if scalar:
            self.var.extend(['phi1'])
            if self.num_scal > 1: self.var.extend(['phi2'])
        if pressure:
            self.var.extend(['pre'])     
        print("--------------------------------------------------")       
        print("Processing 3d-fields   : " + str(self.var) )
        print("In current path        : " + str(self.path))
        self.get_indices()
        return
    #---------------------------------------------------------------------------#
    def get_indices(self):
        ftype_end    = 2 # cut off .[1,2,3] data name suffix
        indices_list = [[] for i in range(len(self.var))]; k=0
        for ikey in self.var:
            ilist=[]
            command = 'find ' + self.path + ' -name '
            key = ikey; suf = '"' + key + '.*"' 
            if key in ['u','v','w']:
                key = 'flow'
                if   ikey=='u':    suf = '"flow.*.1"' 
                elif ikey=='v':    suf = '"flow.*.2"' 
                elif ikey=='w':    suf = '"flow.*.3"' 
            elif key in ['phi1','phi2']: 
                key = 'scal'
                if   ikey=='phi1': suf = '"scal.*.1"' 
                elif ikey=='phi2': suf = '"scal.*.2"' 
            elif key=='pre':        
                key = 'Pressure';  suf = '"Pressure*"' 
            command = command + suf
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE) 
            for file in p.stdout.readlines():
                dummy   = file.strip()
                cstart  = len('{}/{}'.format(self.path, key))
                if key == 'Pressure': cstart = cstart - 1
                cend    = len(dummy) - ftype_end
                if (is_number(dummy[cstart:cend])):
                    niter = int(dummy[cstart:cend])
                    ilist.append([niter])
            indices = np.sort(np.array(ilist, dtype="object")[:,0], axis=0)
            indices_list[k].extend([np.unique(indices)])
            k+=1
        merged_list = []
        for l in indices_list:
            merged_list += l
        self.indices = np.array(merged_list, dtype="object")
        self.indices = self.indices[:,self.skip:]
        for i in range(len(self.var)):
            print('With the indices       :',self.var[i],self.indices[i,:])
            if self.skip != 0:
                print('In total the first {} fields are skipped!'.format(self.skip))
        return
    #-----------------------------------------------------------------------#            
    def avg_x(self):
        if not self.var[0] == 'flow': raise ValueError('Not implemented')
        if not np.diff(self.indices, axis=0).sum() == 0: raise ValueError('Not implemented')
        print("--------------------------------------------------")       
        print('Avg fields in x       : ...')
        fnum   = self.indices.shape[1]
        header = Field(self.path, self.var[0], index=self.indices[0,0])
        ny = header.ny; nz = header.nz
        
        indices = self.indices[0,:]    
        
        var_list = ['u','v','w','uu','vv','ww','uv','uw','vw']
        vel_list = [ ['u','u'],['v','v'],['w','w'],['u','v'],['u','w'],['v','w'] ] 
        if self.scalar: 
            var_list.extend(['phi1','phi1phi1','uphi1','vphi1','wphi1'])
            sc1_list = [ ['phi1','phi1'],['u','phi1'],['v','phi1'],['w','phi1'] ]
            if self.num_scal > 1: 
                var_list.extend(['phi2','phi2phi2','uphi2','vphi2','wphi2'])
                sc2_list = [ ['phi2','phi2'],['u','phi2'],['v','phi2'],['w','phi2'] ]
        if self.pressure: 
            var_list.extend(['pre','prepre','upre','vpre','wpre'])
            pre_list = [ ['pre','pre'],['u','pre'],['v','pre'],['w','pre'] ]

        f2d = np.zeros((ny,nz,len(var_list)))
        for i in range(fnum):
            # read fields
            f = Field(self.path, var='flow', index=indices[i])
            f.read_3d_field(); k = 0
            f2d[:,:,k] = f2d[:,:,k] + f.data['u'].mean(axis=0); k += 1   
            f2d[:,:,k] = f2d[:,:,k] + f.data['v'].mean(axis=0); k += 1
            f2d[:,:,k] = f2d[:,:,k] + f.data['w'].mean(axis=0); k += 1   
            for id1,id2 in vel_list: 
                f2d[:,:,k] = f2d[:,:,k] + (f.data[id1]*f.data[id2]).mean(axis=0); k += 1   
            if self.scalar: 
                s1 = Field(self.path, var='phi1', index=indices[i])
                s1.read_3d_field()
                f.data.update(s1.data); del s1             
                f2d[:,:,k] = f2d[:,:,k] + f.data['phi1'].mean(axis=0); k += 1
                for id1,id2 in sc1_list: 
                    f2d[:,:,k] = f2d[:,:,k] + (f.data[id1]*f.data[id2]).mean(axis=0); k += 1   
                if self.num_scal > 1:   
                    s2 = Field(self.path, var='phi2', index=indices[i])
                    s2.read_3d_field()
                    f.data.update(s2.data); del s2
                    f2d[:,:,k] = f2d[:,:,k] + f.data['phi2'].mean(axis=0); k += 1
                    for id1,id2 in sc2_list: 
                        f2d[:,:,k] = f2d[:,:,k] + (f.data[id1]*f.data[id2]).mean(axis=0); k += 1   
            if self.pressure: 
                p = Field(self.path, var='pre', index=indices[i])
                p.read_3d_field()
                f.data.update(p.data); del p
                f2d[:,:,k] = f2d[:,:,k] + f.data['pre'].mean(axis=0); k += 1   
                for id1,id2 in pre_list: 
                    f2d[:,:,k] = f2d[:,:,k] + (f.data[id1]*f.data[id2]).mean(axis=0); k += 1   
        # avg
        k = 0
        for key in var_list:
            self.data[key] = f2d[:,:,k]; k += 1
        keylist = [key for key in self.data]
        for key in keylist:
            self.data[key] = self.data[key] / fnum
        return
    #-----------------------------------------------------------------------# 
    def avg_z(self):
        if not self.var[0] == 'flow': raise ValueError('Not implemented')
        if not np.diff(self.indices, axis=0).sum() == 0: raise ValueError('Not implemented')
        print("--------------------------------------------------")       
        print('Avg fields in z       : ...')
        fnum   = self.indices.shape[1]
        header = Field(self.path, self.var[0], index=self.indices[0,0])
        nx = header.nx; ny = header.ny
        
        indices = self.indices[0,:]    
        
        var_list = ['u','v','w','uu','vv','ww','uv','uw','vw']
        vel_list = [ ['u','u'],['v','v'],['w','w'],['u','v'],['u','w'],['v','w'] ] 
        if self.scalar: 
            var_list.extend(['phi1','phi1phi1','uphi1','vphi1','wphi1'])
            sc1_list = [ ['phi1','phi1'],['u','phi1'],['v','phi1'],['w','phi1'] ]
            if self.num_scal > 1: 
                var_list.extend(['phi2','phi2phi2','uphi2','vphi2','wphi2'])
                sc2_list = [ ['phi2','phi2'],['u','phi2'],['v','phi2'],['w','phi2'] ]
        if self.pressure: 
            var_list.extend(['pre','prepre','upre','vpre','wpre'])
            pre_list = [ ['pre','pre'],['u','pre'],['v','pre'],['w','pre'] ]

        f2d = np.zeros((nx,ny,len(var_list)))
        for i in range(fnum):
            # read fields
            f = Field(self.path, var='flow', index=indices[i])
            f.read_3d_field(); k = 0
            f2d[:,:,k] = f2d[:,:,k] + f.data['u'].mean(axis=2); k += 1   
            f2d[:,:,k] = f2d[:,:,k] + f.data['v'].mean(axis=2); k += 1
            f2d[:,:,k] = f2d[:,:,k] + f.data['w'].mean(axis=2); k += 1   
            for id1,id2 in vel_list: 
                f2d[:,:,k] = f2d[:,:,k] + (f.data[id1]*f.data[id2]).mean(axis=2); k += 1   
            if self.scalar: 
                s1 = Field(self.path, var='phi1', index=indices[i])
                s1.read_3d_field()
                f.data.update(s1.data); del s1             
                f2d[:,:,k] = f2d[:,:,k] + f.data['phi1'].mean(axis=2); k += 1
                for id1,id2 in sc1_list: 
                    f2d[:,:,k] = f2d[:,:,k] + (f.data[id1]*f.data[id2]).mean(axis=2); k += 1   
                if self.num_scal > 1:   
                    s2 = Field(self.path, var='phi2', index=indices[i])
                    s2.read_3d_field()
                    f.data.update(s2.data); del s2
                    f2d[:,:,k] = f2d[:,:,k] + f.data['phi2'].mean(axis=2); k += 1
                    for id1,id2 in sc2_list: 
                        f2d[:,:,k] = f2d[:,:,k] + (f.data[id1]*f.data[id2]).mean(axis=2); k += 1   
            if self.pressure: 
                p = Field(self.path, var='pre', index=indices[i])
                p.read_3d_field()
                f.data.update(p.data); del p
                f2d[:,:,k] = f2d[:,:,k] + f.data['pre'].mean(axis=2); k += 1   
                for id1,id2 in pre_list: 
                    f2d[:,:,k] = f2d[:,:,k] + (f.data[id1]*f.data[id2]).mean(axis=2); k += 1   
        # avg
        k = 0
        for key in var_list:
            self.data[key] = f2d[:,:,k]; k += 1
        keylist = [key for key in self.data]
        for key in keylist:
            self.data[key] = self.data[key] / fnum
        return
    #-----------------------------------------------------------------------# 
    # average across centerline
    def avg_xcl(self): 
        print("--------------------------------------------------")       
        print('avg centerline fields : ...')

        keylist = [key for key in self.data]

        # new half channel dimension
        ny_old = self.data[keylist[0]].shape[0]
        ny_new = int(np.ceil(ny_old/2))
        if   ny_old % 2 == 0: ny_new_up = ny_new   # even
        elif ny_old % 2 == 1: ny_new_up = ny_new-1 # odd
    
        # average mean quantities 
        vlist = ['u','v','w']
        if self.pressure:
            vlist.extend(['pre'])
        lsign = [ 1, -1,  1,  1 ]
    
        for i in range(len(vlist)):
            self.data[vlist[i]] = 0.5*( self.data[vlist[i]][:ny_new,:]  \
                                        + lsign[i]*np.flip(self.data[vlist[i]][ny_new_up:],axis=0))
        if self.scalar:
            slist = ['phi1','phi1phi1']
            if self.num_scal > 1:
                slist.extend(['phi2','phi2phi2',])
            for i in range(len(slist)):
                self.data[slist[i]] = 0.5*( self.data[slist[i]][:ny_new,:]  \
                                            + np.flip(1.-self.data[slist[i]][ny_new_up:],axis=0))
    
        list1 = ['uu','vv','ww','uw']
        list2 = ['uv','vw']
        if self.scalar:
            list1.extend(['vphi1'])
            # list1.extend(['phi1phi1','vphi1']) # not sure about this
            list2.extend(['uphi1','wphi1'])
            if self.num_scal > 1:
                list1.extend(['vphi2'])
                # list1.extend(['phi2phi2','vphi2']) # not sure about this
                list2.extend(['uphi2','wphi2'])
        if self.pressure:
            list1.extend(['prepre','upre','wpre'])
            list2.extend(['vpre'])
       
        for key in list1:
            self.data[key] = 0.5*( self.data[key][:ny_new,:]  \
                                       + lsign[i]*np.flip(self.data[key][ny_new_up:],axis=0)) 
        for key in list2:
            self.data[key] = 0.5*( self.data[key][:ny_new,:]  \
                                       - lsign[i]*np.flip(self.data[key][ny_new_up:],axis=0)) 
        return      
    #-----------------------------------------------------------------------#              
    def avg_phase(self, nphase): # nphase: number of phases in spanwise direction
        print("--------------------------------------------------")       
        print('avg phase fields      : ...')
        # new spanwise dimensions
        keylist = [key for key in self.data]
        nz_old  = self.data[keylist[0]].shape[1]
        nz_new  = nz_old // nphase
        if (nz_old % nphase) != 0:
            print("Number of nz not dividable by number of phases!")
            return
        for key in keylist:
            self.data[key] =  np.reshape(self.data[key],\
                                        (self.data[key].shape[0],nz_new,nphase), order='F').mean(axis=2)
        return    
    #-----------------------------------------------------------------------# 
    def avg_save_nc(self, name, delete=True, xstats=True, zstats=False):
        print("--------------------------------------------------")       
        print('save avg fields to .nc: ...')
        
        # name of new/destiny file
        fname_dst = self.path + name + '.nc'
        dtype     = 'f8'
        
        # check if file already exists
        file_dst = file_exist(self.path, fname_dst, delete)
        
        # dimensions
        keylist = [key for key in self.data]
        if xstats:
            ny = self.data[keylist[0]].shape[0]
            nz = self.data[keylist[0]].shape[1]
        elif zstats:
            nx = self.data[keylist[0]].shape[0]
            ny = self.data[keylist[0]].shape[1]
            
        it = self.indices[0,:]

        # add dimensions
        file_dst.createDimension('ny', ny)
        if xstats:   file_dst.createDimension('nz', nz)
        elif zstats: file_dst.createDimension('nx', nx)
        file_dst.createDimension('it', it.size)

        # add netcdf variables
        for key in keylist:
            if   xstats: var_dst = file_dst.createVariable(key, dtype, ('ny','nz'))
            elif zstats: var_dst = file_dst.createVariable(key, dtype, ('nx','ny'))
            var_dst[:] = self.data[key][:,:]
        # close
        file_dst.close()  
        return
#############################################################################
#---------------------------------------------------------------------------#
# Read and process spectra files
#
# Variables:
#   - path     : where the fields are located
#   - index    : index of the avg
#   - filetag  : spectra type ('xsp','zsp','rsp','pow','pha')
#   - variables: ('uu', 'vv', ...)

#
#---------------------------------------------------------------------------#
#############################################################################
class Spectra:
    """
        Read and process spectra files
    """
    #-----------------------------------------------------------------------#
    def __init__(self, path, index, filetag, variables, time_avg=False):
        self.path        = path
        self.index_start = index[0]
        self.index_stop  = index[1]
        self.index_step  = index[2]
        self.filetag     = filetag      # no list!
        self.variables   = variables[:] # needs to be a list!
        self.tavg        = time_avg
        #        
        print("--------------------------------------------------")       
        print("Processing spectra     : " + "      ...")
        print("In current path        : " + str(self.path))
        self.grid = DnsGrid(path,'grid') # read grid field
        print('--------------------------------------------------')

        # datatype
        sizeofdata = 4  # in bytes
        etype = "<"     # little-endian # etype = ">" # big-endian
        dtype = "f"     # floating number
        self.ftype = etype + dtype + str(sizeofdata)

        # file names
        if self.filetag == "xsp" or self.filetag == "zsp" or self.filetag == "rsp" or self.filetag == "pow" or self.filetag == "pha" :
            self.nametag = "E"
        if self.filetag == "xcr" or self.filetag == "zcr" or self.filetag == "rcr":
            self.nametag = "C"
        # prerequisites for reading data 
        if   ( self.filetag == 'xsp' ):     # Number of modes
            self.nk = int(self.grid.nx / 2.)
        elif ( self.filetag == 'zsp' ):
            self.nk = int(self.grid.nz / 2.)
        elif ( self.filetag == 'rsp' ):
            if self.grid.nx != self.grid.nz: 
                raise ValueError('Grid in Horizontal needs to be uniform for radial spectra!')
            else:   
                self.nk = int(self.grid.nx / 2.)
        elif ( self.filetag == 'pow' or self.filetag == 'pha' ):
            self.nk  = int(self.grid.nx*self.grid.nz / 2.)
            self.nkx = int(self.grid.nx / 2.)
            self.nkz = int(self.grid.nz / 2.)
        print('Size of each file     :   ' + str(round(self.grid.ny*self.nk*sizeofdata*10**-6, 3)) + ' - Megabytes')
        return
    #---------------------------------------------------------------------------#
    def read_spectra(self): # reading and averaging spectra files
        num = (self.index_stop - self.index_start) / self.index_step + 1  # number of fields 
        print('Processing num steps  :   ' + str(int(num)))
        if ( self.filetag == 'pow' or self.filetag == 'pha' ):
            self.sp = np.zeros((len(self.variables),self.grid.nx//2,self.grid.ny,self.grid.nz)) # read only half of kz
        else:
            self.sp = np.zeros((len(self.variables),self.grid.ny,self.nk))
        for j in range(self.index_start, self.index_stop+self.index_step, self.index_step):
            i = 0
            a = np.zeros((self.grid.ny*self.nk))
            if ( self.filetag == 'pow' or self.filetag == 'pha' ):
                b = np.zeros((len(self.variables),self.grid.nx//2,self.grid.ny,self.grid.nz)) # read only half of kz
            else:
                b = np.zeros((len(self.variables),self.grid.ny,self.nk))
            for var in self.variables:
                if self.tavg:
                    file = self.path + self.filetag + str(self.index_start) + '-' + str(self.index_stop) + '.' + self.nametag + var
                else:
                    file = self.path + self.filetag + str(j) + '.' + self.nametag + var
                print('Processing file       :   ' + file)
                f = open(file,'rb')
                a = np.fromfile(f, np.dtype(self.ftype), self.grid.ny*self.nk)
                if ( self.filetag == 'pow' or self.filetag == 'pha' ):
                    b[i,:,:,:] = a.reshape((self.grid.nx//2,self.grid.ny,self.grid.nz), order='F')[:,:,:]
                    print('Warning - check multiplication factor for 2d spectra in x!!')
                else:
                    b[i,:,:] = a.reshape((self.grid.ny,self.nk)) 
                f.close()
                i += 1
            self.sp = self.sp + b
            if self.tavg: break
        del a, b                                     
        # postprocessing
        if not self.tavg: self.sp = self.sp / num # averaging
        # Files contain only half of the spectra
        if not ( self.filetag == 'pow' or self.filetag == 'pha' ): self.sp = self.sp * 2.
        # wavenumbers/-lengths / mode numbers
        if ( self.filetag == 'pow' or self.filetag == 'pha' ):
            self.kx = np.linspace(0, self.nkx-1, num=self.nkx) # mode number (k.size==nk)
            self.kz = np.linspace(0, self.nkz-1, num=self.nkz)  
            self.lx = np.flipud(self.kx[1:]**-1)               # l[nk] would be infinity; we skip it
            self.lz = np.flipud(self.kz[1:]**-1)                    
        else:
            self.k = np.linspace(0, self.nk-1, num=self.nk)
            self.l = np.flipud(self.k[1:]**-1)              
        print('--------------------------------------------------')
        return
    #-----------------------------------------------------------------------# 
    def spectra_save_nc(self, data, delete=False):
        print("--------------------------------------------------")       
        print('save spec avgs to .nc : ...')
        
        # name of new/destiny file
        fname_dst = self.path + self.filetag + '_avg_spectra' + '.nc'
        
        # check if file already exists
        file_dst = file_exist(self.path, fname_dst, delete)

        # add dimensions
        file_dst.createDimension('ny', self.grid.ny)
        file_dst.createDimension('nl', self.l.size)
        # add netcdf variables
        var_dst = file_dst.createVariable('l', self.ftype, ('nl')); var_dst[:] = self.l

        if ( self.filetag == 'pow' or self.filetag == 'pha' ):
            # add dimensions
            file_dst.createDimension('nkx', self.nkx)      
            file_dst.createDimension('nkz', self.nkz)      
            # add netcdf variables
            var_dst = file_dst.createVariable('kx', self.ftype, ('nkx')); var_dst[:] = self.kx
            var_dst = file_dst.createVariable('kz', self.ftype, ('nkz')); var_dst[:] = self.kz
            i = 0
            for key in self.variables:
                var_dst    = file_dst.createVariable(key, self.ftype, ('nkx','ny','nkz'))
                var_dst[:] = data[i,:,:,:]
                i += 1
        else:
            # add dimensions
            file_dst.createDimension('nk', self.nk)
            # add netcdf variables
            var_dst = file_dst.createVariable('k', self.ftype, ('nk')); var_dst[:] = self.k
            i = 0
            for key in self.variables:
                var_dst    = file_dst.createVariable(key, self.ftype, ('ny','nk'))
                var_dst[:] = data[i,:,:]
                i += 1

        print('Writing .nc file      :  ', str(fname_dst))
                  
        # close new file 
        file_dst.close()
        return
#############################################################################
#---------------------------------------------------------------------------#
# Read and process plane files
#
#   - path      : where the fields are located
#   - index     : index of the planes
#   - planetype : ('yz','yz','xz')
#   - planetype : number of planes in one plane file (levels)
#   - scalar    : planes of salar fields?
#
#---------------------------------------------------------------------------#
#############################################################################
class Planes:
    """
        Read and process planes files
    """
    #-----------------------------------------------------------------------#
    def __init__(self, path='', index=0, planetype='yz', num_pl=1, num_scal=1, pl_log=False):
        self.path        = path
        self.index_start = index[0]
        self.index_stop  = index[1]
        self.index_step  = index[2]
        self.indices     = [j for j in range(index[0], index[1]+index[2], index[2])]
        self.plane_type  = planetype
        self.num_pl      = num_pl 
        self.num_scal    = num_scal 
        self.pl_log      = pl_log 
        self.data        = dict()
        #  
        print("--------------------------------------------------")       
        print("Processing planes     : " + "      ...")
        print("In current path       : " + str(self.path))
        self.grid = DnsGrid(path,'grid') # read grid field
        print('--------------------------------------------------')
        # datatype
        sizeofdata = 4  # in bytes
        etype = "<"     # little-endian # etype = ">" # big-endian
        dtype = "f"     # floating number
        self.ftype = etype + dtype + str(sizeofdata)
        # prerequisites for reading data 
        self.fnames = []
        # order [u, v, w, phi1, .., p, log(entstrophy), log(grad(phi1)), ...]
        num_f = 4 
        if pl_log: 
            num_f += 1
            if num_scal > 0: 
                num_f += 2*num_scal
        else:
            if num_scal > 0: 
                num_f += num_scal
            
        name_tag = 'planes'
        if   self.plane_type == 'xy':
            for ind in self.indices:  self.fnames.extend([name_tag + 'K.' + str(ind)])
            self.nx1 = self.grid.nx; self.nx2 = self.grid.ny
            self.x1  = self.grid.x;  self.x2  = self.grid.y
            self.read = (self.nx1, self.nx2,    self.num_pl, num_f )
        elif self.plane_type == 'xz':
            for ind in self.indices:  self.fnames.extend([name_tag + 'J.' + str(ind)])
            self.nx1 = self.grid.nx; self.nx2 = self.grid.nz
            self.x1  = self.grid.x; self.x2  = self.grid.z
            self.read  = (self.nx1, self.num_pl, num_f,  self.nx2  )
        elif self.plane_type == 'yz':
            for ind in self.indices:  self.fnames.extend([name_tag + 'I.' + str(ind)])
            self.nx1 = self.grid.ny; self.nx2 = self.grid.nz
            self.x1  = self.grid.y;  self.x2  = self.grid.z
            self.read  = (self.nx1, self.num_pl, num_f,  self.nx2  )
        self.fsize = num_f * self.nx1 * self.nx2 * self.num_pl
        print('Size of each file     :   ' + str(round(self.fsize*sizeofdata*10**-6, 3)) + ' - Megabytes')
        return
    #---------------------------------------------------------------------------#
    def read_planes(self): # reading and averaging plane files
        num = (self.index_stop - self.index_start) / self.index_step + 1  # number of fields 
        print('Processing num steps  :   ' + str(int(num)))
        self.list_var = ['u', 'v', 'w']#, 'phi', 'p']
        if self.num_scal == 0:
            self.list_var.append('p')
            if self.pl_log: self.list_var.append('log_ent')
        elif self.num_scal > 0:
            for i in range(self.num_scal):
                tag = 'phi' + str(i+1)
                self.list_var.append(tag) 
            self.list_var.append('p')
            if self.pl_log: 
                self.list_var.append('log_ent')
                if self.num_scal > 0:
                    for i in range(self.num_scal):
                        tag = 'log_' +'phi' + str(i+1)
                        self.list_var.append(tag) 

        b = np.zeros(self.read)

        for l in range(len(self.indices)):
            file = self.path + self.fnames[l]
            print('Processing file       :   ' + file)
            f = open(file, 'rb')
            a = np.fromfile(f, self.ftype, self.fsize)
            a = a.reshape(self.read, order='F')
            f.close()
            b = b + a
        # averaging and assigning data to data.dict
        for i in range(len(self.list_var)):
            self.data[self.list_var[i]] = np.zeros((self.num_pl, self.nx1, self.nx2))
            for l in range(self.num_pl):
                if self.plane_type == 'xy':
                    self.data[self.list_var[i]][l,:,:] = b[:,:,l,i] / num               
                else:
                    self.data[self.list_var[i]][l,:,:] = b[:,l,i,:] / num
        del a,b 
        print('--------------------------------------------------')
        return
    #-----------------------------------------------------------------------# 
    def planes_save_nc(self, delete=True):
        print("--------------------------------------------------")       
        print('save planes avgs to .nc : ...')
        
        # name of new/destiny file
        if   self.plane_type == 'xy': ftag = 'planesK'
        elif self.plane_type == 'xz': ftag = 'planesJ'
        elif self.plane_type == 'yz': ftag = 'planesI'
        self.fname_dst = self.path + ftag + '_avg.nc'
        
        # check if file already exists
        file_dst = file_exist(self.path, self.fname_dst, delete)
        
        # add dimensions
        file_dst.createDimension(self.plane_type[0], self.nx1)
        file_dst.createDimension(self.plane_type[1], self.nx2)
        file_dst.createDimension('t', len(self.indices))
        file_dst.createDimension('num', self.num_pl)
        
        # add netcdf variables
        x1_dst = file_dst.createVariable(self.plane_type[0], self.ftype, (self.plane_type[0],))
        x2_dst = file_dst.createVariable(self.plane_type[1], self.ftype, (self.plane_type[1],))
        it_dst = file_dst.createVariable('it', self.ftype, ('t',))
        
        # write/assign new variables
        x1_dst[:] = self.x1[:]
        x2_dst[:] = self.x2[:]
        it_dst[:] = self.indices[:]
        
        # add netcdf variables
        for key in self.list_var:
            var_dst    = file_dst.createVariable(key, self.ftype, ('num',self.plane_type[0],self.plane_type[1]))
            var_dst[:] = self.data[key][:,:,:]

        print('Writing .nc file      :  ', str(self.fname_dst))
                  
        # close new file 
        file_dst.close()  
        return   
#############################################################################
#---------------------------------------------------------------------------#
# computing Anisotropy Invariance Maps (AIM) from avg.nc files
#
# Variables:
#   - path  : path to file(s)
#   - fname : type of statistic files, needs to be averaged in time 
#             (eg. 'avg', 'avg1s', 'avg2s', ....)
#
#---------------------------------------------------------------------------#
#############################################################################
class Invariance:
    #-----------------------------------------------------------------------#
    def __init__(self, path, fname='avg'):
        self.path  = path 
        self.fname = fname 
        return
    #-----------------------------------------------------------------------#
    def read_nc(self, path, name='avg'):    
        self.file  = path + name + '.nc'            
        self.d  = nc.Dataset(self.file, 'r', format='netCDF-4')
        print("--------------------------------------------------")       
        print('Reading file          :  ', self.file)
        # print(self.data.history) # print(self.data.dimensions)
        return       
    #-----------------------------------------------------------------------#
    def ani(self):
        """
            computes the deviatoric/anisotropic part of Re-stress tensor
        """
        
        # read avg.nc file from tlab (already avg in time)
        self.read_nc(self.path, self.fname)
        
        # skip always y=0 (avoid singularity)
        skip = 1
        
        # vertical coordinate 
        y = self.d['y'][skip:]; ys = y.size
        
        # Reynolds stresses (with the order ['uu','vv','ww','uv','uw','vw'])
        list_rs = ['Rxx', 'Ryy', 'Rzz', 'Rxy', 'Rxz', 'Ryz'] # storage order
        rs = np.zeros((ys,6))
        for i in range(6):
            rs[:,i] = self.d[list_rs[i]][skip:]

        # tke ([k = 1/2*u_ii] (half of the trace))
        # k = self.d['Tke'][skip:]
        k = 0.5 * rs[:,:3].sum(axis=1) # second option

        # normalized deviatoric/anisotropic part of the Re-stress tensor 
        #   [b_ij = (u_ij - 2/3*k*delta_ij) / (2*k)]
        b_ij = np.zeros((rs.shape))
        for i in range(6):
            b_ij[:,i] = rs[:,i] / (2*k[:])
        for i in range(3):
            b_ij[:,i] = b_ij[:,i] - 1/3
            
        return y, rs, b_ij
    #-----------------------------------------------------------------------#
    def aim(self):
        """
            computes anisotropic invariant maps of Re-stress tensor
        """
        
        # read avg.nc file (avg in time, without y=0)
        y, rs, bij = self.ani(); ys = y.size
        
        # build dev. Re-stress tensor
        r = np.zeros((3,3,ys))

        r[0,0,:] = bij[:,0] # uu
        r[1,0,:] = bij[:,3] # uv
        r[2,0,:] = bij[:,4] # uw

        r[0,1,:] = bij[:,3] # uv
        r[1,1,:] = bij[:,1] # vv
        r[2,1,:] = bij[:,5] # vw

        r[0,2,:] = bij[:,4] # uw
        r[1,2,:] = bij[:,5] # vw
        r[2,2,:] = bij[:,2] # ww

        # compute invariants of r tensor
        inv = np.zeros((ys,3)); c   = np.zeros((ys,3))
        xb  = np.zeros(ys);     yb  = np.zeros(ys)
        eta = np.zeros(ys);     xi  = np.zeros(ys)
        eig = np.zeros((3,ys))

        for i in range(ys):
            eig[:,i] = np.sort(LA.eig(r[:,:,i])[0])[::-1]
            # inv[i,0] = np.trace(r[:,:,i]) 
            inv[i,0] = eig[:,i].sum() # second option
            # ------------------------------
            # prefactor 2 is necessary here
            # inv[i,1] = 0.5*( np.trace(r[:,:,i])**2 - np.trace(r[:,:,i]**2)  )            # wrong !
            # inv[i,1] = -2 * (eig[0,i]*eig[1,i] + eig[0,i]*eig[2,i] + eig[1,i]*eig[2,i] ) # second option
            inv[i,1] = 2 * (eig[0,i]**2 + eig[0,i]*eig[1,i] + eig[1,i]**2)                 # third option
            # ------------------------------
            # prefactor 3 is necessary here
            inv[i,2] = 3 * LA.det(r[:,:,i])
            # inv[i,2] = 3* eig[:,i].prod()                            # second option
            # inv[i,2] = - 3 * eig[0,i]*eig[1,i]*(eig[0,i] + eig[1,i]) # third option
            # ------------------------------
            # pope (p.395, eq. 11.31f)
            # can be problematic if 3rd invariant is negativ (happens for ekman flow)!
            eta[i] =   (1/3) * (eig[0,i]**2 + eig[0,i]*eig[1,i] + eig[1,i]**2)
            eta[i] = eta[i] ** (1/2)
            
            xi[i]  = - (1/2) * (eig[0,i]*eig[1,i]*(eig[0,i] + eig[1,i]))
            xi[i]  = xi[i] ** (1/3) 
            # ------------------------------
            # Stiperski, Calaf et al. 2017
            c[i,0] =      eig[0,i] - eig[1,i]
            c[i,1] = 2 * (eig[1,i] - eig[2,i])
            c[i,2] = 3 *  eig[2,i] + 1
            
            xb[i] = c[i,0] + c[i,2]*(1/2)
            yb[i] = c[i,2]*(np.sqrt(3)/2)
        
        return inv , eta,xi  , xb,yb 
#############################################################################    
#---------------------------------------------------------------------------#
# My own miscellaneous functions
#
# Variables:
#   - 
#  
#---------------------------------------------------------------------------#  
#############################################################################
def filter_erf_tlab_1d(a, li=1., p1=3., p2=3.):
    """Spectral-erf filter function in tlab.

    Parameters
    ----------
    a:  array-like
        spectra to be filtered
    li: scale of domain (lx/lz)
    p1: float
        physical frequency for transition (> 0: High-pass, < 0: Low-pass)
    p2: float
        width of transition in log wavenumber space
    {
     p3: missing here, since it is 1D filter and 
         is only used to scale wavenumers if lx!=lz
         }
         
    Returns
    -------
    a: array-like
       filtered spectra
    """
    if p1 >= 0:
        sign_pass =  1.   # HIGHPASS
        off_pass  =  0.
    else:
        sign_pass = -1.   # LOWPASS
        off_pass  =  1. 
    n_ny = a.size # n_ny = n//2+1 since a is alredy kx/2
    # filtering
    print('logarithmic filter ratio = ', np.log(abs(p1))/np.log(n_ny*2/li))
    for i in range(1,n_ny):
        f = i / li
        damp = ( erf( np.log(f/abs(p1)) / p2) + 1.0) / 2.0
        damp = off_pass + sign_pass*damp
        a[i] = a[i] * damp
    return a
#-----------------------------------------------------------------------# 
def set_size(width, fraction=1, subplots=(1, 1), ratio='golden'):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 455
    elif width == 'beamer':
        width_pt = 300
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27
    
    if ratio == 'golden':
        # Golden ratio to set aesthetic figure height
        # https://disq.us/p/2940ij3
        ratio_fac = (5**.5 - 1) / 2
    else:
        ratio_fac = ratio

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * ratio_fac * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)
#-----------------------------------------------------------------------#
# interpolate 2d figure
def interpol_2d(data, idir_old, jdir_old, ipoli, ipolj):
    '''
    interpolates 2d-field on uniform grid for plotting,
    change interpolation type with 'intscheme'

    Parameters
    ----------
    data:      2d-field to be interpolated
    idir_old:  nodes distribution in i
    jdir_old:  nodes distribution in j
    ipoli:     number of nodes in i
    ipolj:     number of nodes in j
    '''
    from scipy.interpolate import interp1d
    intscheme = 'cubic' 

    ni_old = data.shape[0]
    # nj_old = data.shape[1] 

    ni_new = ipoli # ni_old * ipoli
    nj_new = ipolj # nj_old * ipolj
    
    idir_new = np.linspace(idir_old[0], idir_old[-1], ni_new, endpoint=True)
    jdir_new = np.linspace(jdir_old[0], jdir_old[-1], nj_new, endpoint=True)
    #------------------------------------------------------------------------------#
    # interpolating in j-direction
    int1 = np.zeros((ni_old, nj_new))
    print("interpolating in j-direction")
    for i in range(ni_old):
        cs = interp1d(jdir_old, data[i,:], kind=intscheme, fill_value='extrapolate')
        int1[i,:] = cs(jdir_new)
    # write interpolated data to normal data 
    int2 = int1; del int1
    print("interpolating in j-direction: done")
    #------------------------------------------------------------------------------#
    # interpolating in i-direction
    int1 = np.zeros((ni_new, nj_new))
    print("interpolating in i-direction")
    for j in range(nj_new):
        cs = interp1d(idir_old, int2[:,j], kind=intscheme, fill_value='extrapolate')
        int1[:,j] = cs(idir_new)
    # write interpolated data to normal data 
    data_int = int1; del int1, int2
    print("interpolating in i-direction: done")
    return data_int, idir_new, jdir_new
#-----------------------------------------------------------------------# 
# average across centerline
def avg_cl(data, even=True, scalar=False): 
    '''
    averages 1d-data over the center line

    Parameters
    ----------
    data:      1d-array to be averaged
    even:      even/odd variable
    scalar:    1d-mean scalar
    '''
    # new half channel dimension
    nx_old = data.shape[0]
    nx_new = int(np.ceil(nx_old/2))   # upper boundary
    cl     = np.zeros(nx_new)

    # averaging
    if even:
        if   nx_old % 2 == 0: # even
            cl = 0.5 * (data[:nx_new] + np.flipud(data[nx_new  :]))
        elif nx_old % 2 == 1: # odd
            cl = 0.5 * (data[:nx_new] + np.flipud(data[nx_new-1:]))
    else:
        if   nx_old % 2 == 0: # even
            cl = 0.5 * (data[:nx_new] - np.flipud(data[nx_new  :]))
        elif nx_old % 2 == 1: # odd
            cl = 0.5 * (data[:nx_new] - np.flipud(data[nx_new-1:])) 
    if scalar:
        if   nx_old % 2 == 0: # even
            cl = 0.5 * (data[:nx_new] + np.flipud(data[0] - data[nx_new  :]))
        elif nx_old % 2 == 1: # odd
            cl = 0.5 * (data[:nx_new] + np.flipud(data[0] - data[nx_new-1:]))

    return cl
#-----------------------------------------------------------------------# 
# average phase
def avg_phase(data, nphase): 
    '''
    averages 2d-data for phase in spanwise direction

    Parameters
    ----------
    data:      2d-array to be averaged
    nphase:    number of phases in spanwise direction
    '''
    # new spanwise dimensions
    ny      = data.shape[0]
    nz_old  = data.shape[1]
    nz_new  = nz_old // nphase\
        
    phase = np.zeros((ny, nz_new))
    
    if (nz_old % nphase) != 0:
        print("Number of nz not dividable by number of phases!")
        return
    
    phase =  np.reshape(data,(data.shape[0],nz_new,nphase), order='F').mean(axis=2)
    
    return phase 
#---------------------------------------------------------------------------#
def ubulk(u,y,ucl=1):
    '''
    computes the bulk velocity of the flow field u

    Parameters
    ----------
    u:     can be either a 3d-field or u=f(y)
    y:     vertical coordinate
    ucl:   centerline velocity
    '''
    if u.ndim==3:
        um = u.mean(axis=(0,2)) # horizontal average
    elif u.ndim==1:
        if u.size != y.size: # half channel
            if   y.size % 2 == 0: # even
                um = np.concatenate((u,np.flipud(u)))
            elif y.size % 2 == 1: # odd
                um = np.concatenate((u,np.flipud(u)[1:]))
        else:
            um = u
    else:
        print('wrong dimension of u-field')
        sys.exit()
    
    ub = (1 / y.max()) * integrate.simpson(um, y)
    print('--------------------------------------------------')
    print('bulk velocity         :  ', ub)

    # analytical solution
    ycl      = y.max() / 2                            # centerline position
    u_par    = - (ucl / ycl**2 ) * (y - ycl)**2 + ucl # parabolic ini velocity profile
    ub_exact = (2/3) * ucl                            # exact bulk velocity
    print('bulk velocity (exact) :  ', ub_exact)
    print('error [%]             :  ', (ub - ub_exact)/ub_exact * 100)
    
    return ub, u_par
#---------------------------------------------------------------------------#
def read_txt(path,name,header=0,footer=0):
    '''
    reads file.txt to numpy array

    Parameters
    ----------

    '''
    file = path + name
    print('--------------------------------------------------')
    print('read txt file to array:  ', file)
    with open(file, 'r') as f:
        raw = f.readlines()
        raw = np.asanyarray(raw)
        #---------------------------------------------#
        # clean up
        skiprows = header 
        raw = np.delete(raw, np.s_[:skiprows], axis=0)
        if footer > 0:
            skiprows = footer-header-1 
            raw = np.delete(raw, np.s_[skiprows:], axis=0)  
        # split and get array size
        i = 0
        col_out = len(raw[0].split())
        lin_out = len(raw)
        data    = np.zeros((lin_out,col_out))
        for line in raw:
            data_line = line.split()
            for j in range(col_out):
                data[i,j] = float(data_line[j])
            i += 1
    print('data array size       :  ', int(lin_out), 'x', int(col_out))

    return data
#---------------------------------------------------------------------------#
def read_binary_field(path, name, index, numf=1, grid=[1,1,1], prec='single', headsize=0):
    '''
    reads binary file to numpy array

    Parameters
    ----------

    '''
    # file 
    file = path + name + format(index, '0>6d')
    fname = []
    if numf>1: 
        for i in range(1,numf+1):
            fname.extend([file + '.' + str(i)])
    else:
        fname = [file]
    a = np.zeros((numf,grid[0],grid[1],grid[2]))  

    # data types
    if prec == 'single':
        dtype = np.dtype('<f4')
    if prec == 'double':
        dtype = np.dtype('<f8') 

    # load
    print('--------------------------------------------------')
    print('Reading binary file: ', numf, ' x ', grid)
    for i in range(numf):
        a[i,:,:,:] = Field.read_binary(1, fname[i], dtype, np.prod(grid), headsize).reshape((grid[0],grid[1],grid[2]), order='F')
        print('file: ', fname[i])
    return a
#---------------------------------------------------------------------------#
def file_exist(path, fname_dst, delete):
    '''
    checks if file exists already
    
    Parameters
    ----------

    '''
    if os.path.exists(fname_dst):
        print("--------------------------------------------------")
        print('ERROR File exists already')
        if delete:
            print('Delete file           :  ', fname_dst)
            os.remove(fname_dst)
        else: 
            print('STOP')
            sys.exit()
    if delete or not os.path.exists(fname_dst):
        print("--------------------------------------------------")
        print('Create new file       :  ', fname_dst)
        file_dst = nc.Dataset(fname_dst, 'w', format='NETCDF4')
    
    return file_dst
#---------------------------------------------------------------------------#
def write_dict2nc(dict, path, name):
    '''
    writes dict with values to nc file
    
    Parameters
    ----------

    '''
    # name of new/destiny file
    fname_dst = path + name + '.nc'
    dtype     = 'f8'
    
    # check if file already exists
    file_dst = file_exist(path, fname_dst, delete=True)
    
    # write dict to nc file
    i = 0
    file_dst.createDimension('val', 1)
    for varname in dict.keys():
        vardata  = dict[varname]
        array = isinstance(vardata, (np.ndarray)) 
        if not array:
            var_name = file_dst.createVariable(varname, dtype, ('val',))
        else:
            siz = vardata[:].size
            name = 't'+str(i); i+=1
            file_dst.createDimension(name,siz)
            var_name = file_dst.createVariable(varname, dtype, (name,))
        var_name[:] = vardata
    file_dst.close()
    return
#---------------------------------------------------------------------------#
# def error(a1, a2):    
#     delta = a1 - a2
#     mean  = delta.mean()
#     er    = ((delta-mean)**2).sum()/a1.size
#     print(er)
#     return er
#---------------------------------------------------------------------------#  
#---------------------------------------------------------------------------#      
# functions for Thomas algorithm
#---------------------------------------------------------------------------#  
def tdma(a,b,c,rhs): # own implementtation
    p, q = forward(a,b,c,rhs) 
    x    = backward(p,q)
    return x
#---------------------------------------------------------------------------#
def forward(a,b,c,r):
    p = np.zeros(a.size); q = np.zeros(a.size)
    p[0] = - c[0] / b[0]
    q[0] =   r[0] / b[0]
    for i in range(1,a.size):
        denom = b[i] + a[i]*p[i-1]
        p[i]  = - c[i] / denom
        q[i]  = (r[i] - a[i]*q[i-1]) / denom
    return p, q
#---------------------------------------------------------------------------#      
def backward(p,q):
    x = np.zeros(q.size)
    x[-1] = q[-1]
    for i in range(q.size-2,-1,-1):
        x[i] = q[i] + p[i]*x[i+1]
    return x
#---------------------------------------------------------------------------#      
# functions to convert/modify eps-field (IBM geometry) [bitwise <--> int1]
#---------------------------------------------------------------------------#      
def int2bit_1(out,data): # option 1
    bsize = data.size
    for i in range(bsize):
        ip = i * 8
        binary = f'{data[i]:08b}'
        if binary[0] == '-':
            binary = f'{data[i]+256:08b}'
        j = 0
        for k in range(-1,-9,-1):
            out[j+ip] = int(binary[k])
            j += 1
    return out
#---------------------------------------------------------------------------#      
def int2bit_2(out,data): # option 2 (bit faster then option 1)
    import struct
    bsize = data.size
    for i in range(bsize):
        ip = i * 8
        by   = struct.pack('b',data[i])
        by2b = ''.join(format(ord(by), '08b') for byte in by)
        j = 0
        for k in range(-1,-9,-1):
            out[j+ip] = int(str(by2b)[k])
            j += 1
    return out
#---------------------------------------------------------------------------#      
def bit2int_2(out,data):
    bsize = out.size
    for i in range(bsize):
        by = []
        ip = i * 8
        j  = 0
        for k in range(7,-1,-1):
            by.append(int(data[k+ip]))
            j += 1
        bit = "".join(str(i) for i in by)
        if bit[0] == '1':
            out[i] = int(bit,2)-256
        else:    
            out[i] = int(bit,2)
    return out
#############################################################################           
#---------------------------------------------------------------------------#
# Functions from Cedrick Ansorge
#
# Variables:
#   - 
#  
#---------------------------------------------------------------------------#   
#############################################################################        
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False    
#---------------------------------------------------------------------------#
def rotate(a,b,alpha) : 
   c= np.cos(alpha)*a + np.sin(alpha)*b 
   d=-np.sin(alpha)*a + np.cos(alpha)*b
   return [c,d] 
#---------------------------------------------------------------------------#
def rotate_tensor(a,alpha) :
   #check shape 
   if ( a.shape != (3,3) ) : 
      print ('my_pylib: ERROR Shape {} not supported'.format(a.shape)) 
      quit()
   c=np.cos(alpha) 
   s=np.sin(alpha) 
   [m,c2,s2] = [c*s, c*c, s*s] 
   b=np.zeros([3,3])  

   b[0,0] = c2*a[0,0] + s2*a[2,2] + 2*m*a[0,2]  
   b[0,1] = c* a[0,1] + s*a[1,2] 
   b[0,2] = m*(a[2,2]-a[0,0]) +a[0,2]*(c2-s2) 
   b[1,0] = b[0,1] 
   b[1,1] = a[1,1] 
   b[1,2] =-s*a[0,1]  + c*a[1,2] 
   b[2,0] = b[0,2] 
   b[2,1] = b[1,2]
   b[2,2] = s2*a[0,0] + c2*a[2,2] - 2*m*a[0,2] 
   return b
#---------------------------------------------------------------------------#   
def geti(arr,val): 
    for i in range(len(arr)): 
        if ( arr[i] >= val ): 
            break
    return i 
#---------------------------------------------------------------------------#
def smooth(x,y,n,s):  
    xo =np.zeros(n-(2*s+1))
    smt=np.zeros(n-(2*s+1))
    for i in range(s,n-(s+1)): 
        for ss in range(-s,s+1): 
            smt[i-s] = smt[i-s] + y[i+ss]
            xo[i-s]  = xo[i-s]  + x[i+ss]
    smt=smt/(2*s+1)
    xo =xo/ (2*s+1)
    return [xo,smt]
#---------------------------------------------------------------------------#
def smooth_new(a,ns):
    n = len(a) 
    a_s = np.zeros(n) 
    if ns > n-1: 
        ns = n-1
    # print 2*ns+1, range(-ns,ns+1) 
    for i in range(ns):  
        for j in range(i+ns+1): 
            a_s[i] += a_s[j]/float(i+ns+1)
    
    for i in range(ns,n-ns):  
        for j in range(-ns,ns+1):
            a_s[i] += a[i+j]/float(2*ns+1)

    for i in range(n-ns,n): 
        for j in range(i-ns,n): 
            a_s[i] += a[j]/float(n-i+ns) 
    return a_s 
#---------------------------------------------------------------------------#
def integrate_function(x,y,n):
    sum=0
    for i in range(1,n): 
       sum= sum+ (x[i]-x[i-1])*(y[i] + y[i-1])/2 
    return sum
#---------------------------------------------------------------------------#
def derivative1( x,y,n ):
   der=np.zeros(n) 
   der[0] = (y[1]-y[0])/(x[1]-x[0]) 
   for idum in range(n-2): 
      i=idum+1
      der[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1]) 
   der[n-1]=(y[n-1]-y[n-2])/(x[n-1]-x[n-2])     
   return der
#---------------------------------------------------------------------------#
def derivative2 ( x,y,n ):
   der=np.zeros(n) 
   der[0] = 4*(y[2]-2*y[1]-y[0])/(x[2]-x[0])**2 
   for idum in range(n-2): 
      i=idum+1
      der[i] = 4*(y[i+1]+y[i-1]-2*y[i])/(x[i+1]-x[i-1])**2 
   der[n-1]=4*(y[n-1]-2*y[n-2]+y[n-3])/(x[n-1]-x[n-3])**2     
   return der