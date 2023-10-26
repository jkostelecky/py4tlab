#!/usr/bin/env python3
import my_pylib as mp
import os
import matplotlib.pyplot as plt
import netCDF4 as nc
import warnings; warnings.filterwarnings("ignore", category=DeprecationWarning) 
#-----------------------------------------------------------------------------#
# path to flow fields
current_path = os.getcwd() + '/'
path         = current_path + '../example_data/planes/'

# read and plot planes
index      = 10
num_pl     = 1      # number of slices
num_scal   = 2      # number of scalars
pl_log     = True   # log of enstrophy (if scalar, log of scalar gradient)
ind        = [index,index,1]

# plane_type = 'xy' # planesK
plane_type = 'xz'; num_pl = 3 # planesJ
# plane_type = 'yz' # planesI

#-----------------------------------------------------------------------------#
# processing data
pl = mp.Planes(path=path, index=ind, planetype=plane_type, num_pl=num_pl, num_scal=num_scal, pl_log=pl_log)
pl.read_planes()
pl.planes_save_nc()

# plot
plt.close('all')
for l in range(num_pl):
    for key in pl.list_var:
        plt.figure(figsize=(8,4))
        plt.title(key + str(l))
        if   plane_type == 'xy':
            plt.pcolormesh(pl.x1, pl.x2, pl.data[key][l,:,:].T, cmap='RdBu_r')
        elif plane_type == 'xz':
            plt.pcolormesh(pl.x1, pl.x2, pl.data[key][l,:,:].T, cmap='RdBu_r')
        elif plane_type == 'yz':
            plt.pcolormesh(pl.x2, pl.x1, pl.data[key][l,:,:],   cmap='RdBu_r')
        plt.colorbar()
        plt.show()
        
# # check nc file
# f = nc.Dataset(pl.fname_dst, 'r', format='netCDF-4')
# print(f.variables.keys())
# print(f.variables['u'][:,:,:].shape)
# for l in range(num_pl):
#     for key in list_var:
#         plt.figure(figsize=(8,4))
#         plt.title(key + str(l))
#         if   plane_type == 'xy':
#             plt.pcolormesh(x1, x2, f.variables[key][l,:,:].T, cmap='RdBu_r')
#         elif plane_type == 'xz':
#             plt.pcolormesh(x1, x2, f.variables[key][l,:,:].T, cmap='RdBu_r')
#         elif plane_type == 'yz':
#             plt.pcolormesh(x2, x1, f.variables[key][l,:,:],   cmap='RdBu_r')
#         plt.colorbar()
#         plt.show()