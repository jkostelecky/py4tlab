#!/usr/bin/env python3
import matplotlib.pyplot as plt
from   matplotlib import rc
import os
import my_pylib as mp
import copy
import numpy as np
#-----------------------------------------------------------------------------#
# file
index        = 10
forcing_flow = 'cpg' # constant pressure gradient forcing
current_path = os.getcwd() + '/'
path_fields  = current_path + '../example_data/'
path_grid    = path_fields
path2fig     = current_path +  'figs/'
os.makedirs(path2fig, exist_ok=True)
name         = 'channel'

# plot settings 
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
font = {'family':'monospace', 'weight':'bold', 'size':14}; rc('font', **font) # font
plt.rcParams['figure.dpi'] = 210
size    = (8,6)
shading = 'nearest'#'gouraud'
cmap    = copy.copy(plt.cm.get_cmap("RdBu_r")) 
plt.close('all')
#---------------------------------------------------------------------------#
# read grid
grid = mp.DnsGrid(path_grid, 'grid')

# read flow fields
fld = mp.Field(path_fields,var='flow',index=index, forcing=forcing_flow)
fld.read_3d_field()
ub = mp.ubulk(fld.data['u'], grid.y) 

# # read eps fields
# eps = mp.Field(path_fields,var='eps_int',index=0)
# eps.read_3d_field()

# read scalar fields
scal = mp.Field(path_fields,var='phi1',index=index, forcing=forcing_flow)
scal.read_3d_field()
fld.data.update(scal.data)
del scal.data

# # read pressure fields
# pre = mp.Field(path_fields,var='pre',index=index, forcing=forcing_flow)
# pre.read_3d_field()
# fld.data.update(pre.data)
# del pre.data
#---------------------------------------------------------------------------#
# plot
list_var = [key for key in fld.data]
# list_var.remove('eps')

xpos = grid.nx//2
zpos = grid.nz//2

for key in list_var:
    ysize = grid.ly/grid.lx
    size = (8,8*ysize+1.0)
    fig = plt.figure(figsize=size)
    ax  = plt.axes()
    ax.set_title(r'instantaneous '+ key)
    ax.set_xlabel(r'$x/\delta$')
    ax.set_ylabel(r'$y/\delta$')
    ax.set_aspect('equal')
    data = fld.data[key][:,:,zpos]
    cmp_max = max(abs(np.nanmin(data)), np.nanmax(data))
    cmp_min = - cmp_max
    if key in ['phi1', 'phi2']:
        idx = np.where(data == 1.); data[idx] = np.nan
    else:
        idx = np.where(data == 0.); data[idx] = np.nan
    cmap.set_bad(color='black')
    if key in ['u', 'phi1', 'phi2']: cmp_min = 0
    pcm = ax.pcolormesh(grid.x,grid.y,data.T, shading=shading, cmap=cmap, vmin=cmp_min, vmax=cmp_max)
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    fig.colorbar(pcm, cax=cax)
    # ax.add_patch(Rectangle((z[36],0), 0.1, y[19], fill=False, edgecolor='grey', lw=4))
    output_name = name + '_' + key + '_xy.jpg'
    plt.savefig(path2fig+output_name, dpi=400, format='jpg')
    plt.show()
    # ---------------------------------------------------------------------- # 
    ysize = grid.ly/grid.lz
    size = (8,8*ysize+1.0)
    fig = plt.figure(figsize=size)
    ax  = plt.axes()
    ax.set_title(r'instantaneous '+ key)
    ax.set_xlabel(r'$z/\delta$')
    ax.set_ylabel(r'$y/\delta$')
    ax.set_aspect('equal')
    data = fld.data[key][xpos,:,:]
    cmp_max = max(abs(np.nanmin(data)), np.nanmax(data))
    cmp_min = - cmp_max
    if key in ['phi1', 'phi2']:
        idx = np.where(data == 1.); data[idx] = np.nan
    else:
        idx = np.where(data == 0.); data[idx] = np.nan    
    cmap.set_bad(color='black')
    if key in ['u', 'phi1', 'phi2']: cmp_min = 0
    pcm = ax.pcolormesh(grid.z, grid.y, data, shading=shading, cmap=cmap, vmin=cmp_min, vmax=cmp_max)
    cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
    fig.colorbar(pcm, cax=cax)
    # ax.add_patch(Rectangle((z[36],0), 0.1, y[19], fill=False, edgecolor='grey', lw=4))
    output_name = name + '_' + key + '_yz.jpg'
    plt.savefig(path2fig+output_name, dpi=400, format='jpg')
    plt.show()
    #
    plt.close('all')