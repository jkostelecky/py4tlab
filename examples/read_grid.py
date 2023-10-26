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
# in case of channel flow
re_tau = 180

# read grid and flow fields
grid = mp.DnsGrid(path_grid, 'grid')

# grid spacing
dy = np.diff(grid.y)

# positions of mid nodes (dy is plotted here)
ym = (grid.y[:-1] + grid.y[1:]) / 2 

# print information
print('--------------------------------------------------')
print('stretched grid information in vertical direction')
print('origin     :', grid.y[0])
print('end        :', grid.y[-1])
print('min step   :', dy.min())
print('max step   :', dy.max())
print('--------------------------------------------------')
print('stretching in viscous units for re_tau = ', re_tau)
print('min step + :', dy.min()*re_tau)
print('max step + :', dy.max()*re_tau)
print('--------------------------------------------------')

#---------------------------------------------------------------------------#
# plot settings 
plt.rcParams['figure.dpi'] = 210 
size    = (8,6)
shading = 'nearest'#'gouraud'
figs    = 'figs' 
plt.close('all')

#-----------------------------------------------------------------------------#
# plot vertical grid spacing
plt.figure(figsize=size)
plt.grid(True)
plt.xlim(0,ym.max())
plt.ylim(0,1.2*dy.max())
plt.xlabel('y mid-node postions')
plt.ylabel('delta y')
plt.plot(grid.y[0],dy[0],   'o', color='red',   label='origin')
plt.plot(grid.y[-1],dy[-1], 'o', color='black', label='end')
plt.plot(ym, dy, marker='.',label='dy')
plt.legend(loc=1)
plt.show()

# plot vertical grid spacing in visc. units (channel flow)
plt.figure(figsize=size)
plt.grid(True)
plt.xlim(0,ym.max())
plt.ylim(0,1.2*dy.max()*re_tau)
plt.xlabel('y mid-node postions')
plt.ylabel('delta y+')
plt.plot(grid.y[0],dy[0]*re_tau,   'o', color='red',   label='origin')
plt.plot(grid.y[-1],dy[-1]*re_tau, 'o', color='black', label='end')
plt.plot(ym, dy*re_tau, marker='.',label='dy')
plt.legend(loc=1)
plt.show()