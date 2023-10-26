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

# read dns_out
out = mp.DnsOut(path_fields, 'dns.out', 1)
var = ['iteration_step', 'time_total','time_step','cfl_number','dif_number','visc','dil_min','dil_max']
#---------------------------------------------------------------------------#
# plots
plt.figure(figsize=size)
plt.title(r'dilatation')
plt.xlabel(r"iterations")
plt.ylabel(r'$dil_{max}$')
# plt.ylim(0.64,0.69)
plt.scatter(out.data['it'], out.data['dil_max'], color='blue', marker='x', label=r'$dil_{max}$')
plt.legend()
plt.grid(True)
plt.show()