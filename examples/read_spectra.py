#!/usr/bin/env python3
import matplotlib.pyplot as plt
from   matplotlib import rc
import os
import my_pylib as mp
import copy
import numpy as np
#-----------------------------------------------------------------------------#
# file
forcing_flow = 'cpg' # constant pressure gradient forcing
current_path = os.getcwd() + '/'
path_spec    = current_path + '../example_data/spec/'
path_grid    = current_path + '../example_data/'
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
colors  = ['red','blue','orange','black']
prefix_fig = "spectra"
plt.close('all')
#---------------------------------------------------------------------------#
# average spectras to nc fields
index        = [10,10,1] 
# index        = [0,10,10] 
variables    = ["uu", "vv", "ww", "pp"]
filetag      = 'xsp'
nametag      = "E"
re_tau       = 180 
y_pos        = [1,5,10,50,100] # vertical position

# read grid field
grid = mp.DnsGrid(path_grid, 'grid')

# read spectra binaries
sp = mp.Spectra(path_spec, index, filetag, variables)
sp.read_spectra()
#---------------------------------------------------------------------------#
# lines
# xl    = np.arange(2, 100, 1)
# ylvel = 3e-3 * xl**(-5/3)
# ylpre = 4e-4 * xl**(-7/3)
for j in y_pos:
    # plot single spectra
    plt.figure(figsize=size)
    i = 0
    for var in variables:
        plt.plot(sp.sp[i,j,:], color=colors[i], linestyle='solid', alpha=0.7, label=r'{}'.format(var))    
        i += 1
    plt.title(r'Single Spectra {} at vertical postion $y^+={}^+$'.format("in {}-direction".format(filetag[0]), round(grid.y[j]*re_tau,1)))
    plt.xlabel(r'wave number')
    plt.ylabel(r'mode energy')
    plt.ylim(1e-24,1e-0)
    plt.xlim(1, 1000)
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout(pad=0.1)
    plt.legend(loc=1)
    plt.grid('True')
    # output_prefix_fig = prefix_fig +'_x_'+str(j)+'.svg'
    # plt.savefig(path2fig+output_prefix_fig)
    plt.show()
    #---------------------------------------------------------------------------#