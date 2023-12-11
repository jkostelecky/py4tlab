import os
import my_pylib as mp
import matplotlib.pyplot as plt
import copy
import seaborn as sns
import sys
import glob

#---------------------------------------------------------------------------#
# files
#---------------------------------------------------------------------------#
current_path = os.getcwd() + '/'
paths        = [
                current_path + '../example_data/'
                ] 
# name         = [
                # ]
path_grid    = paths[0]
# path2param   = paths[0] + 'param_data/'
path2fig     = current_path +  'figs/'
os.makedirs(path2fig, exist_ok=True)

# plot settings 
plt.style.use('tex')
plt.rcParams['figure.dpi'] = 300
shading = 'nearest' # 'gouraud'
width  = 427
alpha  = 0.7 # opacity
pad    = 0.3
colors = ['red','orange','green','blue']
cmap_rb   = copy.copy(plt.cm.get_cmap("RdBu_r")) 
plt.close('all')
cmap_sns = sns.color_palette("viridis_r", as_cmap=True)

#---------------------------------------------------------------------------#
# grid
grid = mp.DnsGrid(paths[0], 'grid')

# exclude ini-step (sometimes weird values)
exclude = True
mean    = True

name_out = 'dns.out'
name_obs = 'dns.obs'

#%% 
'''
# concatenate all files
def cat(path, otype, delete_cat=True):
    print('--------------------------------------------------')
    os.chdir(path)
    name = 'dns'
    # delete alrady cat file
    cat_file = name + '.' + otype + '_cat'
    del_file = path + cat_file
    print(del_file)
    if os.path.exists(del_file):
        if not delete_cat:
            print('ERROR File exists already, no merging possible')
            sys.exit()
        else:
            print('Delete file           :  ', del_file)
            os.remove(del_file)    
    # sort
    olist = sorted([n[len(path):] for n in glob.glob(path + name + '.' + str(otype) + '*')])
    if 'dns.out' in olist: olist += [olist.pop(0)] # move to the end of lis
    # concatenate
    command = 'cat '
    for i in olist: command += ' ' + str(i)
    command += ' > ' + name + '.' + otype + '_cat' 
    os.system(command)
    # add cat file to the end of olist
    olist.append(cat_file)
    return olist

out_list = cat(paths[0], otype='out')
obs_list = cat(paths[0], otype='obs')
'''
out_list = [name_out]
#%% 
if exclude:
    e = 1
else:
    e = 0

# plot settings
plt.close('all')
width = 1000
ratio = 0.8
ma = '.'
ms = 1
lw = 0.3
   
# plotting
for lout in out_list:
    # read dns_out
    out = mp.DnsOut(paths[0], lout)
    restart = out.read_dnso()

    
    fig, axs = plt.subplots(3, 2, sharex=True, sharey=False, figsize=mp.set_size(width=width, ratio=ratio), layout="constrained", gridspec_kw={"hspace": 0.1})#, layout="constrained")
    plt.suptitle("Quicklook for '{}'".format(lout))
    
    ii = 0; jj = 0; key='dt'; xvar='it'
    axs[ii,jj].set_title(r'$\Delta t_{sim}\cdot (1/f)$')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),7))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 0; jj += 1; key='dt_it'
    if lout == 'dns.out_cat': axs[ii,jj].set_ylim(25,35)
    axs[ii,jj].set_title(r'$\frac{\Delta t_{real}}{it}$ (time for one iteration in [s])')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    k=0
    for i in restart:
        if k==0: axs[ii,jj].axvline(x=i, color='red',label='restart'); k+=1
        else:    axs[ii,jj].axvline(x=i, color='red')
    k=0
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),2))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 1; jj = 0 ; key='cfl'
    axs[ii,jj].set_title(r'$CFL$-number')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),2))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 1; jj += 1 ; key='dif'
    axs[ii,jj].set_title(r'$Dif$-number')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),3))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 2; jj = 0 ; key='dil_min'
    axs[ii,jj].set_title(r'$Dil_{min}$')
    axs[ii,jj].set_xlabel(r'iteration step')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid(True)
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean()))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 2; jj += 1 ; key='dil_max'
    axs[ii,jj].set_title(r'$Dil_{max}$')
    axs[ii,jj].set_xlabel(r'iteration step')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean()))), zorder=-1)
        axs[ii,jj].legend()
    
    output_name = 'quicklook_out' + lout[7:] +'.pdf'
    plt.savefig(path2fig+output_name, dpi=800, format='pdf')
    plt.show()
    # plt.close()
'''    
#%%
for lobs in obs_list:
    # read dns.obs
    out = mp.DnsOut(paths[0], lobs)
    restart = out.read_dnso()
    
    plt.close('all')
    
    width = 1000
    ratio = 1
    
    fig, axs = plt.subplots(4, 2, sharex=True, sharey=False, figsize=mp.set_size(width=width, ratio=ratio), layout="constrained", gridspec_kw={"hspace": 0.1})#, layout="constrained")
    plt.suptitle("Quicklook for '{}'".format(lobs))
    
    ii = 0; jj = 0; key='ub'; xvar='time'
    axs[ii,jj].set_title(r'$u_{bulk}$')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),4))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 0; jj += 1; key='wb'
    axs[ii,jj].set_title(r'$w_{bulk}$')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),4))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 1; jj = 0 ; key='al1'
    axs[ii,jj].set_title(r'$\alpha(n_i=1)$ [deg]')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),2))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 1; jj += 1 ; key='alny'
    axs[ii,jj].set_title(r'$\alpha(n_i=ny)$ [deg]')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),2))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 2; jj = 0 ; key='uy1'
    axs[ii,jj].set_title(r'$\partial_y u$')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid(True)
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),0))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 2; jj += 1 ; key='wy1'
    axs[ii,jj].set_title(r'$\partial_y w$')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),0))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 3; jj = 0 ; key='ent'
    axs[ii,jj].set_xlabel(r'simulation time')
    axs[ii,jj].set_title(r'$\int\Omega^2dy$ (integrated entstophy)')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid()
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),0))), zorder=-1)
        axs[ii,jj].legend()
    
    ii = 3; jj += 1 ; key='s1y'
    axs[ii,jj].set_xlabel(r'simulation time')
    axs[ii,jj].set_title(r'$\partial_y s_1$ (scalar gradient)')
    axs[ii,jj].margins(0.1,0.1)
    axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
    axs[ii,jj].plot(out.data[xvar][e:], out.data[key][e:], color='blue', marker=ma, linewidth=lw, markersize=ms)
    axs[ii,jj].grid(True)
    if mean:
        axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean {}'.format(str(round(out.data[key][e:].mean(),0))), zorder=-1)
        axs[ii,jj].legend()
    
    output_name = 'quicklook_obs' + lobs[7:] + '.pdf'
    plt.savefig(path2fig+output_name, dpi=800, format='pdf')
    # plt.show()
    plt.close()
'''