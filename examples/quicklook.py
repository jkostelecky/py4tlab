import os
import my_pylib as mp
import matplotlib.pyplot as plt
from   matplotlib import rc
import copy
import seaborn as sns
#---------------------------------------------------------------------------#
# files
#---------------------------------------------------------------------------#
current_path = os.getcwd() + '/'
paths        = [
                current_path
                ] 
# name         = [
                # ]
path_grid    = paths[0]
# path2param   = paths[0] + 'param_data/'
path2fig     = current_path +  'figs_quicklook/'
os.makedirs(path2fig, exist_ok=True)

current_path = os.getcwd() + '/'
paths        = [current_path + '../example_data/quicklook/']
path2fig     = current_path +  'figs/'
os.makedirs(path2fig, exist_ok=True)

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
# exclude ini-step (sometimes weird values)
exclude = True
mean    = True

#%% PLOTS
if exclude:
    e = 1
else:
    e = 0
    
# read dns_out
out = mp.DnsOut(paths[0], 'dns.out')

plt.close('all')
width = 1000
ratio = 0.8
fig, axs = plt.subplots(3, 2, sharex=True, sharey=False, figsize=mp.set_size(width=width, ratio=ratio), layout="constrained", gridspec_kw={"hspace": 0.1})#, layout="constrained")
plt.suptitle('Quicklook for dns.out')

ii = 0; jj = 0; key='dt'
axs[ii,jj].set_title(r'$\Delta t_{sim}$')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 0; jj += 1; key='dt_it'
axs[ii,jj].set_title(r'$\frac{\Delta t_{real}}{it}$ (time for one iteration in [s])')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 1; jj = 0 ; key='cfl'
axs[ii,jj].set_title(r'$CFL$-number')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 1; jj += 1 ; key='dif'
axs[ii,jj].set_title(r'$Dif$-number')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 2; jj = 0 ; key='dil_min'
axs[ii,jj].set_title(r'$Dil_{min}$')
axs[ii,jj].set_xlabel(r'iteration step')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid(True)
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 2; jj += 1 ; key='dil_max'
axs[ii,jj].set_title(r'$Dil_{max}$')
axs[ii,jj].set_xlabel(r'iteration step')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

output_name = 'quicklook_out.pdf'
plt.savefig(path2fig+output_name, dpi=800, format='pdf')
plt.close()


#%%
# read dns.obs
out = mp.DnsOut(paths[0], 'dns.obs')

plt.close('all')
width = 1000
ratio = 1
fig, axs = plt.subplots(4, 2, sharex=True, sharey=False, figsize=mp.set_size(width=width, ratio=ratio), layout="constrained", gridspec_kw={"hspace": 0.1})#, layout="constrained")
plt.suptitle('Quicklook for dns.obs')

ii = 0; jj = 0; key='ub'
axs[ii,jj].set_title(r'$u_{bulk}$')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 0; jj += 1; key='wb'
axs[ii,jj].set_title(r'$w_{bulk}$')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 1; jj = 0 ; key='al1'
axs[ii,jj].set_title(r'$\alpha(n_i=1)$ [deg]')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 1; jj += 1 ; key='alny'
axs[ii,jj].set_title(r'$\alpha(n_i=ny)$ [deg]')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 2; jj = 0 ; key='uy1'
axs[ii,jj].set_title(r'$\partial_y u$')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid(True)
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 2; jj += 1 ; key='wy1'
axs[ii,jj].set_title(r'$\partial_y w$')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 3; jj = 0 ; key='ent'
axs[ii,jj].set_xlabel(r'iteration step')
axs[ii,jj].set_title(r'$\int\Omega^2dy$ (integrated entstophy)')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid()
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

ii = 3; jj += 1 ; key='s1y'
axs[ii,jj].set_xlabel(r'iteration step')
axs[ii,jj].set_title(r'$\partial_y s_1$ (scalar gradient)')
axs[ii,jj].margins(0.1,0.1)
axs[ii,jj].ticklabel_format(style='sci', useOffset=False)
axs[ii,jj].plot(out.data['it'][e:], out.data[key][e:], color='blue', marker='x')
axs[ii,jj].grid(True)
if mean:
    axs[ii,jj].axhline(y=out.data[key][e:].mean(), color='darkred', label='mean')
    axs[ii,jj].legend()

output_name = 'quicklook_obs.pdf'
plt.savefig(path2fig+output_name, dpi=800, format='pdf')
plt.close()