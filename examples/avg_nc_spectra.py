#!/usr/bin/env python3
import os
import my_pylib as mp
#---------------------------------------------------------------------------#
# path to spectra fields
current_path = os.getcwd() + '/'
path         = current_path + '../example_data/spec/'
index        = [0,10,10]
filetag      = ["xsp", "zsp"]#, "pow", "pha"]  # ['xsp', 'zsp', 'pow', 'pha', ... ]
variables    = ["uu", "vv", "ww", "pp"]#, "11", "22"] #, 'uv','uw','vw']

for tag in filetag:
    sp_smooth = mp.Spectra(path, index, tag, variables)
    sp_smooth.read_spectra()
    sp_smooth.spectra_save_nc(sp_smooth.sp, delete=True)
    del sp_smooth