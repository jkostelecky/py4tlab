#!/usr/bin/env python3
import os
import my_pylib as mp
#-----------------------------------------------------------------------------#
# file
current_path = os.getcwd() + '/'
path2ini     = current_path + '../example_data/'
path2param   = current_path +  'param_data/'
os.makedirs(path2param, exist_ok=True)
#---------------------------------------------------------------------------#

# create and save dict with params
p = mp.InputParam(path2ini, path2param, 'param', sim_type='not specified')

# print dict
for k in p.param.keys():
    print(k, p.param[k])