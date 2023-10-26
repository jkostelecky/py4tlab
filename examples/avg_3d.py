import os
import my_pylib as mp
#-----------------------------------------------------------------------------#
current_path = os.getcwd() + '/'
path         = current_path + '../example_data/'
# nphase       = 4
#---------------------------------------------------------------------------#
# average in x, center line, phase average and save
fields = mp.AVG_Fields(path, ftype='flow', scalar=True, pressure=False, num_scal=2, skip=0)
fields.avg_x()
fields.avg_save_nc('avg_x')
fields.avg_xcl()
fields.avg_save_nc('avg_cl')
# fields.avg_phase(nphase)
# fields.avg_save_nc('avg_pa')
fields.avg_save_nc('avg_cl')

# average in z and save
fields = mp.AVG_Fields(path, ftype='flow', scalar=True, pressure=False, num_scal=2, skip=0)
fields.avg_z()
fields.avg_save_nc('avg_z')