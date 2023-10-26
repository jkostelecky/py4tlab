import my_pylib as mp
import os
#-----------------------------------------------------------------------------#
# files
current_path = os.getcwd() + '/'
path         = current_path + '../example_data/'

skip  = 0

name1 = 'avg'
name2 = 'avg1s'
name3 = 'avg2s'

path1 = path + 'avg/'
path2 = path + 'avgs/'
path3 = path + 'avgs/'
#-----------------------------------------------------------------------------#
# avgerage all .nc fiels in directory
avg = mp.Statistics(path1, fname=name1)
avg.merge_nc_all(delete=True)
avg.avg_all(fin='avg',    fout='avg_t',  skip=skip,   ensemble=True,  twa=False,   delete=True)
del avg
#-----------------------------------------------------------------------------#
# avgerage all .nc fiels in directory
avg = mp.Statistics(path2, fname=name2)
avg.merge_nc_all(delete=True)
avg.avg_all(fin='avg1s',    fout='avg1s_t',  skip=skip,   ensemble=True, twa=False,   delete=True)
del avg
#-----------------------------------------------------------------------------#
# avgerage all .nc fiels in directory
avg = mp.Statistics(path3, fname=name3)
avg.merge_nc_all(delete=True)
avg.avg_all(fin='avg2s',    fout='avg2s_t',  skip=skip,   ensemble=True, twa=False,   delete=True)
del avg