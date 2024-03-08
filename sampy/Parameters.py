########################################################
#File to define the direction of the SAM runs.
########################################################

### My own files

##Adress of the computer
from     files_direction  import *

import   sam_python.var_files.var_to_load_ccpp  as ccpp

########################################################
#Files of SAM to load.
########################################################


var1d             = ['shf','lhf', 'tprcp_accum']
var2d             = ['T','ql','w_ls']
vars_diurnal      = ['T','ql','w_ls']


#Make the first dictionary 
name       = 'goamazon_iop1'
file1      = '%s/git_repositories/ccpp-scm/scm/run/output_goamazon_iop1_A_SCM_GFS_v16/output.nc'%(computer)
data         = [(2014,2,15,18),(2014,3,26,0)]
data_diurnal = [(2014,2,19,8),(2014,2,19,18)]
cal          = ["seconds  since 2014-2-15 00:00:00 +00:00:00",'gregorian']
goa1       = ccpp.ncload(name,data,file1,cal,var1d,var2d,vars_diurnal,data_diurnal) 

name       = 'goamazon_iop2'
file2      = '%s/git_repositories/ccpp-scm/scm/run/output_goamazon_iop2_A_SCM_GFS_v16/output.nc'%(computer)
data       = [(2014,9,1,0),(2014,10,10,0)]
cal        = ["seconds  since 2014-9-01 00:00:00 +00:00:00",'gregorian']
goa2       = ccpp.ncload(name,data,file2,cal,var1d,var2d,vars_diurnal,data_diurnal) 

