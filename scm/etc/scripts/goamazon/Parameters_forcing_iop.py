#To load the variable  in the datafile
import  goamazon.variabletoload_go_iop_scm   as scm

computer='/home/jhona'

file_out    ="%s/repositories/shallow_2/fig"%(computer)

file_scm_iop1  = '../../data/raw_case_input/GOAMAZON/forcing_goamazon_IOP1.nc'
calendar = ['seconds since 2014-2-15T00:00:00 +04:00:00','gregorian']
iop1 = scm.ncload(file_scm_iop1,calendar)

file_scm_iop2  = '../../data/raw_case_input/GOAMAZON/forcing_goamazon_IOP2.nc'
calendar = ['seconds since 2014-9-1T00:00:00 +04:00:00','gregorian']
iop2 = scm.ncload(file_scm_iop2,calendar)


trop_atm  = '../../data/raw_case_input/GOAMAZON/tro.atm'
mls_atm  = '../../data/raw_case_input/GOAMAZON/mls.atm'

hgt_std  = '../../data/raw_case_input/GOAMAZON/hgt_std.atm'


