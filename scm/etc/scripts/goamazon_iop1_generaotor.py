#!/usr/bin/env python

from   netCDF4 import Dataset
import numpy                 as np
import forcing_file_common   as ffc
import scipy.interpolate
import scm_plotting_routines as spr
import datetime

import matplotlib.pyplot     as plt
import matplotlib

#jhona
from   goamazon.Parameters_forcing_iop import iop1,iop2,trop_atm,mls_atm,hgt_std

import goamazon.calculate_variables_forcing as cal 

import goamazon.read_atm                    as std 

import goamazon.data_own  as down

#to do not use datetimeleap
import  cftime  as cf
import xarray as xr

from    netCDF4         import date2index


matplotlib.use("TkAgg")


#reload(ffc)

#read in raw ARM SGP 1997 summer data input file


case_period_labels = ['A']

case_period_start_datetimes = [cf.datetime(2014,2,15 ,0,0,0),
                              ]
case_period_end_datetimes =   [cf.datetime(2014,3,26 ,0,0,0),
                              ]

lat   = -3.15000 #degrees north
lon   = -60.0    #degrees E
z_sfc = 63 # [m]

#lat   = iop1.lat[0]#degrees north
#lon   = iop1.lon[0]


#ncdump to look at raw input file
#nc_attrs, nc_dims, nc_vars = ffc.ncdump(nc_fid, False)

time=iop1.tsec[:]

T_abs,thetal,qt,ql,qi,tke,\
h_advec_T ,v_advec_T ,dT_dt,h_advec_thil, \
h_advec_qt,v_advec_qt,dqt_dt,\
h_advec_s ,v_advec_s ,ds_dt,rad_heating,\
u_wind,v_wind,u_g,v_g,v_advec_thil,w_sub,omega,\
phi_s,T_surf,p_surf,height,\
levels,sh_flux_sfc,lh_flux_sfc=cal.variables_metpy_scm_ccpp(iop1)


#Call atm trop to get the O3
o3=std.call_read_atm(trop_atm,levels,T_abs)
#o3=std.call_read_atm(mls_atm,levels.magnitude)

data_datetimes = iop1.data


#for t in timestamps:
#    data_datetimes.append(datetime.datetime.utcfromtimestamp(t))
#
case_period_start_indices = []
case_period_end_indices   = []


for i in range(len(case_period_labels)):

    idi,idf=down.date_n(case_period_start_datetimes[i],case_period_end_datetimes[i],data_datetimes[:])

    case_period_start_indices.append(idi)
    case_period_end_indices.append(idf)

    print('Time indices for case {} = [{},{}]. This corresponds to {} - {}'.format(case_period_labels[i], case_period_start_indices[i], case_period_end_indices[i],data_datetimes[case_period_start_indices[i]],data_datetimes[case_period_end_indices[i]]))

    start_t_index = case_period_start_indices[i]
    end_t_index = case_period_end_indices[i]

    #subtract the initial time_offset from all values to get elapsed time since the start of the simulation
    #time = time_offset[start_t_index:end_t_index] - time_offset[start_t_index]
    
    #open processed input file for writing

    writefile_fid = Dataset('../../data/processed_case_input/goamazon_iop1_{}.nc'.format(case_period_labels[i]), 'w', format='NETCDF4')
    writefile_fid.description = "CCPP SCM forcing file for the GOAMAZON 2014 IOP1 case (Period {})".format(case_period_labels[i])

    #create groups for scalars, intitialization, and forcing

    writefile_scalar_grp = writefile_fid.createGroup("scalars")
    writefile_initial_grp = writefile_fid.createGroup("initial")
    writefile_forcing_grp = writefile_fid.createGroup("forcing")
    
    #create dimensions and write them out

    writefile_time_dim = writefile_fid.createDimension('time', None)
    writefile_time_var = writefile_fid.createVariable('time', 'f4', ('time',))
    writefile_time_var[:] = time
    writefile_time_var.units = 's'
    writefile_time_var.description = 'elapsed time since the beginning of the simulation'

    writefile_levels_dim = writefile_fid.createDimension('levels', None)
    writefile_levels_var = writefile_fid.createVariable('levels', 'f4', ('levels',))
    writefile_levels_var[:] = levels
    writefile_levels_var.units = 'Pa'
    writefile_levels_var.description = 'pressure levels'

    #create variables and write them out

    #scalar group

    writefile_lat_var = writefile_scalar_grp.createVariable('lat', 'f4')
    writefile_lat_var[:] = lat
    writefile_lat_var.units = 'degrees N'
    writefile_lat_var.description = 'latitude of column'

    writefile_lon_var = writefile_scalar_grp.createVariable('lon', 'f4')
    writefile_lon_var[:] = lon
    writefile_lon_var.units = 'degrees E'
    writefile_lon_var.description = 'longitude of column'
    
    #initial group


    writefile_height_var = writefile_initial_grp.createVariable('height', 'f4', ('levels',))
    writefile_height_var[:] = height
    writefile_height_var.units = 'm'
    writefile_height_var.description = 'physical height at pressure levels'


    writefile_thetail_var = writefile_initial_grp.createVariable('thetail', 'f4', ('levels',))
    writefile_thetail_var[:] = thetal[:,start_t_index]
    writefile_thetail_var.units = 'K'
    writefile_thetail_var.description = 'initial profile of ice-liquid water potential temperature'





    writefile_qt_var = writefile_initial_grp.createVariable('qt', 'f4', ('levels',))
    writefile_qt_var[:] = qt[:,start_t_index]
    writefile_qt_var.units = 'kg kg^-1'
    writefile_qt_var.description = 'initial profile of total water specific humidity'

    writefile_ql_var = writefile_initial_grp.createVariable('ql', 'f4', ('levels',))
    writefile_ql_var[:] = ql[:,start_t_index]
    writefile_ql_var.units = 'kg kg^-1'
    writefile_ql_var.description = 'initial profile of liquid water specific humidity'

    writefile_qi_var = writefile_initial_grp.createVariable('qi', 'f4', ('levels',))
    writefile_qi_var[:] = qi[:,start_t_index]
    writefile_qi_var.units = 'kg kg^-1'
    writefile_qi_var.description = 'initial profile of ice water specific humidity'

    writefile_u_var = writefile_initial_grp.createVariable('u', 'f4', ('levels',))
    writefile_u_var[:] = u_wind[:,start_t_index]
    writefile_u_var.units = 'm s^-1'
    writefile_u_var.description = 'initial profile of E-W horizontal wind'








    writefile_v_var = writefile_initial_grp.createVariable('v', 'f4', ('levels',))
    writefile_v_var[:] = v_wind[:,start_t_index]
    writefile_v_var.units = 'm s^-1'
    writefile_v_var.description = 'initial profile of N-S horizontal wind'

    writefile_tke_var = writefile_initial_grp.createVariable('tke', 'f4', ('levels',))
    writefile_tke_var[:] = tke[:,start_t_index]
    writefile_tke_var.units = 'm^2 s^-2'
    writefile_tke_var.description = 'initial profile of turbulence kinetic energy'

    writefile_ozone_var = writefile_initial_grp.createVariable('ozone', 'f4', ('levels',))
    writefile_ozone_var[:] = o3
    writefile_ozone_var.units = 'kg kg^-1'
    writefile_ozone_var.description = 'initial profile of ozone mass mixing ratio'
    
    #forcing group

    writefile_p_surf_var = writefile_forcing_grp.createVariable('p_surf', 'f4', ('time',))
    writefile_p_surf_var[:] = p_surf[start_t_index:end_t_index]
    writefile_p_surf_var.units = 'Pa'
    writefile_p_surf_var.description = 'surface pressure'

    writefile_T_surf_var = writefile_forcing_grp.createVariable('T_surf', 'f4', ('time',))
    writefile_T_surf_var[:] = T_surf[start_t_index:end_t_index]
    writefile_T_surf_var.units = 'K'
    writefile_T_surf_var.description = 'surface absolute temperature'

    writefile_sh_flux_sfc_var = writefile_forcing_grp.createVariable('sh_flux_sfc', 'f4', ('time',))
    writefile_sh_flux_sfc_var[:] = sh_flux_sfc[start_t_index:end_t_index]
    writefile_sh_flux_sfc_var.units = 'K m s^-1'
    writefile_sh_flux_sfc_var.description = 'surface sensible heat flux'

    writefile_lh_flux_sfc_var = writefile_forcing_grp.createVariable('lh_flux_sfc', 'f4', ('time',))
    writefile_lh_flux_sfc_var[:] = lh_flux_sfc[start_t_index:end_t_index]
    writefile_lh_flux_sfc_var.units = 'kg kg^-1 m s^-1'
    writefile_lh_flux_sfc_var.description = 'surface latent heat flux'

    #fig   = plt.figure()
    #ax    = plt.axes()
    #plt.plot(u_wind,height)

    #fig   = plt.figure()
    #ax    = plt.axes()
    #plt.plot(v_wind,height)

    #fig   = plt.figure()
    #ax    = plt.axes()
    #plt.plot(w_sub,height)
    #plt.show()

    fig   = plt.figure()
    ax    = plt.axes()
    plt.plot(time,lh_flux_sfc)

    fig   = plt.figure()
    ax    = plt.axes()
    plt.plot(time,sh_flux_sfc)
    plt.show()
    plt.show()

    exit()









    writefile_w_ls_var = writefile_forcing_grp.createVariable('w_ls', 'f4', ('levels','time',))
    writefile_w_ls_var[:] = w_sub[:,start_t_index:end_t_index]
    writefile_w_ls_var.units = 'm s^-1'
    writefile_w_ls_var.description = 'large scale vertical velocity'

    writefile_omega_var = writefile_forcing_grp.createVariable('omega', 'f4', ('levels','time',))
    writefile_omega_var[:] = omega[:,start_t_index:end_t_index]
    writefile_omega_var.units = 'Pa s^-1'
    writefile_omega_var.description = 'large scale pressure vertical velocity'


    writefile_u_g_var = writefile_forcing_grp.createVariable('u_g', 'f4', ('levels','time',))
    writefile_u_g_var[:] = u_g[:,start_t_index:end_t_index]
    writefile_u_g_var.units = 'm s^-1'
    writefile_u_g_var.description = 'large scale geostrophic E-W wind'

    writefile_v_g_var = writefile_forcing_grp.createVariable('v_g', 'f4', ('levels','time',))
    writefile_v_g_var[:] = v_g[:,start_t_index:end_t_index]
    writefile_v_g_var.units = 'm s^-1'
    writefile_v_g_var.description = 'large scale geostrophic N-S wind'

    writefile_u_nudge_var = writefile_forcing_grp.createVariable('u_nudge', 'f4', ('levels','time',))
    writefile_u_nudge_var[:] = u_wind[:,start_t_index:end_t_index]
    writefile_u_nudge_var.units = 'm s^-1'
    writefile_u_nudge_var.description = 'E-W wind to nudge toward'

    writefile_v_nudge_var = writefile_forcing_grp.createVariable('v_nudge', 'f4', ('levels','time',))
    writefile_v_nudge_var[:] = v_wind[:,start_t_index:end_t_index]
    writefile_v_nudge_var.units = 'm s^-1'
    writefile_v_nudge_var.description = 'N-S wind to nudge toward'

    writefile_T_nudge_var = writefile_forcing_grp.createVariable('T_nudge', 'f4', ('levels','time',))
    writefile_T_nudge_var[:] = T_abs[:,start_t_index:end_t_index]
    writefile_T_nudge_var.units = 'K'
    writefile_T_nudge_var.description = 'absolute temperature to nudge toward'

    writefile_thil_nudge_var = writefile_forcing_grp.createVariable('thil_nudge', 'f4', ('levels','time',))
    writefile_thil_nudge_var[:] = thetal[:,start_t_index:end_t_index]
    writefile_thil_nudge_var.units = 'K'
    writefile_thil_nudge_var.description = 'potential temperature to nudge toward'

    writefile_qt_nudge_var = writefile_forcing_grp.createVariable('qt_nudge', 'f4', ('levels','time',))
    writefile_qt_nudge_var[:] = qt[:,start_t_index:end_t_index]
    writefile_qt_nudge_var.units = 'kg kg^-1'
    writefile_qt_nudge_var.description = 'q_t to nudge toward'

    writefile_rad_heating_var = writefile_forcing_grp.createVariable('dT_dt_rad', 'f4', ('levels','time',))
    writefile_rad_heating_var[:] = rad_heating[:,start_t_index:end_t_index]
    writefile_rad_heating_var.units = 'K s^-1'
    writefile_rad_heating_var.description = 'prescribed radiative heating rate'

    writefile_h_advec_thil_var = writefile_forcing_grp.createVariable('h_advec_thetail', 'f4', ('levels','time',))
    writefile_h_advec_thil_var[:] = h_advec_thil[:,start_t_index:end_t_index]
    writefile_h_advec_thil_var.units = 'K s^-1'
    writefile_h_advec_thil_var.description = 'prescribed theta_il tendency due to horizontal advection'

    writefile_v_advec_thil_var = writefile_forcing_grp.createVariable('v_advec_thetail', 'f4', ('levels','time',))
    writefile_v_advec_thil_var[:] = v_advec_thil[:,start_t_index:end_t_index]
    writefile_v_advec_thil_var.units = 'K s^-1'
    writefile_v_advec_thil_var.description = 'prescribed theta_il tendency due to vertical advection'

    writefile_h_advec_qt_var = writefile_forcing_grp.createVariable('h_advec_qt', 'f4', ('levels','time',))
    writefile_h_advec_qt_var[:] = h_advec_qt[:,start_t_index:end_t_index]
    writefile_h_advec_qt_var.units = 'kg kg^-1 s^-1'
    writefile_h_advec_qt_var.description = 'prescribed q_t tendency due to horizontal advection'

    writefile_v_advec_qt_var = writefile_forcing_grp.createVariable('v_advec_qt', 'f4', ('levels','time',))
    writefile_v_advec_qt_var[:] = v_advec_qt[:,start_t_index:end_t_index]
    writefile_v_advec_qt_var.units = 'kg kg^-1 s^-1'
    writefile_v_advec_qt_var.description = 'prescribed q_t tendency due to vertical advection'

    #close processed input file
    writefile_fid.close()
    

#close raw input file
#writefile_fid.close()
