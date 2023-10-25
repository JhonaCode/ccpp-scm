################################################# 
# Program to read variable of a nc file 
# using python with NetCdf
# Create by: Jhonatan Aguirre 
# Date:06/02/2020
# working:yes
#################################################

# Python library to work with Netcdf4 
from    netCDF4         import Dataset

# To save the time coordinate in specific format 
from    netCDF4         import num2date, date2num

#to do not use datetimeleap
import  cftime as  cft

#cft.num2date, cft.num2pydate


class variables(object):

    def __init__(self):

	    self.bdate              ='bdate'
	    self.lev                ='lev'
	    self.lon                ='lon'
	    self.lat                ='lat'
	    self.tsec               ='tsec'
	    self.phis               ='phis'
	    self.prec_srf           ='prec_srf'
	    self.lhflx              ='lhflx'
	    self.shflx              ='shflx'
	    self.Ps                 ='Ps'
	    self.Tsair              ='Tsair'
	    self.RH_srf             ='RH_srf'
	    self.wspd_srf           ='wspd_srf'        
	    self.usrf               ='usrf'
	    self.vsrf               ='vsrf'
	    self.rad_net_srf        ='rad_net_srf'
	    self.lw_net_toa         ='lw_net_toa'
	    self.sw_net_toa         ='sw_net_toa'
	    self.sw_dn_toa          ='sw_dn_toa'
	    self.cld_low            ='cld_low'
	    self.cld_mid            ='cld_mid'
	    self.cld_high           ='cld_high'
	    self.cld_tot            ='cld_tot'
	    self.cld_thick          ='cld_thick'
	    self.cld_top            ='cld_top'
	    self.LWP                ='LWP'
	    self.dh2odt_col         ='dh2odt_col'
	    self.h2o_adv_col        ='h2o_adv_col'
	    self.evap_srf           ='evap_srf'
	    self.dsdt_col           ='dsdt_col'
	    self.s_adv_col          ='s_adv_col'
	    self.rad_heat_col       ='rad_heat_col'
	    self.qs                 ='qs'
	    self.s_srf              ='s_srf'
	    self.PW                 ='PW'
	    self.lw_up_srf          ='lw_up_srf'
	    self.lw_dn_srf          ='lw_dn_srf'
	    self.sw_up_srf          ='sw_up_srf'
	    self.sw_dn_srf          ='sw_dn_srf'
	    self.Tg                 ='Tg'
	    self.Ptend              ='Ptend'
	    self.T                  ='T'
	    self.q                  ='q'
	    self.u                  ='u'
	    self.v                  ='v'
	    self.omega              ='omega'
	    self.div                ='div'
	    self.divT               ='divT'
	    self.vertdivT           ='vertdivT'
	    self.divq               ='divq'
	    self.vertdivq           ='vertdivq'
	    self.s                  ='s'
	    self.divs               ='divs'
	    self.vertdivs           ='vertdivs'
	    self.dsdt               ='dsdt'
	    self.dTdt               ='dTdt'
	    self.dqdt               ='dqdt'
	    self.q1                 ='q1'
	    self.q2                 ='q2'

    def __iter__(self):
        for each in self.__dict__.keys():
            yield self.__getattribute__(each)


#To assint the label of the family of 
#variables

def ncload(file_l,calendar):

    label=variables() 

    # Your filename
    nc_file    = '%s'%(file_l)  

    # Dataset is the class behavior to open the file
    # and create an instance of the ncCDF4 class
    nc_v = Dataset(nc_file,'r+')

    label.bdate              =nc_v.variables['bdate']
    label.lev                =nc_v.variables['lev']
    label.lon                =nc_v.variables['lon']
    label.lat                =nc_v.variables['lat']
    label.tsec               =nc_v.variables['tsec']
    label.phis               =nc_v.variables['phis']
    label.prec_srf           =nc_v.variables['prec_srf']
    label.lhflx              =nc_v.variables['lhflx']
    label.shflx              =nc_v.variables['shflx']
    label.Ps                 =nc_v.variables['Ps']
    label.Tsair              =nc_v.variables['Tsair']
    label.RH_srf             =nc_v.variables['RH_srf']
    label.wspd_srf           =nc_v.variables['wspd_srf']
    label.usrf               =nc_v.variables['usrf']
    label.vsrf               =nc_v.variables['vsrf']
    label.rad_net_srf        =nc_v.variables['rad_net_srf']
    label.lw_net_toa         =nc_v.variables['lw_net_toa']
    label.sw_net_toa         =nc_v.variables['sw_net_toa']
    label.sw_dn_toa          =nc_v.variables['sw_dn_toa']
    label.cld_low            =nc_v.variables['cld_low']
    label.cld_mid            =nc_v.variables['cld_mid']
    label.cld_high           =nc_v.variables['cld_high']
    label.cld_tot            =nc_v.variables['cld_tot']
    label.cld_thick          =nc_v.variables['cld_thick']
    label.cld_top            =nc_v.variables['cld_top']
    label.LWP                =nc_v.variables['LWP']
    label.dh2odt_col         =nc_v.variables['dh2odt_col']
    label.h2o_adv_col        =nc_v.variables['h2o_adv_col']
    label.evap_srf           =nc_v.variables['evap_srf']
    label.dsdt_col           =nc_v.variables['dsdt_col']
    label.s_adv_col          =nc_v.variables['s_adv_col']
    label.rad_heat_col       =nc_v.variables['rad_heat_col']
    label.qs                 =nc_v.variables['qs']
    label.s_srf              =nc_v.variables['s_srf']
    label.PW                 =nc_v.variables['PW']
    label.lw_up_srf          =nc_v.variables['lw_up_srf']
    label.lw_dn_srf          =nc_v.variables['lw_dn_srf']
    label.sw_up_srf          =nc_v.variables['sw_up_srf']
    label.sw_dn_srf          =nc_v.variables['sw_dn_srf']
    label.Tg                 =nc_v.variables['Tg']
    label.Ptend              =nc_v.variables['Ptend']
    label.T                  =nc_v.variables['T']
    label.q                  =nc_v.variables['q']
    label.u                  =nc_v.variables['u']
    label.v                  =nc_v.variables['v']
    label.omega              =nc_v.variables['omega']
    label.div                =nc_v.variables['div']
    label.divT               =nc_v.variables['divT']
    label.vertdivT           =nc_v.variables['vertdivT']
    label.divq               =nc_v.variables['divq']
    label.vertdivq           =nc_v.variables['vertdivq']
    label.s                  =nc_v.variables['s']
    label.divs               =nc_v.variables['divs']
    label.vertdivs           =nc_v.variables['vertdivs']
    label.dsdt               =nc_v.variables['dsdt']
    label.dTdt               =nc_v.variables['dTdt']
    label.dqdt               =nc_v.variables['dqdt']
    label.q1                 =nc_v.variables['q1']
    label.q2                 =nc_v.variables['q2']


    tu = calendar[0]#"seconds since 2014-2-15 0:00:00 0:00"
    tc = calendar[1]#"gregorian"

    #label.data       = num2date(label.tsec[:],units=tu,calendar=tc, only_use_cftime_datetimes=False)
    label.data       = cft.num2date(label.tsec[:],units=tu,calendar=tc, only_use_cftime_datetimes=True)
    #label.data.astype("datetime64[ns]")

    #solution without used Datetimenoleap, and using  cftime
    #label.data       = cft.num2pydate(label.tsec[:],units=tu,calendar=tc)

    return label 
