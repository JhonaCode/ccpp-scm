##**************************
#Program to created the forcing files 
#to SAM using the single column model forcing. 
##################################################
#The data from the SCM program was used as forcing. 
##################################################
#Note: The FORCING ARE CRETED IN UTC HOUR
#Note: The plots are in local Manuaus time UTC-4
#Data:20-04-22
#Created by: Jhonatan Aguirre
#**************************

import  os,sys
import  inspect

import  numpy as np
import  matplotlib as mpl
import  matplotlib.pyplot as plt

#       To change the plot parameter 
from   	plotparameters 		import 	*

#To load the variable  in the datafile
from    variabletoload_go   import *

import  datetime as dt

from    Parameters_forcing_iop import *

#import   metpy.calc as mpy 
import   metpy.calc

from    metpy.units import units

#To find data in the array
import  data_own        as down

#Calculate necessary varibles to use as 
#foncing in SAM, usint METPY
from    calculate_variables_forcing import *

#Mean of days with diurnal profiles
import  diurnal_days as d_days 

#hours colors
import  hour_color as h_color 

today = dt.datetime.today()


folder='scm_lsf_experiments'

ano=2014
#exp     =   scm_iop1
#label     =   'iop1_scm_sam' 
#exp     =   scm_iop2
#label     =   'iop2_scm_sam' 
exp       =   scm_2014
label     =   'scm_2014' 

#hour
h=[0,19]
#days
days=[1,31]
#month
month=[1,12]


idi     = dt.datetime(ano, month[0] ,days[0], h[0])
idf     = dt.datetime(ano, month[1] ,days[1], h[1])

ni,nf= down.data_n(idi,idf,exp.data[:])

dayinit = exp.tsec[ni]
dayfin  = exp.tsec[nf]


tsec= exp.tsec
data= exp.data
lev = exp.lev[::-1]

ndtp   = len(tsec)
ndlev  = len(lev)

##########################
#To transform time to days since date_inicial(2014/1/1)
#format to SAM

#date_inicial=dt.datetime(2014,1,1,0)
#Equal to used in the iop2 forcing 
date_inicial=dt.datetime(2013,12,31,0)
timedelta= data[:]-date_inicial

time = np.zeros(len(tsec))

j=0
for i in timedelta:
    time[j]=i.total_seconds()/86400.0
    j=j+1
##########################

T,theta,q,t_ls,q_ls,s_ls,u,v,w,omega,T_skin,z,pressure,SH,LH=variables_metpy_scm(exp)


#Out files 
file1  = open('../../snd_%s'%(label),'w+')
#Static energy h
file2  = open('../../lsf_%s_h'%(label),'w+')
#potencial temperature p
file21 = open('../../lsf_%s_p'%(label),'w+')
file3  = open('../../sfc_%s'%(label),'w+')

file11= open('%s/snd_%s'  %(folder,label),  'w+')
file22= open('%s/lsf_%s_h'%(folder,label),'w+')
file23= open('%s/lsf_%s_p'%(folder,label),'w+')
file33= open('%s/sfc_%s'  %(folder,label)  ,'w+')

file1.write( "z[m]\tp[mb]\ttp[K]\tq[g/kg]\tu[m/s]\tv[m/s]\n")
file2.write( "z[m]\tp[mb]\tsls(h)[K/s] \tqls[kg/kg/s]\tu[m/s]\tv[m/s]\tw[m/s]\n")
file21.write("z[m]\tp[mb]\ttls[K/s]    \tqls[kg/kg/s]\tu[m/s]\tv[m/s]\tw[m/s]\n")
file3.write( "day\t\tsst(K)\t\tH(W/m2)\t\tLE(W/m2)\tTAU(m2/s2) \n")

file11.write("z[m]\tp[mb]\ttp[K]\tq[g/kg]\tu[m/s]\tv[m/s]\n")
file22.write("z[m]\tp[mb]\tsls(h)[K/s] \tqls[kg/kg/s]\tu[m/s]\tv[m/s]\tw[m/s]\n")
file23.write("z[m]\tp[mb]\ttls(h)[K/s] \tqls[kg/kg/s]\tu[m/s]\tv[m/s]\tw[m/s]\n")
file33.write("day\t\tsst(K)\t\tH(W/m2)\t\tLE(W/m2)\tTAU(m2/s2) \n")



for k in range(0,ndtp): 

    #if(time[k]>=dayinit and time[k]<=dayfin):

    file1.write("%f\t%d\t%f\tday\tlevels\tpres0\n"%( time[k],ndlev,pressure[0]))
    file11.write("%f\t%d\t%f\tday\tlevels\tpres0\n"%(time[k],ndlev,pressure[0]))

    file2.write( "%f\t%d\t%f\tday\tlevels\tpres0\n"%(time[k],ndlev,pressure[0]))
    file21.write("%f\t%d\t%f\tday\tlevels\tpres0\n"%(time[k],ndlev,pressure[0]))
    file22.write("%f\t%d\t%f\tday\tlevels\tpres0\n"%(time[k],ndlev,pressure[0]))
    file23.write("%f\t%d\t%f\tday\tlevels\tpres0\n"%(time[k],ndlev,pressure[0]))

    Tau         =   velw_sfc_flx(u[k,0],v[k,0],w[k,0]) 

    file3.write("%f\t%f\t%f\t%f\t%f\n"%( time[k],T_skin[k],SH[k],LH[k],Tau[0])) 
    file33.write("%f\t%f\t%f\t%f\t%f\n"%(time[k],T_skin[k],SH[k],LH[k],Tau[0])) 

    for i in range(0,ndlev):

        file1.write("%f\t%f\t%f\t%f\t%f\t%f\n"%(z[0,i],pressure[i],theta[k,i],q[k,i],u[k,i],v[k,i])) 
        file11.write("%f\t%f\t%f\t%f\t%f\t%f\n"%(z[0,i],pressure[i],theta[k,i],q[k,i],u[k,i],v[k,i])) 

        file2.write( "%f\t%f\t%e\t%e\t%f\t%f\t%f\n"%(z[0,i],pressure[i],s_ls[k,i],q_ls[k,i],u[k,i],v[k,i],w[k,i])) 
        file21.write("%f\t%f\t%e\t%e\t%f\t%f\t%f\n"%(z[0,i],pressure[i],t_ls[k,i],q_ls[k,i],u[k,i],v[k,i],w[k,i])) 
        file22.write("%f\t%f\t%e\t%e\t%f\t%f\t%f\n"%(z[0,i],pressure[i],s_ls[k,i],q_ls[k,i],u[k,i],v[k,i],w[k,i])) 
        file23.write("%f\t%f\t%e\t%e\t%f\t%f\t%f\n"%(z[0,i],pressure[i],t_ls[k,i],q_ls[k,i],u[k,i],v[k,i],w[k,i])) 

file1.close()
file2.close()
file21.close()
file3.close()

file11.close()
file22.close()
file23.close()
file33.close()


plt.close()


print()

exit()

## Writing NetCDF files to scam 
#
## Open a new NetCDF file to write the data to. For format, you can choose from
## 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'

nc_out              = Dataset('goamazon_iop2_scm_forcing.nc', 'w', format='NETCDF4')
nc_out.description  = " FORCING TO REPRESENT A COMPOSITE OF  \
PURE SHALLOW CONVECTION COLECTED IN THE OBSERTATIONAL CAMPAIGN \
GOAMAZON THROUGHT THE IOP1 AND IOP2. THE DATA IS REPRESENTED IN THE FIRST DAY OF 2014"


nc_out.history = "Created " + today.strftime("%d/%m/%y")
nc_out.author = "Jhonatan A. A manco"


time_nc = nc_out.createVariable('bdate', 'i4',())
nc_out.variables['bdate'][:] =  20140101 

# we can create the new dimension
## Even though we know the size, we are going to set the size to unknown
nc_out.createDimension('lon', 1)
nc_out.createDimension('lat', 1)

dtime =len(time_sec)
dlev=len(lev)
nc_out.createDimension('lev',dlev )
nc_out.createDimension('time', dtime)

##########################################
#time_nc = nc_out.createVariable('time', 'f4',('time',))
#
### Assign the dimension data to the new NetCDF file.
#nc_out.variables['time'][:] = time[:]


tsec_nc = nc_out.createVariable('tsec', 'f4',('time',))
tsec_nc.setncatts({'long_name': u"Time in seconds after 0Z on base date",\
                    'units': u"s",\
                    'statistic': u'none'})


nc_out.variables['tsec'][:] = time_sec[:]

##############################

lon_nc = nc_out.createVariable('lon', 'f4', ('lon'))
lon_nc.setncatts({'long_name': u"Longitude",\
                    'units': u"degrees",\
                    'var_desc': u"Longitude of ARM SITE MANACAPURU AMAZONAS",\
                    'statistic': u'none'})

nc_out.variables['lon'][:] = -60.00

##############################

lat_nc = nc_out.createVariable('lat', 'f4', ('lat'))
lat_nc.setncatts({'long_name': u"Latitude",\
                    'units': u"degrees",\
                    'var_desc': u"latidude of ARM SITE MANACAPURU AMAZONAS",\
                    'statistic': u'none'})

nc_out.variables['lat'][:] = -3.15


##############################

p_nc = nc_out.createVariable('lev', 'f4', ('lev',))
p_nc.setncatts({'long_name': u"Pressure",\
                    'units': u"Pa", 'level_desc': u'Surface',\
                    'var_desc': u"pressure from varanal analysis",\
                    'statistic': u'none'})

nc_out.variables['lev'][:] = pa[::-1].to('Pa').magnitude
#nc_out.variables['lev'][:] = pa[:].to('Pa')


##############################

sh_nc = nc_out.createVariable('shflx', 'f4', ('time','lat','lon'))
sh_nc.setncatts({'long_name': u"Surface sensible heat flux",\
                    'units': u"W/m2", 'level_desc': u'Surface',\
                    'var_desc': u"Surface Sensible heat flux from varanal analysis",\
                    'statistic': u'none\nM'})


sh_3d=np.zeros([SH_mean.shape[0],1,1])

sh_3d[:,0,0]=SH_mean[:]


nc_out.variables['shflx'][:] = sh_3d[:,0,0]


##############################
lh_nc = nc_out.createVariable('lhflx', 'f4', ('time','lat','lon'))
lh_nc.setncatts({'long_name': u"Surface latent heat flux",\
                    'units': u"W/m2", 'level_desc': u'Surface',\
                    'var_desc': u"Surface Latent  heat flux from varanal analysis",\
                    'statistic': u'none\nM'})
#nc_out.variables['lhflx'][:] = LH_mean[:]

lh_3d=np.zeros([LH_mean.shape[0],1,1])

lh_3d[:,0,0]=LH_mean[:]

nc_out.variables['lhflx'][:] = lh_3d[:]

##############################
phis_nc = nc_out.createVariable('phis', 'f4', ('lon','lat',))
phis_nc.setncatts({'long_name': u"Surface geopotential",\
                    'units': u"m2/s2", 'level_desc': u'Surface',\
                    'var_desc': u"Surface geopotential bomex has ",\
                    'statistic': u'none\nM'})

phis_3d=np.zeros([1,1])

phis_3d[0,0]=0.0

nc_out.variables['phis'][:] = phis_3d[:]

################################
ustar_nc = nc_out.createVariable('ustar', 'f4', ('time','lat','lon'))

ustar_nc.setncatts({'long_name': u"Surface frictional velocity",\
                    'units': u"m/s",\
                    'var_desc': u"Surface frictional velocity, bomex has",\
                    'statistic': u'none'})

ustar_3d=np.zeros([theta_inv.shape[0],1,1])

#bomex value
ustar_3d[:,0,0]=0.28

nc_out.variables['ustar'][:] = ustar_3d[:]


##############################

Tg_nc = nc_out.createVariable('Tg', 'f4', ('time','lat','lon'))

Tg_nc.setncatts({'long_name': u"Ground Temperature",\
                    'units': u"K",\
                    'var_desc': u"Ground Temperature from varanal analysis",\
                    'statistic': u'none'})

Tg_3d=np.zeros([theta_inv.shape[0],1,1])

Tg_3d[:,0,0]=T_skin_mean[:]

nc_out.variables['Tg'][:] = Tg_3d[:]

################################
##

Ps_nc = nc_out.createVariable('Ps', 'f4', ('time','lat','lon',))

Ps_nc.setncatts({'long_name': u"Surface pressure",\
                    'units': u"Pa",\
                    'var_desc': u"Surface pressure from varanal analysis",\
                    'statistic': u'none'})

Ps_3d=np.zeros([theta_inv.shape[0],1,1])

Ps_3d[:,0,0]=pa[0].to('Pa')

nc_out.variables['Ps'][:] = Ps_3d[:]

################################
##

Pt_nc = nc_out.createVariable('Ptend', 'f4', ('time','lat','lon',))

Pt_nc.setncatts({'long_name': u"Surface pressure tendency",\
                    'units': u"Pa/s",\
                    'var_desc': u"Surface pressure tendency, bomex has",\
                    'statistic': u'none'})

Pt_3d=np.zeros([theta_inv.shape[0],1,1])

Pt_3d[:,0,0]= 0.0

nc_out.variables['Ptend'][:] = Pt_3d[:]



##############################

T_nc = nc_out.createVariable('T', 'f4', ('time','lev','lat','lon'))
T_nc.setncatts({'long_name': u"Air Temperature",\
                    'units': u"K", 'level_desc': u'Surface',\
                    'var_desc': u"Air Temperature from varanal analysis",\
                    'statistic': u'none\nM'})

theta_4d=np.zeros([theta_inv.shape[0],theta_inv.shape[1],1,1])

#theta_4d[:,:,0,0]=theta[:,:]
theta_4d[:,:,0,0]=T_mean[:,::-1]

nc_out.variables['T'][:] = theta_4d[:]

##############################
q_nc = nc_out.createVariable('q', 'f4', ('time','lev','lat','lon'))
q_nc.setncatts({'long_name': u"Water vapor mixing ratio",\
                    'units': u"kg/kg", 'level_desc': u'Surface',\
                    'var_desc': u"Water vapor mixin ratio from varanal analysis",\
                    'statistic': u'none\nM'})

q_4d=np.zeros([theta_inv.shape[0],theta_inv.shape[1],1,1])

q_4d[:,:,0,0]=q[:,::-1].to('kg/kg')

nc_out.variables['q'][:] = q_4d[:]

##############################

u_nc = nc_out.createVariable('u', 'f4', ('time','lev','lat','lon'))
u_nc.setncatts({'long_name': u"u wind velocity ",\
                    'units': u"m/s", 'level_desc': u'Surface',\
                    'var_desc': u"u zonal  wind velocity  from varanal analysis",\
                    'statistic': u'none\nM'})

u_4d=np.zeros([theta_inv.shape[0],theta_inv.shape[1],1,1])

u_4d[:,:,0,0]=u[:,::-1]

nc_out.variables['u'][:] = u_4d[:]

##############################

v_nc = nc_out.createVariable('v', 'f4', ('time','lev','lat','lon'))
v_nc.setncatts({'long_name': u"v wind velocity ",\
                    'units': u"m/s", 'level_desc': u'Surface',\
                    'var_desc': u"v meridional  wind velocity  from varanal analysis",\
                    'statistic': u'none\nM'})

v_4d=np.zeros([theta_inv.shape[0],theta_inv.shape[1],1,1])

v_4d[:,:,0,0]=v[:,::-1]

nc_out.variables['v'][:] = v_4d[:]

##############################

o_nc = nc_out.createVariable('omega', 'f4', ('time','lev','lat','lon'))
o_nc.setncatts({'long_name': u"omega",\
                    'units': u"Pa/s", 'level_desc': u'Surface',\
                    'var_desc': u" vertical pressure velocity  from varanal analysis",\
                    'statistic': u'none\nM'})

o_4d=np.zeros([theta_inv.shape[0],theta_inv.shape[1],1,1])

o_4d[:,:,0,0]=omega_mean[:,::-1].to('Pa/s')

nc_out.variables['omega'][:] = o_4d[:]

##############################

tls_nc = nc_out.createVariable('divT', 'f4', ('time','lev','lat','lon'))
tls_nc.setncatts({'long_name': u"Large scale temperature forcing",\
                      'units': u"K/s", 'level_desc': u'Surface',\
                      'var_desc': u" Large scale temperature forcing from varanal analysis",\
                      'statistic': u'none\nM'})
   
tls_4d=np.zeros([theta_inv.shape[0],theta_inv.shape[1],1,1])

tls_4d[:,:,0,0]=t_ls[:,::-1]

nc_out.variables['divT'][:] = tls_4d[:]

##############################

qls_nc = nc_out.createVariable('divq', 'f4', ('time','lev','lat','lon'))
qls_nc.setncatts({'long_name': u" Large scale water vapor mixing ratio forcing",\
                      'units': u"kg/kg/s", 'level_desc': u'Surface',\
                      'var_desc': u" Large scale water vapor mixing ratio forcing from varanal analysis",\
                      'statistic': u'none\nM'})
   
qls_4d=np.zeros([theta_inv.shape[0],theta_inv.shape[1],1,1])

qls_4d[:,:,0,0]=q_ls[:,::-1]

nc_out.variables['divq'][:] = qls_4d[:]

##############################



nc_out.close()  # close the new file
print('Files were created' )
print('lsf_%s'%(explabel[0]) )
print('snd_%s'%(explabel[0]) )
print('sfd_%s'%(explabel[0]) )
print('%s.nc'%(explabel[0]) )


exit()

#used user parameter to plot(plotparameter.py
mpl.rcParams.update(params3)
###########################################3

fig   = plt.figure()
ax1   = plt.axes()

color = 'tab:red'
color2 = 'tab:green'
ax1.set_xlabel(r'$\mathrm{\theta_{inv}\ [K]}$', color=color)
ax1.set_ylabel('z [km]')

ax2 = ax1.twiny()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'
ax2.set_xlabel(r'q [g/kg]', color=color)

ax1.axis([296,320,0,5])
ax2.axis([3,18.0,0,5])

#fig.tight_layout()  # otherwise the right y-label is slightly clipped



for k in range(0,len(time)): 

        for i in range(0,len(hp)): 

            if hp[i]==data[k].hour: 

                print data[k]

                line,color=h_color.color_hours(hp[i])

                #-4 MANAUS LOCAL TIME
                ax1.plot(theta_inv[k,:], z[k,:]/1000.0,color=color,label='%d:00'%(hp[i]-4),dashes=line)

                ax2.plot(q[k,:], z[k,:]/1000.0,color=color,dashes=line)

#ax1.legend(frameon=False,loc='center right')

plt.savefig('%s/forcing_mean_tq.pdf'%file_out,format='pdf',bbox_inches='tight' , dpi=1000)
################################################################3

fig   = plt.figure()
ax   = plt.axes()

#ax.legend(frameon=False)

ax.axis([-12,5,0,4])

ax.text(-10, 3, r'u $\rightarrow$ ')
ax.text( 1, 3, r'$\leftarrow$ v')

#plt.ylabel(r'z [km]') 
plt.ylabel(r'') 
plt.xlabel(r' $\mathrm{ [ms^{-1}]}$ ') 

for k in range(0,len(time)): 

        for i in range(0,len(hp)): 

            if hp[i]==data[k].hour: 

                line,color=h_color.color_hours(hp[i])

                #-4 MANAUS LOCAL TIME
                ax.plot(u[k,:] ,z[k,:]/1000.0,label='%d:00'%(hp[i]-4),dashes=line,color=color)
                ax.plot(v[k,:] ,z[k,:]/1000.0,dashes=line,color=color)

#ax.legend(frameon=False,loc='lower right')


plt.savefig('%s/forcing_mean_uv.pdf'%file_out, bbox_inches='tight' ,format='pdf', dpi=1000)
plt.close()
############################################################

fig   = plt.figure()
ax    = plt.axes()

#ax.legend(frameon=False)

ax.axis([-0.98,0.4,0,4])

#plt.ylabel(r'z [km]') 
plt.xlabel(r'w  $\mathrm{ [cms^{-1}]}$ ') 
plt.ylabel(r'') 

for k in range(0,len(time)): 

        for i in range(0,len(hp)): 

            if hp[i]==data[k].hour: 

                line,color=h_color.color_hours(hp[i])

                #-4 MANAUS LOCAL TIME
                ax.plot(w[k,:].to('cm/s') ,z[k,:]/1000.0,label='%d:00'%(hp[i]-4),dashes=line,color=color)

#ax.legend(frameon=False,loc='lower right')

plt.savefig('%s/forcing_mean_w.pdf'%file_out, bbox_inches='tight' ,format='pdf', dpi=1000)
############################################################


################################
fig   = plt.figure()
ax   = plt.axes()

color = 'black'
ax.set_ylabel('')
ax.set_xlabel(r'$\mathrm{\theta_{ls}\ [Ks^{-1}]}$', color=color)
#ax.set_ylabel('z [km]')

ax.axis([-3e-5, 2e-5,0,4])

#fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.ylabel(r'') 
plt.xlabel(r'$\mathrm{ \theta_{ls}[Ks^{-1}]}$ ') 


##to plot with scientific notation 
ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))

for k in range(0,len(time)): 

        for i in range(0,len(hp)): 

            if hp[i]==data[k].hour: 

                line,color=h_color.color_hours(hp[i])

                #-4 MANAUS LOCAL TIME
                ax.plot(t_ls[k,:] ,z[k,:]/1000.0,dashes=line,color=color)


plt.savefig('%s/forcing_mean_tls.pdf'%file_out, format='pdf',bbox_inches='tight' , dpi=1000)

################################
fig   = plt.figure()
ax   = plt.axes()

color = 'black'
ax.set_xlabel(r'$\mathrm{q_{ls} [g kg^{-1} s^{-1}]}$', color=color)
#ax.set_ylabel('z [km]')
plt.ylabel(r'') 

color = 'black'

ax.axis([-4e-8,2.5e-8,0,4])

#plt.tight_layout()  # otherwise the right y-label is slightly clipped


##to plot with scientific notation 
ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))


for k in range(0,len(time)): 

        for i in range(0,len(hp)): 

            if hp[i]==data[k].hour: 

                line,color=h_color.color_hours(hp[i])

                #-4 MANAUS LOCAL TIME
                ax.plot(q_ls[k,:] ,z[k,:]/1000.0,label='%d:00'%(hp[i]-4),dashes=line,color=color)

ax.legend(frameon=False,loc='upper right')

plt.savefig('%s/forcing_mean_qls.pdf'%file_out, format='pdf',bbox_inches='tight' , dpi=1000)



plt.show()





