################################################################
#Program to plot meteorological date
#of the OUT_STAT files of SAM with
#python and  Netcdf library.

###########################
#Modified:23/01/23
# To run using the same files
# python_src files

###########################
#Create by: Jhonatan Aguirre
#Date:07/04/2022
#working
#python 3.9

################################################################
# To activate this environment, use
#
# $ conda activate py37
# Went panda is use
# To deactivate an active environment, use
#
# $ conda deactivate
#
# path of the ncfile
#
from    Parameters import *

################################################################ 
# to defined fig out direction 
from    files_direction import *

# Load function to make the diurnal cycle and profiles figures.  
import  sam_python.diurnal   as dc

import  sam_python.temporal_plot   as tp


#separate with colon

#ex           =  [f1['data']] 
#
#di           =  (f1['date'][0],f1['date'][1],f1['date'][2],9,f1['date'][0],f1['date'][1],f1['date'][2],18)
#
#name        =  [f1['name']]

ex          =  [goa1  ] 
name        =  ['iop1']
di          =  (2014,2,15,0,2014,3,15,18)



ex_date    =   [
                    di, di, di, 
                    di, di, di, 
                    di# di, di, 
               ]

var         =  [
                   ex[0].MCUP  , ex[0].CLD ,ex[0].TKE , 
                   ex[0].TVFLUX, ex[0].WOBS,ex[0].RELH,  
                   ex[0].QC    #, ex[0].Q1C ,ex[0].Q2,  
               ]


lim         =  [
                (0,0.1),(0,15),(0,3),
                (-0.7E-3,5.5E-3),(-1.5,2.0),(20,90),
                (0, 0.04),(-20,20),(-80,90),
               ]

var_to      =  [
                    1, 100, 1,  
                    9.81/(ex[0].THETA[10,10]*1005.0*ex[0].RHO[10,10]), 100, 1,  
                    1, 1 , 1,                    
               ]

#figures name
exp_label   =  [ 
                 'uMF', 'CF'  , 'TKE',
                 'B'  , 'Wobs', 'RH',
                 'CLW', 'Q1'  , 'Q2',
               ]

l1           =  ( 0.00,4.0,'lower right',True,True)
lrh          =  ( 70.0,4.0,'lower right',True,True)

leg_loc      =  [ 
                l1,l1,l1,
                l1,l1,lrh,
                l1,l1,l1,
                ]

alt         =  [
                5.0, 5.0, 5.0, 
                5.0, 5.0, 5.0, 
                5.0, 5.0, 5.0, 
               ]

color       =  [
                'black','black','black',
                'black','black','black',
                'black','black','black',
               ]

show       =  [
                True,True,True,
                True,True,True,
                True,True,True,
               ]

diurnal     =  [
                True,True,True,
                True,True,True,
                True,True,True,
              ]

#Minimum calll
#fig,ax      =   dc.diurnal_hours(ex,name,var,ex_date)
#complete call

dc.diurnal_hours(ex,name,var,ex_date,alt, lim, var_to,color,exp_label,leg_loc,diurnal, show)
