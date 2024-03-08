######################################
######################################
#Program to plot meteorological date 
#of the OUT_STAT files 
#of SAM with python and  Netcdf library.
#
#Create by: Jhonatan Aguirre 
#Date:23/03/1010
#working 
######################################
######################################

import  sam_python.two_dimensional   as cf

from    files_direction import *

from    Parameters import *

import  matplotlib.pyplot as plt

#Defined the name of the experiments defined in Parameters files 
#separate with colon

exp          =  [goa1] 

contour      =[
                (240,310,21),
                (0.000001,0.001,51),
                (-0.2,0.2,21),
              ]


var_to      =  [
                 1, 1, 1,  
               ]

#figures name
exp_label   =  [ 
                 'T', 'Q_l','Wobs'
               ]

l1           =  ( 0.00,14.0)
l2           =  ( 0.00,14.0)
l3           =  ( 0.00,14.0)

leg_loc      =  [ 
                l1,l1,l1,
                ]

alt         =  [
                15.0, 15.0, 15.0, 
               ]

color       =  [
                'RdBu_r','RdBu_r','RdBu_r','RdBu_r'
               ]

show       =  [
                True,True,True,
                True,True,True,
                True,True,True,
               ]

#Minimum calll
#fig,ax      =   dc.diurnal_hours_dict(ex,name,var,ex_date)
#complete call

cf.plot2d_contour_ccpp(exp,var_to=var_to,explabel=exp_label,leg_loc=leg_loc,contour=contour,color=color, alt=alt,show=show)
