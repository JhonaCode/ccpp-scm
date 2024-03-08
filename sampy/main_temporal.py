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
#import  sam_python.diurnal   as dc

import  sam_python.temporal_plot   as tp

import  datetime as dt 

#separate with colon

#print(goa1.pres)
#print(goa1.vert_dim_layer)

exp         =  [goa1,goa2] 

c1          =  ['green','magenta','blue']

color       =  [
                c1,c1,c1
               ]
lim         =  [
                [(0,160,10),(0,440,10),(0,0.6,10)],
                [(0,160,10),(0,440,10),(0,0.6,10)],
                [(0,160,10),(0,440,10),(0,0.6,10)],
               ]

var=['shf','lhf']

#X           =  'Hours LT (UTC-4)'
#Y           =  'Sensible Heat FLux $\mathrm{ [W m^{-2}]}$'
#Y2          =  'Latent   Heat FLux $\mathrm{ [W m^{-2}]}$'
#Y3          =  'Precipitation  [mm]'
X           =  'Horas HL (UTC-4)'
Y           =  'Fluxo de Calor Sensible$\mathrm{ [W m^{-2}]}$'
Y2          =  'Fluxo de Calor Latente $\mathrm{ [W m^{-2}]}$'
Y3          =  'Precipitação [mm]'
aformat     =  ['default']
#aformat     =  ['Day','%d']

l1          =  ([X,Y], ['a)',(2014,2,15,0),140],[False,'upper right'] ,[0.75,0],aformat)
l2          =  ([X,Y2],['b)',(2014,2,15,0),400],[False,'upper right'] ,[0.75,0],aformat)
l3          =  ([X,Y3],['c)',(2014,2,15,0),0.55],[False,'upper right'],[0.75,0],aformat)

l4          =  ([X,Y], ['a)',(2014,9,1,0),140],[False,'upper right']  ,[0.75,0],aformat)
l5          =  ([X,Y2],['b)',(2014,9,1,0),400],[False,'upper right']  ,[0.75,0],aformat)
l6          =  ([X,Y3],['c)',(2014,9,1,0),0.55],[False,'upper right'] ,[0.75,0],aformat)


plot_def    =  [ 
                [l1,l2,l3],
                [l4,l5,l6],
                ]


#To transform the units of the variable to plot 
var_to      =  [
                  [1,1,1000],  [1,1,1000] ,  [1,1,1000]  
               ]

#figures name
exp_label   =  [ 
                 ['SHF_IOP1','LHF_IOP1','PREC_IOP1'],
                 ['SHF_IOP2','LHF_IOP2','PREC_IOP2'],
                 ['SHF_IOP2','LHF_IOP2','PREC_IOP2'],
               ]



show        =  [
                [True,True,True],
                [True,True,True],
                [True,True,True],
                [False,False,False]
                #True, True
               ]


#for ex in exp:
#
#    tp.temporal_plot_dict(ex,lim=lim,var_to=var_to,exp_label=exp_label,plot_def=plot_def,color=color)

#tp.temporal_plot_exp(exp,var_to=var_to,lim=lim,exp_label=exp_label,plot_def=plot_def,color=color,show=True)

tp.temporal_plot_exp_var(exp,var=var,var_to=var_to,explabel1=exp_label,plot_def=plot_def,alt=lim,color=color,show=show)





