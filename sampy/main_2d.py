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


import  sam_python.two_dimensional   as two


import  datetime as dt 

#separate with colon

#print(goa1.pres)
#print(goa1.vert_dim_layer)

exp         =  [goa1,goa2] 


explabel1     = [
                 ['a)FMA IOP1',r'a)T IOP1','a)q$_l$ IOP1'],
                 ['b)FMA IOP2',r'b)T IOP2','b)q$_l$ IOP2'],
                ]



cmf=(00.0,100.0,21,r'kg/m/s$^{2}$')
ccf=(0.0,50,21,'%')
crh=(0,100,21,'%')
cw =(-3,3,21,'cm/s')
#cq =(0,16,21,'g/kg')
cq =(0.0,1.0,21,'g/kg')

cT =(190,310,21,'K')


c1=[cmf,cq,cT]

contour      =[
                c1,c1,c1,c1,c1,
                c1,c1,c1,c1,c1,
                c1,c1,c1,c1,c1,
                c1,c1,c1,c1,c1,

                ]


var=['upd_mf','ql','T',]



a0=(0,15.0)
a1=[a0,a0,a0,a0]
alt         =  [
               a1,a1,a1,a1,a1,
               a1,a1,a1,a1,a1,
               a1,a1,a1,a1,a1,
               a1,a1,a1,a1,a1,
               ]

#co1='RdBu_r'
co1         =  ['whbuyl','whbuyl','whbuyl','RdBu_r']

color       =  [
                co1,co1,co1,co1,co1,co1,
                co1,co1,co1,co1,co1,co1,
                co1,co1,co1,co1,co1,co1,
                co1,co1,co1,co1,co1,co1,
               ]




#To transform the units of the variable to plot 
var_to      =  [
                  [1,1000,1],  [1,1000,1] ,  [1,1000,1]  
               ]

s1=[True,True,True,False,]
show       =  [ s1,s1,s1,s1,s1,s1, 
                s1,s1,s1,s1,s1,s1,
                s1,s1,s1,s1,s1,s1,
               ]


                #barra,x,y,lfc..
a0          =  (False,True,False,False,0.35,0.00,['default'])
ab          =  (True ,True,False,False,0.35,1.28,['default'])
abz         =  (True ,True,False,False,0.35,1.34,['default'])
az          =  (False,True,False,False,0.35,0.06,['default'])

a1          =  [a0,a0,a0,a0]
a2          =  [az,az,az,az]
a3          =  [ab,ab,ab,ab]
a4          =  [abz,abz,abz,abz]
axis_on     =  [
                a4,
                a4,
                a4,
                a4,
               ]




two.plot2d_contour_ccpp(exp,var=var,contour=contour,explabel1=explabel1,alt=alt,color=color,var_to=var_to,axis_on=axis_on,show=show)





