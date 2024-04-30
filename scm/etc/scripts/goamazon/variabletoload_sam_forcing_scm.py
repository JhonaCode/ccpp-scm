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

    label.q2                 =nc_v.variables['q2']


    tu = calendar[0]#"seconds since 2014-2-15 0:00:00 0:00"
    tc = calendar[1]#"gregorian"

    #label.data       = num2date(label.tsec[:],units=tu,calendar=tc, only_use_cftime_datetimes=False)
    label.data       = cft.num2date(label.tsec[:],units=tu,calendar=tc, only_use_cftime_datetimes=True)
    #label.data.astype("datetime64[ns]")

    #solution without used Datetimenoleap, and using  cftime
    #label.data       = cft.num2pydate(label.tsec[:],units=tu,calendar=tc)

    return label 
