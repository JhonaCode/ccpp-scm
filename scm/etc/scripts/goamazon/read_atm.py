import numpy as np
import scipy.interpolate


def call_read_atm(file,levels,T):

    atm  =read_atm(file)

    #print(atm.keys())
    
    #nc_fid_o3 = Dataset("../../data/raw_case_input/mid_lat_summer_std.nc", 'r')
    #oz_pres = nc_fid_o3.variables['pressure']
    #oz_data = nc_fid_o3.variables['o3']

    #convert levels in mb to Pa *100
    oz_pres = atm['pre'][:]*100 # to transform em Pa
    oz_data = atm['o3'][:]
    
    oz_f = scipy.interpolate.interp1d(oz_pres, oz_data)
    #o3 = oz_f(levels[1:])
    o3 = oz_f(levels[0:].tolist())


    #https://en.wikipedia.org/wiki/Air_pollutant_concentrations
    #ppmv to kg/kg
    #mg/m^3=ppmv*Mw/(0.082057338*T[k])
    for i in range(0,levels.size):
    
        o3[i]=o3[i]*48/(1e6*2.14*(T[i,0])*0.082057338)
       
    return o3
    

def read_atm ( filename ):
  """ Return contents of RFM .atm file as dictionary
  
  VERSION
    18FEB23 AD Original

  PARAMETERS
    filename Str : Name of .atm file

  RETURNS
    dictionary : { 'nlev': Int,      No.of profile levels
                   'prf1': Flt(nlev) Numpy array of nlev prf1 values
                   'prf2': Flt(nlev)   "       "         prf2 
                   etc. } 

  DESCRIPTION
    The dictionary always contains a key 'nlev' with a value of the 
    number of levels in all the profiles (necessarily the same in an 
    .atm file), and then a variable set of keys constructed from the 
    lower-case labels of each profile with a value consisting of a 
    numpy array of profile values
    See http://eodg.atm.ox.ac.uk/RFM/sum/atmfil.html for file format

  USAGE
    For example
      atm  = read_atm ( 'day.atm' )
      print(atm['nlev'])    # print size of profiles
      print(atm.keys())     # print list of supplied profiles in file
      print(atm['so2'])     # print array of SO2 profile values
  """

  with open(filename) as f:
    rec = '!'
    #print('k')
    while rec[0] == '!': rec = f.readline()  # skip initial comments
    flds = rec.split()
    print(flds[:])
    #nlev = int(flds[0])
    nlev = flds[0]
    #print('k',nlev)
    atm = { 'nlev':nlev } 
    #print(nlev)
    rec = f.readline()
    i=0
    while rec[0:4].lower() != '*end':  # repeat until end marker record
           if rec[0] == '!': continue 
           flds = rec.split()
           #print('m',1, rec)
           #print('m',2, flds[0][0])
           #print('m',3, flds[0][1])
           #print('m',4, flds[0][1:])
           key = flds[0][1:].lower() # remove '*' and change to lower case
           prf = np.fromfile(f,sep=",",count=int(nlev))
           atm[key] = prf
           rec = f.readline()
           if rec == '': break   # also exit if end-of-file without *END
           i+=1
  return atm

###nc_fid_o3 = Dataset("../../data/raw_case_input/mid_lat_summer_std.nc", 'r')
###
###oz_pres = nc_fid_o3.variables['pressure']
###oz_data = nc_fid_o3.variables['o3']
###
###oz_f = scipy.interpolate.interp1d(oz_pres, oz_data)
####o3 = oz_f(levels[1:])
###o3_2 = oz_f(levels[2:].magnitude.tolist())
###
###matplotlib.use("TkAgg")
###
###fig   = plt.figure()
###ax    = plt.axes()
###
###ax.plot(o3[2:]  ,levels[2:],color='red' )
###ax.plot(o3_2,levels[2:],color='blue')
###
###plt.show()
###exit()
###
###exit()


