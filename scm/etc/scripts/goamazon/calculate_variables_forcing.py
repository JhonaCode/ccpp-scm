
import  numpy as np
#import  matplotlib as mpl
#import  matplotlib.pyplot as plt

#import   metpy.calc as mpy 
import   metpy.calc

from    metpy.units import units
#       To change the plot parameter 

import forcing_file_common as  ffc

def variables_metpy_sam_scm_ccpp(exp):


    #Para listar 
    #print(list(exp.keys()))

    ##bdate', 'tsec', 'shflx', 'lhflx', 'phis', 'ustar', 'Tg', 'Ps', 'Ptend', 'T', 'q', 'u', 'v', 'omega', 'divT', 'divq'

    #Lengh of the time array to search
    ndtp    = len(exp.tsec) 
    #Lengh of the time array to search
    ndlev   = len(exp.lev) 

    #Mixing ratio
    #divq    #"kg/kg/s""Horizontal q advection"
    mr_ls_u     =   units.Quantity(exp.divq[:,::-1,0,0],'kg/kg/s')

    #vertdivq#"Vertical q advection""kg/kg/s" ;
    v_adv_mr_u  =   units.Quantity(exp.vertdivq[:,::-1,0,0],'kg/kg/s')

    #T tendency "K/s"
    dmrdt_u     =   units.Quantity(exp.dqdt[:,::-1,0,0],'kg/kg/s')

    mr_u        =   units.Quantity(exp.q[:,::-1,0,0],'kg/kg')


    #convert to specific humidity from mixing ratio
    q_u           =   mr_u/(1.0 + mr_u) 

    #mixin ratiion=w
    #dqdz=dqdwdwdz
    #dqdz=1/(1+w)^2dwdz

    #convert to specific humidity from mixing ratio
    q_ls_u           =     mr_ls_u/(1.0 + mr_u)**2 
    vqadv_u          =  v_adv_mr_u/(1.0 + mr_u)**2 

    
    dqdt_u           =    dmrdt_u/(1.0 + mr_u)**2 


    #divT    #"Horizontal Temp advection" units = "K/s";
    t_ls_u     =   units.Quantity(exp.divT[:,::-1,0,0],'kelvin/s')

    #T tendency "K/s"
    dTdt_u     =   units.Quantity(exp.dTdt[:,::-1,0,0],'kelvin/s')

    #vertdivT#"Vertical Temp advection""K/s"
    vTadv_u    =   units.Quantity(exp.vertdivT[:,::-1,0,0],'kelvin/s')
    
    #s       #"Dry satic energy/Cp" units = "K" ;
    s_ls_u     =   units.Quantity(exp.divs[:,::-1,0,0],'kelvin/s')

    #vertdivs#vertdivs"Vert. dry static energy adv./Cp"units = "K/s" ;
    vsadv_u       =   units.Quantity(exp.vertdivs[:,::-1,0,0],'kelvin/s')

    dsdt_u       =   units.Quantity(exp.dsdt[:,::-1,0,0],'kelvin/s')


    
    T_u          =   units.Quantity(exp.T[:,::-1,0,0],'kelvin')
    pressure_u   =   units.Quantity(exp.lev[::-1], 'Pa')

    theta_u       =   metpy.calc.potential_temperature(pressure_u, T_u)


    prec_srf_u    =   units.Quantity(exp.prec_srf[:,0,0],'mm/hour')

    p_sur_u   =   units.Quantity(exp.Ps[:,0,0], 'Pa')


    phi_u   =   units.Quantity(exp.phis[0,0],'m^2/s^2')


    u_u           =   units.Quantity(exp.u[:,::-1,0,0],'m/s')
    v_u           =   units.Quantity(exp.v[:,::-1,0,0],'m/s')

    u_surf           =   units.Quantity(exp.usrf[:,0,0],'m/s')
    v_surf           =   units.Quantity(exp.vsrf[:,0,0],'m/s')


    omega_u       =   units.Quantity(exp.omega[:,::-1,0,0],'Pa/s')

    w_u           =   metpy.calc.vertical_velocity(omega_u, pressure_u, T_u, mr_u)




    #unidades estao erradas par apoder fazer a transformacao é K/s
    divT    =   units.Quantity(exp.divT[:,::-1,0,0],'kelvin')

    h_advec_thil   =   metpy.calc.potential_temperature(pressure_u, divT).magnitude

    h_advec_thil_u =   units.Quantity(h_advec_thil,'kelvin/s') 


    #unidades estao erradas par apoder fazer a transformacao é K/s
    vdivT =   units.Quantity(exp.vertdivT[:,::-1,0,0],'kelvin')

    v_advec_thil   =   metpy.calc.potential_temperature(pressure_u,vdivT).magnitude

    v_advec_thil_u =   units.Quantity(v_advec_thil,'kelvin/s') 


    #
    Tg_u          =   units.Quantity(exp.Tsair[:,0,0],'kelvin')


    g       = 9.81 *units('m/s^2')   # [m/s^2]
    R_d     = 287.0*units('J/kg/K')   # [J/kg/K]
    c_p     = 1004.0*units('J/K/kg')

    L_v     = 2.5E6*units('J/kg')

    #Sensible upwardheat flux 
    SH          =   units.Quantity(exp.shflx[:,0,0],'W/m^2') 


    #convert to Km/s
    #SH1 = SH*ffc.R_dry*Tg_u/(ffc.c_p*p_sur_u) #convert to Km/s
    #print(SH1[0:50])
    
    SH_u = SH*R_d*Tg_u/(c_p*p_sur_u) #convert to Kelvin m/s

    #Latent upward heat flux 
    LH          =   units.Quantity(exp.lhflx[:,0,0],'W/m^2') 
    #convert to Km/s
    #LH          = LH*ffc.R_dry*Tg_u/(ffc.L_v*p_sur_u) #convert to Km/s
    LH_u          = LH*R_d*Tg_u/(L_v*p_sur_u) #convert to m/s

    #data listed in the case instructions
    z_sfc = 63.0     #m above sea level
    
    
    #level function of time
    z_u       = np.zeros((ndtp,ndlev))*units('m')
    z_u[0]    = 63*units('m')#alt    


    height = ffc.get_height_from_pres(T_u[0,:].magnitude,pressure_u.magnitude,z_sfc)


    Tv_u      =   metpy.calc.virtual_temperature(T_u, q_u).to('kg*K/kg')

    for k in range(0,1):
        for i in range(1,ndlev):

            z_u[k,i]   =  R_d*Tv_u[k,i]/g*np.log(pressure_u[i-1]/pressure_u[i])+z_u[k,i-1]    


    #ql and tke are not specified; set to zero
    ql          = np.zeros((ndlev,ndtp),dtype=float)
    qi          = np.zeros((ndlev,ndtp),dtype=float)
    ql          = np.zeros((ndlev,ndtp),dtype=float)
    qi          = np.zeros((ndlev,ndtp),dtype=float)
    tke         = np.zeros((ndlev,ndtp),dtype=float)
    rad_heating = np.zeros((ndlev,ndtp),dtype=float)
    u_g          = np.zeros((ndlev,ndtp),dtype=float)
    v_g          = np.zeros((ndlev,ndtp),dtype=float)


    T           =np.swapaxes(           T_u[::1,0:ndlev].magnitude,0,1) 
    theta       =np.swapaxes(       theta_u[::1,0:ndlev].magnitude,0,1)
    q           =np.swapaxes(           q_u[::1,0:ndlev].magnitude,0,1)
    t_ls        =np.swapaxes(        t_ls_u[::1,0:ndlev].magnitude,0,1)
    v_advec_thil=np.swapaxes(v_advec_thil_u[::1,0:ndlev].magnitude,0,1)
    dTdt        =np.swapaxes(        dTdt_u[::1,0:ndlev].magnitude,0,1)
    h_advec_thil=np.swapaxes(h_advec_thil_u[::1,0:ndlev].magnitude,0,1)
    q_ls        =np.swapaxes(        q_ls_u[::1,0:ndlev].magnitude,0,1)
    vqadv       =np.swapaxes(       vqadv_u[::1,0:ndlev].magnitude,0,1)
    dqdt        =np.swapaxes(        dqdt_u[::1,0:ndlev].magnitude,0,1)
    s_ls        =np.swapaxes(        s_ls_u[::1,0:ndlev].magnitude,0,1) 
    vsadv       =np.swapaxes(       vsadv_u[::1,0:ndlev].magnitude,0,1)
    dsdt        =np.swapaxes(        dsdt_u[::1,0:ndlev].magnitude,0,1) 
    u           =np.swapaxes(           u_u[::1,0:ndlev].magnitude,0,1)
    v           =np.swapaxes(           v_u[::1,0:ndlev].magnitude,0,1)
    w           =np.swapaxes(           w_u[::1,0:ndlev].magnitude,0,1)
    omega       =np.swapaxes(       omega_u[::1,0:ndlev].magnitude,0,1) 
    phi         =phi_u.magnitude
    pressure    =pressure_u[0:ndlev].magnitude
    Tg          =Tg_u.magnitude
    p_sur       =p_sur_u.magnitude
    z           =z_u[0,0:ndlev].magnitude
    SH          =SH_u.magnitude
    LH          =LH_u.magnitude

    #vTadv       =np.swapaxes(       vTadv_u.magnitude,0,1)
    #v_advec_thil=np.swapaxes(v_advec_thil_u.magnitude,0,1)

    return T ,theta,q,ql,qi,tke,\
           t_ls,v_advec_thil,dTdt,h_advec_thil,\
           q_ls,vqadv,dqdt, \
           s_ls,vsadv,dsdt,rad_heating, \
           u,v,u_g,v_g,w,omega, \
           phi,Tg,p_sur,z, \
           pressure,SH,LH

def variables_metpy_scm_ccpp(exp):


    #Lengh of the time array to search
    ndtp    = len(exp.tsec) 
    #Lengh of the time array to search
    ndlev   = len(exp.lev) 



    #Mixing ratio
    #divq    #"kg/kg/s""Horizontal q advection"
    mr_ls_u     =   units.Quantity(exp.divq[:,::-1,0,0],'kg/kg/s')

    #vertdivq#"Vertical q advection""kg/kg/s" ;
    v_adv_mr_u  =   units.Quantity(exp.vertdivq[:,::-1,0,0],'kg/kg/s')

    #T tendency "K/s"
    dmrdt_u     =   units.Quantity(exp.dqdt[:,::-1,0,0],'kg/kg/s')

    mr_u        =   units.Quantity(exp.q[:,::-1,0,0],'kg/kg')


    #convert to specific humidity from mixing ratio
    q_u           =   mr_u/(1.0 + mr_u) 

    #mixin ratiion=w
    #dqdz=dqdwdwdz
    #dqdz=1/(1+w)^2dwdz

    #convert to specific humidity from mixing ratio
    q_ls_u           =     mr_ls_u/(1.0 + mr_u)**2 
    vqadv_u          =  v_adv_mr_u/(1.0 + mr_u)**2 

    
    dqdt_u           =    dmrdt_u/(1.0 + mr_u)**2 


    #divT    #"Horizontal Temp advection" units = "K/s";
    t_ls_u     =   units.Quantity(exp.divT[:,::-1,0,0],'kelvin/s')

    #T tendency "K/s"
    dTdt_u     =   units.Quantity(exp.dTdt[:,::-1,0,0],'kelvin/s')

    #vertdivT#"Vertical Temp advection""K/s"
    vTadv_u    =   units.Quantity(exp.vertdivT[:,::-1,0,0],'kelvin/s')
    
    #s       #"Dry satic energy/Cp" units = "K" ;
    s_ls_u     =   units.Quantity(exp.divs[:,::-1,0,0],'kelvin/s')

    #vertdivs#vertdivs"Vert. dry static energy adv./Cp"units = "K/s" ;
    vsadv_u       =   units.Quantity(exp.vertdivs[:,::-1,0,0],'kelvin/s')

    dsdt_u       =   units.Quantity(exp.dsdt[:,::-1,0,0],'kelvin/s')


    
    T_u          =   units.Quantity(exp.T[:,::-1,0,0],'kelvin')
    pressure_u   =   units.Quantity(exp.lev[::-1], 'Pa')

    theta_u       =   metpy.calc.potential_temperature(pressure_u, T_u)


    prec_srf_u    =   units.Quantity(exp.prec_srf[:,0,0],'mm/hour')

    p_sur_u   =   units.Quantity(exp.Ps[:,0,0], 'Pa')


    phi_u   =   units.Quantity(exp.phis[0,0],'m^2/s^2')


    u_u           =   units.Quantity(exp.u[:,::-1,0,0],'m/s')
    v_u           =   units.Quantity(exp.v[:,::-1,0,0],'m/s')

    u_surf           =   units.Quantity(exp.usrf[:,0,0],'m/s')
    v_surf           =   units.Quantity(exp.vsrf[:,0,0],'m/s')


    omega_u       =   units.Quantity(exp.omega[:,::-1,0,0],'Pa/s')

    w_u           =   metpy.calc.vertical_velocity(omega_u, pressure_u, T_u, mr_u)




    #unidades estao erradas par apoder fazer a transformacao é K/s
    divT    =   units.Quantity(exp.divT[:,::-1,0,0],'kelvin')

    h_advec_thil   =   metpy.calc.potential_temperature(pressure_u, divT).magnitude

    h_advec_thil_u =   units.Quantity(h_advec_thil,'kelvin/s') 


    #unidades estao erradas par apoder fazer a transformacao é K/s
    vdivT =   units.Quantity(exp.vertdivT[:,::-1,0,0],'kelvin')

    v_advec_thil   =   metpy.calc.potential_temperature(pressure_u,vdivT).magnitude

    v_advec_thil_u =   units.Quantity(v_advec_thil,'kelvin/s') 


    #
    Tg_u          =   units.Quantity(exp.Tsair[:,0,0],'kelvin')


    g       = 9.81 *units('m/s^2')   # [m/s^2]
    R_d     = 287.0*units('J/kg/K')   # [J/kg/K]
    c_p     = 1004.0*units('J/K/kg')

    L_v     = 2.5E6*units('J/kg')

    #Sensible upwardheat flux 
    SH          =   units.Quantity(exp.shflx[:,0,0],'W/m^2') 


    #convert to Km/s
    #SH1 = SH*ffc.R_dry*Tg_u/(ffc.c_p*p_sur_u) #convert to Km/s
    #print(SH1[0:50])
    
    SH_u = SH*R_d*Tg_u/(c_p*p_sur_u) #convert to Kelvin m/s

    #Latent upward heat flux 
    LH          =   units.Quantity(exp.lhflx[:,0,0],'W/m^2') 
    #convert to Km/s
    #LH          = LH*ffc.R_dry*Tg_u/(ffc.L_v*p_sur_u) #convert to Km/s
    LH_u          = LH*R_d*Tg_u/(L_v*p_sur_u) #convert to m/s

    #data listed in the case instructions
    z_sfc = 63.0     #m above sea level
    
    
    #level function of time
    z_u       = np.zeros((ndtp,ndlev))*units('m')
    z_u[0]    = 63*units('m')#alt    


    height = ffc.get_height_from_pres(T_u[0,:].magnitude,pressure_u.magnitude,z_sfc)


    Tv_u      =   metpy.calc.virtual_temperature(T_u, q_u).to('kg*K/kg')

    for k in range(0,1):
        for i in range(1,ndlev):

            z_u[k,i]   =  R_d*Tv_u[k,i]/g*np.log(pressure_u[i-1]/pressure_u[i])+z_u[k,i-1]    


    #ql and tke are not specified; set to zero
    ql          = np.zeros((ndlev,ndtp),dtype=float)
    qi          = np.zeros((ndlev,ndtp),dtype=float)
    ql          = np.zeros((ndlev,ndtp),dtype=float)
    qi          = np.zeros((ndlev,ndtp),dtype=float)
    tke         = np.zeros((ndlev,ndtp),dtype=float)
    rad_heating = np.zeros((ndlev,ndtp),dtype=float)
    u_g          = np.zeros((ndlev,ndtp),dtype=float)
    v_g          = np.zeros((ndlev,ndtp),dtype=float)


    T           =np.swapaxes(           T_u[::1,0:ndlev].magnitude,0,1) 
    theta       =np.swapaxes(       theta_u[::1,0:ndlev].magnitude,0,1)
    q           =np.swapaxes(           q_u[::1,0:ndlev].magnitude,0,1)
    t_ls        =np.swapaxes(        t_ls_u[::1,0:ndlev].magnitude,0,1)
    v_advec_thil=np.swapaxes(v_advec_thil_u[::1,0:ndlev].magnitude,0,1)
    dTdt        =np.swapaxes(        dTdt_u[::1,0:ndlev].magnitude,0,1)
    h_advec_thil=np.swapaxes(h_advec_thil_u[::1,0:ndlev].magnitude,0,1)
    q_ls        =np.swapaxes(        q_ls_u[::1,0:ndlev].magnitude,0,1)
    vqadv       =np.swapaxes(       vqadv_u[::1,0:ndlev].magnitude,0,1)
    dqdt        =np.swapaxes(        dqdt_u[::1,0:ndlev].magnitude,0,1)
    s_ls        =np.swapaxes(        s_ls_u[::1,0:ndlev].magnitude,0,1) 
    vsadv       =np.swapaxes(       vsadv_u[::1,0:ndlev].magnitude,0,1)
    dsdt        =np.swapaxes(        dsdt_u[::1,0:ndlev].magnitude,0,1) 
    u           =np.swapaxes(           u_u[::1,0:ndlev].magnitude,0,1)
    v           =np.swapaxes(           v_u[::1,0:ndlev].magnitude,0,1)
    w           =np.swapaxes(           w_u[::1,0:ndlev].magnitude,0,1)
    omega       =np.swapaxes(       omega_u[::1,0:ndlev].magnitude,0,1) 
    phi         =phi_u.magnitude
    pressure    =pressure_u[0:ndlev].magnitude
    Tg          =Tg_u.magnitude
    p_sur       =p_sur_u.magnitude
    z           =z_u[0,0:ndlev].magnitude
    SH          =SH_u.magnitude
    LH          =LH_u.magnitude

    #vTadv       =np.swapaxes(       vTadv_u.magnitude,0,1)
    #v_advec_thil=np.swapaxes(v_advec_thil_u.magnitude,0,1)

    return T ,theta,q,ql,qi,tke,\
           t_ls,v_advec_thil,dTdt,h_advec_thil,\
           q_ls,vqadv,dqdt, \
           s_ls,vsadv,dsdt,rad_heating, \
           u,v,u_g,v_g,w,omega, \
           phi,Tg,p_sur,z, \
           pressure,SH,LH

