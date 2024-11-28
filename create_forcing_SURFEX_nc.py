import numpy as np
import os
import yaml
import pandas as pd
import math

from netCDF4 import Dataset

from tqdm import tqdm
from yaml.loader import SafeLoader

# ------------------------------------------------------------------------------

def create_forcing_SURFEX_nc():

    # 1 - Read namelist
    namelist = read_namelist('../nam/namelist_create_forcing_SURFEX_nc.yaml')

    # 2 - Read variables info
    var = read_dataframe('../nam/namelist_variables_SURFEX_forcing')

    for pgw_case in namelist['PGW']:

        # 3 - Read grid info
        namelist = rewrite_output_file(namelist, pgw_case)

        # 4 - Read grid info
        time, longitude, latitude = read_grid_info_time(namelist)

        # 5 - Save grid information
        save_grid_info(namelist, pgw_case, time, longitude, latitude)

        # 6 - Include time step
        include_timestep(namelist)

        # 7 - Loop over all variables
        for varname in tqdm(var.SFXname.values, desc='Variables:', leave=True):
	
            include_var(namelist, var[var.SFXname==varname], pgw_case)

def d2m_2_q2m(dew,p):

    """ 
    Transform dew temperature to specific humidity

    Parameters:
    dew(float) : dew temperature (K)
    p (float)  : surface pressure (Pa)

    Returns:
    float: specific humidity (kg/kg)

    """

    e0 = 611.3 # saturation vapor pressure in Pa
    # e0 and Pressure have to be in same units
    c_water = 5423 # L/R for water in Kelvin
    T0 = 273.15 # Kelvin

    #calculating specific humidity, q directly from dew point temperature
    #using equation 4.24, Pg 96 Practical Meteorolgy (Roland Stull)
    q = (622 * e0 * np.exp(c_water * (dew - T0)/(dew * T0)))/p # g/kg 
    # 622 is the ratio of Rd/Rv in g/kg

    # return in kg/kg
    return q/1000.


def m2s2_2_m(z):

    """ 
    Transform geopotential height to altitude

    Parameters:
    z(float) : geopotential height (m2/s2)

    Returns:
    float: orography height (m)

    """

    g0 = 9.81 # m/s2

    return z/g0


def m_2_kgm2s(pr, tt=3.):

    """ 
    Transform precipitation from meters to flux

    Parameters:
    pr(float)             : precipitation (m)
    tt (float, optional)  : time step of the interval (hours)

    Returns:
    float: precipitation flux (kg/m2/s)

    """

    rho0 = 1000 # kg/m3 

    tstep = tt*3600 # tt in hours

    return pr*rho0/tstep


def Jm2_2_Wm2(f,tt=3.):

    """ 
    Transform radiation from energy to flux

    Parameters:
    f(float)              : radiation (J)
    tt (float, optional)  : time step of the interval (hours)

    Returns:
    float: radiation flux (W/m2)

    """

    tstep = tt*3600 # tt in hours

    return f/tstep


def windspeed(u,v):
 
    """ 
    Compute windspeed from u and v

    Parameters:
    u(float)  : zonal wind (m/s)
    v(float)  : meridional wind (m/s)

    Returns:
    float: wind speed (m/s)

    """

    field = np.sqrt(np.square(u) + np.square(v))
    np.where(field < 0.005, field, 0)

    return field


def winddir(u,v):

    """ 
    Compute windspeed from u and v

    Parameters:
    u(float)  : zonal wind (m/s)
    v(float)  : meridional wind (m/s)

    Returns:
    float: wind direction

    """

    return np.mod(np.rad2deg(np.arctan2(v, u)) + 180, 360)



def ppm_2_con(co2value):

    """
    Tranform ppmv to kg/m3

    Parameters:
    co2value (float) : co2 concentration in ppmv

    Returns:
    float : co2 concentration in kg/m3

    """

    rho_co2 = 1.87 # kg/m3

    return rho_co2*co2value*1e-6


def read_dataframe(filename):

    """
    Read dataframe

    Parameters:
    filename (str) : path to the namelist file

    Returns:
    dataframe : including all the information regarding the variables

    """

    return pd.read_csv(filename, sep=',',skipinitialspace = True)


def read_namelist(filename):

    """
    Read yaml namelist file

    Parameters:
    filename (str) : path to the namelist file

    Returns:
    dictionary : dictionary with all the namelist information

    """

    with open(filename) as f:
        namelist = yaml.load(f, Loader=SafeLoader)

    return namelist


def rewrite_output_file(namelist, pgw_case):

    """
    Read yaml namelist file

    Parameters:
    namelist (dict) : information regarding the domain and the year

    Returns:
    dictionary: updated version of the namelist

    """

    namelist['forcingfile'] = 'FORCING_ERA5_'+namelist['domain']+'_'+str(namelist['year'])+'_'+pgw_case+'.nc'

    return namelist

def read_grid_info_time(namelist):

    """ 
    Read time and lat/lon info from file

    Parameters:
    namelist (dict) : dictionary containing all the variables names

    Returns:
    time0 (float): time variables
    longitude (float): longitude information
    latitude (float): latitude information

    """

    nc = Dataset(namelist['nonaccu_file'],'r')

    # Longitude, latitude and time variables
    longitude0 = nc.variables['lon'][:-11,:-11]
    latitude0  = nc.variables['lat'][:-11,:-11]
    time0      = nc.variables['time'][:]

    nc.close()

    return time0, longitude0, latitude0


def save_grid_info(namelist, pgw_case, time0, longitude0, latitude0):

    """
    Create file with grid information and other variables

    Parameters:
    namelist (dict)    : dictionary containing all the variables names
    pgw_case (str)     : PGW case
    time0 (float)      : time array
    longitude0 (float) : longitude array
    latitude0 (float)  : latitude array

    """

    new_file = Dataset(namelist['out_dir']+'tmp_forcing.nc','w')

    # Number of points
    num_points = longitude0.shape[0]*longitude0.shape[1]

    # Dimensions
    time = new_file.createDimension('time', None)
    numpoints = new_file.createDimension('Number_of_points', num_points)
   
    namelist['numpoints'] = num_points

    # time 
    time = new_file.createVariable('time', 'f4', ('time'))

    time[:]=time0-time0[0]

    time.standard_name = "time" ;
    time.units = "hours since "+str(namelist['year'])+"-1-1 00:00:00" ;
    time.calendar = "standard" ;
    time.axis = "T" ;

    namelist['numtime'] = len(time)
    
    # Latitude / Longitude
    latitude  = new_file.createVariable('LAT','f4',('Number_of_points'))
    longitude = new_file.createVariable('LON','f4',('Number_of_points'))

    latitude[:] = np.reshape(latitude0, (num_points))
    latitude.long_name = "Latitude" ;

    longitude[:] = np.reshape(longitude0, (num_points))
    longitude.long_name = "Longitude" ;

    # Reference heights
    zref  = new_file.createVariable('ZREF','f4',('Number_of_points'))
    uref  = new_file.createVariable('UREF','f4',('Number_of_points'))

    zref[:] = 2.0
    zref.long_name = "Reference_Height" ;
    zref.units = "m"

    uref[:] = 10.0
    uref.long_name = "Reference_Height_for_Wind" ;
    uref.units = "m"

    # CO2   
    co2  = new_file.createVariable('CO2air','f4',('time','Number_of_points'))
    co2[:,:] = ppm_2_con(namelist['CO2'][pgw_case])

    co2.units = "kg/m3"
    co2.long_name = "Near_Surface_CO2_Concentration"

    new_file.close()


def include_timestep(namelist):

    """
    Include timestep information into forcing file

    Parameters:
    namelist (dict) : dictionary containing all the variables names

    """

    # Time step
    timestep = 3*3600.0
 
    cmd1 = "ncap2 -O -s 'FRC_TIME_STP=double("+str(timestep)+")' "+namelist['out_dir']+"tmp_forcing.nc "+namelist['out_dir']+namelist['forcingfile']
    cmd2 = "ncatted -a longname,FRC_TIME_STP,o,c,'Forcing_Time_Step' "+namelist['out_dir']+namelist['forcingfile']
    cmd3 = "rm -rf "+namelist['out_dir']+"tmp_forcing.nc "

    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)


def extract_rawdata(filename, var1, var2):

    """
    Extract data from input netcdf file

    Parameters:
    filename (str) : path to file
    var1 (str)      : name of the main variable
    var2 (str)      : name of secondary variable if needed

    Returns
    float : 1/2-D fields

    """

    ncfile = Dataset(filename,'r')

    var1_nc = ncfile.variables[var1][:,:-11,:-11].data

    if var2 is not np.nan:  
       var2_nc = ncfile.variables[var2][:,:-11,:-11].data
    else:
       var2_nc = None

    ncfile.close()
   
    return var1_nc, var2_nc


def reshape_var(namelist,var,dim):

    """
    Reduce dimensions by reshaping

    Parameters:
    namelist (dict) : dictionary containing all the variables names
    var (float)     : variable to be reduced
    dim (integer)   : dimension of the out field

    Returns
    float : 1/2-D field

    """

    if dim == 1:
       return np.reshape(var[0,:,:], (namelist['numpoints']))

    elif dim == 2:
       return np.reshape(var, (var.shape[0],namelist['numpoints']))


def accumulate_fields(namelist, var):

    """
    Reduce dimensions by reshaping

    Parameters:
    namelist (dict) : dictionary containing all the variables names
    var (float)     : variable to be reduced
    
    Returns
    float : 1/2-D field

    """
    
    accu_fld = np.zeros((namelist['numtime'],namelist['numpoints']))
    for ii in range(namelist['numpoints']):
        accu_fld[:,ii] = np.sum(var[:,ii].reshape(-1,3),axis=1)

    return  accu_fld

def apply_fun(fun, var1, var2, sub = False):

    """
    Apply function to variable

    Parameters:
    fun (str)       : name of function
    var1 (float)    : main variable
    var2 (float)    : secondary variable
    sub(logical)    : is it a variable that has been substracted?
    
    Returns
    float : 1/2-D field

    """    

    if fun is np.nan:
        return var1
    elif var2 is None or sub:
        return eval(fun+"(var1)")
    else:
        return eval(fun+"(var1, var2)")


def extract_data(namelist, var1, var2, fun, accu, sub, dim):

    """
    Extract data from input netcdf file

    Parameters:
    namelist (dict) : dictionary containing all the variables names
    var1 (str)      : name of the main variable
    var2 (str)      : name of secondary variable if needed
    fun (str)       : function to be applied
    accu (logical)  : accumulated field
    sub (logical)   : subtract between var1 and var2
    dim (integer)   : dimension of the out field

    Returns
    float : 1/2-D field

    """

    filename = namelist['nonaccu_file']
    if accu[0]: filename = namelist['accu_file']

    var1_nc, var2_nc = extract_rawdata(filename, var1[0], var2[0])
 
    var1_nc = reshape_var(namelist, var1_nc, dim)
    if var2_nc is not None: var2_nc = reshape_var(namelist, var2_nc, dim)

    if sub[0]: var1_nc = var1_nc - var2_nc

    if accu[0]: var1_nc = accumulate_fields(namelist,var1_nc)

    var_out = apply_fun(fun[0],var1_nc, var2_nc, sub[0])

    return var_out


def include_var_1D(namelist, var):

    """
    Include 1D variables into forcing file

    Parameters:
    namelist (dict) : dictionary containing all the information about the input files
    var(dict)       : dictionary containing all the information regarding the variables

    """

    var_values = extract_data(namelist, var.var1.values, var.var2.values, var.function.values, var.accu.values, var.subst.values, 1)

    new_file = Dataset(namelist['out_dir']+namelist['forcingfile'],'r+')

    nc_var = new_file.createVariable(var.SFXname.values[0],'f4',('Number_of_points'))

    nc_var[:] = var_values
    nc_var.long_name = var.long_name.values[0].strip()
    nc_var.units = var.units.values[0].strip()

    new_file.close()

def include_var_2D(namelist, var, pgw_case):

    """
    Include 2D variables into forcing file

    Parameters:
    namelist (dict) : dictionary containing all the information about the input files
    var(dict)       : dictionary containing all the information regarding the variables

    """
    
    var_values = extract_data(namelist, var.var1.values, var.var2.values, var.function.values, var.accu.values, var.subst.values, 2)

    new_file = Dataset(namelist['out_dir']+namelist['forcingfile'],'r+')

    nc_var = new_file.createVariable(var.SFXname.values[0],'f4',('time','Number_of_points'))

    nc_var[:,:] = pgw_correction(var_values, pgw_case, namelist['PGW_path'], var.delta_var.values[0], var.delta_method.values[0], namelist['leapyear'])
    nc_var.long_name = var.long_name.values[0].strip()
    nc_var.units = var.units.values[0].strip()

    new_file.close()


def pgw_correction(var_value, pgw_case, pgw_path, delta_var, delta_method, leapyear):

    """
    Modify variable value incorporating a PGW delta

    Parameters:
    var_values (float) : variable of interest
    pgw_case (str)     : pgw case    
    pgw_path (str)     : path to the deltas
    delta_var (str)    : name of the delta variable
    delta_method (int) : method type
    leapyear(logical)  : is it a leap year?

    Returns:
    float: variable corrected

    """

    if pgw_case == 'PGW1' or delta_method == 0:
        # No delta
        return var_value
    elif delta_method == 1:
        # Add delta
        return add_delta(var_value, pgw_case, pgw_path, delta_var, leapyear)
    elif delta_method == 2:
        # Delta multiplication
        return multiply_delta(var_value, pgw_case, pgw_path, delta_var, leapyear)
    elif delta_method == 3:
        # Combined method to apply delta
        return combined_delta(var_value, pgw_case, pgw_path, delta_var, leapyear)


def delta_leapday(delta,nt):

    """
    Adjust delta to the leap day, repeat delta 

    Parameters:
    delta (float)    : variable of interest
    nt (int)         : time steps

    Returns:
    float: delta value adjusted
    int  : number of time steps adjusted

    """

    # 3h time step
    if nt == 2921:
         nshift = 8
         ntstep = 472

    # Daily data
    elif nt == 366:
         nshift = 1
         ntstep = 58

    else:
         return delta, nt
    
    delta_tmp = np.zeros((nt+nshift,delta.shape[1])) 

    delta_tmp[:ntstep,:] = delta[:ntstep,:]
    delta_tmp[ntstep:ntstep+nshift] = delta[ntstep-nshift:ntstep,:]
    delta_tmp[ntstep+nshift:,:] = delta[ntstep:,:]

    return delta_tmp, nt+nshift


def monthly_datastep(leapyear):

    """
    Return the time step in each month

    Parameters:
    leapyear (logical) : is it a leap year?

    Returns:
    list: time step corresponding to the last day of each month
    list: number of days per month

    """

    if leapyear:
         month_tstep = [0,744, 1444, 2184, 2904, 3648, 4368, 5112, 5856, 6556, 7320, 8040, 8784]
         month_days  = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
         month_tstep = [0,744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6532, 7296, 8016, 8760]
         month_days  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    return month_tstep, month_days


def read_delta(deltafile, delta_var, leapyear):

    """
    Read PGW delta from file

    Parameters:
    deltafile (str)    : variable of interest
    delta_var (str)    : name of the delta variable
    leapyear (logical) : is it a leap year?

    Returns:
    float : value of delta
    int   : number of time steps 

    """

    ncdelta = Dataset(deltafile,'r')
    delta = ncdelta.variables[delta_var][:,:-11,:-11]
    ncdelta.close()

    nt, nx, ny = delta.shape

    delta = np.reshape(delta, (nt,nx*ny))

    if leapyear: delta, nt = delta_leapday(delta, nt)

    return delta, nt

    
def add_delta(var_value, pgw_case, pgw_path, delta_var, leapyear):

    """
    Add PGW delta

    Parameters:
    var_values (float) : variable of interest
    pgw_case (int)     : pgw case    
    pgw_path (str)     : path to the deltas
    delta_var (str)    : name of the delta variable

    Returns:
    float : corrected variable value

    """

    deltafile = pgw_path+pgw_case+'/ensemble/'+delta_var+'_'+pgw_case+'.nc'

    delta, _ = read_delta(deltafile, delta_var, leapyear)

    var_value = var_value + delta

    var_value[var_value<0.] = 0.
    
    return var_value 


def multiply_delta(var_value, pgw_case, pgw_path, delta_var, leapyear):

    """
    Multiply PGW delta

    Parameters:
    var_values (float) : variable of interest
    pgw_case (int)     : pgw case    
    pgw_path (str)     : path to the deltas
    delta_var (str)    : name of the delta variable

    Returns:
    float : corrected variable value

    """

    deltafile = pgw_path+pgw_case+'/ensemble/'+delta_var+'_'+pgw_case+'.nc'

    delta, nt = read_delta(deltafile, delta_var, leapyear)

    for tt in range(nt):  
        var_value[8*tt:8*(tt+1),:]  = var_value[8*tt:8*(tt+1),:]*delta[tt,:]

    var_value[-1,:] = var_value[-1,:]*delta[0,:]

    var_value[var_value<0.] = 0.

    return var_value

def combined_delta(var_value, pgw_case, pgw_path, delta_var, leapyear):

    """
    Combined PGW delta

    Parameters:
    var_values (float) : variable of interest
    pgw_case (int)     : pgw case    
    pgw_path (str)     : path to the deltas
    delta_var (str)    : name of the delta variable

    Returns:
    float : corrected variable value

    """

    deltafile = pgw_path+pgw_case+'/ensemble/'+delta_var+'_'+pgw_case+'.nc'

    delta, nt = read_delta(deltafile, delta_var, leapyear)

    month_tstep, month_days = monthly_datastep(leapyear)

    for tt in range(nt):

        accumulate_val = np.sum(var_value[month_tstep[tt]:month_tstep[tt+1],:],0)/(8*month_days[tt])
        accumulate_val[accumulate_val<0.] = 0.

        new_accumulate_val = accumulate_val + delta[tt,:]
        new_accumulate_val[new_accumulate_val<0.] = 0.

        new_delta = new_accumulate_val/accumulate_val
        new_delta[accumulate_val==0.] = 0.
        
        var_value[month_tstep[tt]:month_tstep[tt+1],:] = var_value[month_tstep[tt]:month_tstep[tt+1],:]*new_delta

    var_value[-1,:] = var_value[-1,:]*new_delta

    var_value[var_value<0.] = 0.

    return var_value


def include_var(namelist, var, pgw_case):

    """
    Includes variables into forcing file

    Parameters:
    namelist (dict) : dictionary containing all the information about the input files
    var(dict)       : dictionary containing all the information regarding the variables

    """

    if var.dim.values == 1:
       include_var_1D(namelist, var)
    elif var.dim.values == 2:
       include_var_2D(namelist, var, pgw_case)


# ------------------------------------------------------------------------------

create_forcing_SURFEX_nc()

# ----------------------------------------------------------------------------------
