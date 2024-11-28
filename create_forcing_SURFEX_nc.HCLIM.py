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
    namelist = read_namelist('../nam/namelist_create_forcing_SURFEX_nc.HCLIM.yaml')

    # 2 - Read variables info
    var = read_dataframe('../nam/namelist_variables_SURFEX_forcing.HCLIM')

    for pgw_case in namelist['PGW']:


        print('3 - Read grid info')
        namelist = rewrite_output_file(namelist, pgw_case)

        # 4 - Read grid info
        print('4 - Read grid info')
        time, longitude, latitude = read_grid_info_time(namelist)

        # 5 - Save grid information
        print('5 - Save grid information')
        save_grid_info(namelist, pgw_case, time, longitude, latitude)

        # 6 - Include time step
        print('6 - Include time step')
        include_timestep(namelist)

        # 7 - Loop over all variables
        #for varname in tqdm(var.SFXname.values, desc='Variables:', leave=True):
        for varname in var.SFXname.values:
	
            include_var(namelist, var[var.SFXname==varname], pgw_case)

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

    return np.mod(np.rad2deg(np.arctan2(-u, -v)) + 180, 360)
    #return np.mod(np.rad2deg(np.arctan2(v, u)) + 180, 360)
    #return np.rad2deg(np.arctan2(u,v))

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

def read_dataframe(filename):

    """
    Read dataframe

    Parameters:
    filename (str) : path to the namelist file

    Returns:
    dataframe : including all the information regarding the variables

    """

    return pd.read_csv(filename, sep=',',skipinitialspace = True)

def rewrite_output_file(namelist, pgw_case):

    """
    Read yaml namelist file

    Parameters:
    namelist (dict) : information regarding the domain and the year

    Returns:
    dictionary: updated version of the namelist

    """

    namelist['forcingfile'] = 'FORCING_HCLIM_'+namelist['domain']+'_'+str(namelist['year'])+'_'+pgw_case+'.nc'

    return namelist

def ref_height():

    f  = Dataset('/nobackup/rossby27/users/sm_aital/analysis/surfex_online_vs_offline/forcing_files/FORCING_L65_2018.nc','r')
    t65 = f.variables['Tair'][:,:]
    f.close()

    R = 8.3144598
    g = 9.81
    M =  0.0289644
    L = 0.0065	

    C1 = 0.5*(0.997039230000000+1)
    C2 = M*g/(R*L)

    C3 = np.exp(np.log(C1)/C2)

    z65 = t65*((1/(C3)-1))/L

    return np.mean(z65,0)

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

    print(num_points)
    # Dimensions
    time = new_file.createDimension('time', None)
    numpoints = new_file.createDimension('Number_of_points', num_points)
   
    namelist['numpoints'] = num_points

    # time
    print('Define time')
    time = new_file.createVariable('time', 'f4', ('time'))

    time[:]=time0-time0[0]

    time.standard_name = "time" ;
    time.units = "hours since "+str(namelist['year'])+"-05-25 00:00:00" ;
    #time.units = "hours since "+str(namelist['year'])+"-07-02 00:00:00" ;
    time.calendar = "standard" ;
    time.axis = "T" ;

    namelist['numtime'] = len(time)
    
    # Latitude / Longitude
    print('Define lat/lon')
    latitude  = new_file.createVariable('LAT','f4',('Number_of_points'))
    longitude = new_file.createVariable('LON','f4',('Number_of_points'))

    latitude[:] = np.reshape(latitude0, (num_points))
    latitude.long_name = "Latitude" ;

    longitude[:] = np.reshape(longitude0, (num_points))
    longitude.long_name = "Longitude" ;

    # Reference heights
    print('Define reference height')
    zref  = new_file.createVariable('ZREF','f4',('Number_of_points'))
    uref  = new_file.createVariable('UREF','f4',('Number_of_points'))

    #zref[:] = ref_height()
    #zref[:] = 2.0
    #zref[:] = 12.5
    zref[:] = 50.0

    zref.long_name = "Reference_Height" ;
    zref.units = "m"

    #uref[:] = ref_height()
    #uref[:] = 10.0
    #uref[:] = 12.5
    uref[:] = 50.0

    uref.long_name = "Reference_Height_for_Wind" ;
    uref.units = "m"

    # CO2   
    print('Define CO2')
    co2  = new_file.createVariable('CO2air','f4',('time','Number_of_points'))
    co2[:,:] = ppm_2_con(namelist['CO2'][pgw_case])

    co2.units = "kg/m3"
    co2.long_name = "Near_Surface_CO2_Concentration"

    print('Close file')
    new_file.close()

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

    nc = Dataset(namelist['ini_file'],'r')

    # Longitude, latitude and time variables
    longitude0 = nc.variables['lon'][:-11,:-11]
    latitude0  = nc.variables['lat'][:-11,:-11]
    time0      = nc.variables['time'][:]

    nc.close()

    return time0, longitude0, latitude0

def include_timestep(namelist):

    """
    Include timestep information into forcing file

    Parameters:
    namelist (dict) : dictionary containing all the variables names

    """

    # Time step
    timestep = 3600.0
 
    cmd1 = "ncap2 -O -s 'FRC_TIME_STP=double("+str(timestep)+")' "+namelist['out_dir']+"tmp_forcing.nc "+namelist['out_dir']+namelist['forcingfile']
    cmd2 = "ncatted -a longname,FRC_TIME_STP,o,c,'Forcing_Time_Step' "+namelist['out_dir']+namelist['forcingfile']
    cmd3 = "rm -rf "+namelist['out_dir']+"tmp_forcing.nc "

    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

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

    try:
        var1_nc = ncfile.variables[var1][:,:-11,:-11].data
    except:
        var1_nc = ncfile.variables[var1][:-11,:-11].data

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
       return np.reshape(var[:,:], (namelist['numpoints']))

    elif dim == 2:
       return np.reshape(var, (var.shape[0],namelist['numpoints']))
    
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

    filename = namelist['ini_file']

    var1_nc, var2_nc = extract_rawdata(filename, var1[0], var2[0])
 
    var1_nc = reshape_var(namelist, var1_nc, dim)
    if var2_nc is not None: var2_nc = reshape_var(namelist, var2_nc, dim)

    if sub[0]: var1_nc = var1_nc - var2_nc

    var_out = apply_fun(fun[0],var1_nc, var2_nc, sub[0])

    return var_out

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

    var_values[np.where(var_values<0)] = 0.
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

    var_values[np.where(var_values<0)] = 0.
    nc_var[:,:] = var_values
    nc_var.long_name = var.long_name.values[0].strip()
    nc_var.units = var.units.values[0].strip()

    new_file.close()
# ------------------------------------------------------------------------------

create_forcing_SURFEX_nc()
