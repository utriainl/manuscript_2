import numpy as np
import os
import yaml
import pandas as pd
import math
import sys

from netCDF4 import Dataset
from tqdm import tqdm
from yaml.loader import SafeLoader

# ------------------------------------------------------------------------------



def is_year_leap(year):

    """
    Check if a specific year is leap

    Parameters:
    year (int) : year of interest

    Returns:
    Logical : True (False) for leap (non-leap year)

    """

    if (year % 400 == 0) and (year % 100 == 0):
       return True
    elif (year % 4 ==0) and (year % 100 != 0):
       return True
    else:
       return False


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


def read_file_pandas(filename, col_list, names_list):

    """
    Read data from file and transferred it into a dataframe

    Parameters:
    filename (str)          : path to the file of interest
    col_list (list of int)  : columns to be extracted
    names_list(list of str) : names for the columns

    Returns:
    dataframe : extracted information

    """

    df = pd.read_csv(filename, header=None, skipinitialspace = True,
            usecols=col_list, names=names_list, sep=" ")

    return df


def read_SWL_info(namelist,caselist):

    """
    Read specific warming level information from files

    Parameters:
    namelist (dict): dictionary containing paths and names
    caselist (str) : SWL case

    Returns:
    dictionary : dictionary containing SWL information

    """

    dir_file = namelist['SWL_info_dir']

    SWL_info = {}
    for case in caselist:

        name_file = namelist[case]['info_file']

        df = read_file_pandas(dir_file+name_file, [2], ['year'])

        SWL_info[case] = df.year.values
    
    return SWL_info 


def create_workdir(namelist, member=None, cases=None):

    """
    Create work directories 

    Parameters:
    namelist (dict): dictionary containing paths and names
    member (int)   : member number
    cases (str)    : PGW information

    """

    if member == None:
        os.system('mkdir -p '+namelist['output_dir']+'tmpDIR')

    else:

        for case in cases:
            os.system('mkdir -p '+namelist['output_dir']+'HCLIM_grib/'+case+'/m'+str(member+101))


def remove_workdir(namelist):

    """
    Remove working directory

    Parameters:
    namelist (dict): dictionary containing path information

    """

    os.system('rm -rf '+namelist['output_dir']+'tmpDIR')


def clean_workdir(namelist,member):

    """
    Clean working directory by removing all the temporary netCDF files

    Parameters:
    namelist (dict): dictionary containing path information
    
    """


    os.system('rm -rf '+namelist['output_dir']+'tmpDIR/*m'+str(member+101)+'*')


def extract_LESM_files(namelist, SWL_info, member):

    """
    Extract LESM files by remapping to a new domain and save them into the member file. 
    Remove the 29th February from leap years.

    Parameters:
    namelist (dict): dictionary containing paths, names and variables names
    SWL_info (dict): SWL information
    member (int)   : member information

    """

    main_path = namelist['main_dir']
    yearshift = namelist['PGW_years']

    # Join all variables except last one
    varlist = namelist['varnames_sub'][:]
    varlist.remove('var134')
    varlist.remove('var157')
    
    listvars  = ','.join(varlist)

    for case in ['SWL09']: #SWL_info.keys():

        yearref=SWL_info[case][member]
        yearlist=np.arange(yearref-yearshift+1,yearref+1)

        str1 = main_path.replace("CASECASE", str(namelist[case]['case']))
        str2 = str1.replace('EXPEXP', namelist[case]['exp'])

        path_files = str2.replace('MMM',str(101+member))

        for year in yearlist:
            for month in range(1,13):

                os.system('cdo -s -daymean -remapbil,'+namelist['grid_info']+' -selvar,'+listvars+' -mergetime '+path_files+str(year)+str(month).zfill(2)+'* '
                           +namelist['output_dir']+'tmpDIR/m'+str(101+member)+str(year)+str(month).zfill(2)+case+'.grb')

        os.system('cdo -s mergetime -del29feb '+namelist['output_dir']+'tmpDIR/m'+str(101+member)+'*'+case+'* '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member))


def delta_substract_LESM(namelist, SWL_info, member):

      """
      Compute PGW deltas by substracting two SWL

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      member (int)   : member number

      """

      case0 = list(SWL_info.keys())[0]

      for case in list(list(SWL_info.keys())[1:]):

           for varname in namelist['varnames_sub']:
      
                if varname == 'var152': continue
                if varname == 'var133': continue

                os.system('cdo -s sub '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_'+varname+'.grb '
                          +namelist['output_dir']+'tmpDIR/'+case0+'_m'+str(101+member)+'_'+varname+'.grb '
                          +namelist['output_dir']+'tmpDIR/'+case+'-'+case0+'_m'+str(101+member)+'_'+varname+'.grb ')


def delta_archive(namelist, SWL_info, member):

      """
      Archive PGW delta information to a specific member folder

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      member (int)   : member number

      """

      case0 = list(SWL_info.keys())[0]

      for case in list(list(SWL_info.keys())[1:]):
            for varname in namelist['varnames_sub']:

                if varname == 'var152': continue
                if varname == 'var133': continue

                os.system('cdo -s -O ydaymean -del29feb '+namelist['output_dir']+'tmpDIR/'+case+'-'+case0+'_m'+str(101+member)+'_'+varname+'.grb '
                          +namelist['output_dir']+'HCLIM_grib/PGW'+case[3]+'/m'+str(member+101)+'/'+varname+'_'+str(101+member)
                          +'_'+case+'-'+case0+'.grb')
               

def separate_LESM_files_variables(namelist, SWL_info, member):

      """
      Separate grib files into a file per variable

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      member (int)   : member number

      """

      for case in SWL_info.keys(): 
    
         os.system('cdo -s splitname '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+' '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_')


def transform_lnsp2sp_LESM_files(namelist, SWL_info, member):
      
      """
      Transform logarithm of surface pressure into surface pressure file

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      member (int)   : member number

      """

      for case in SWL_info.keys():
          
          os.system('cdo -s expr,"var134=exp(var152)" '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var152.grb '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var134.grb.tmp ')
          os.system('grib_set -s table2Version=128 '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var134.grb.tmp '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var134.grb ')


def create_pres_LESM_files(namelist, SWL_info, member):
      
      """
      Create 3D pressure field using surface pressure and hybrid levels

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      member (int)   : member number

      """

      for case in SWL_info.keys():
          
          os.system('cdo -s -pressure_fl '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+' '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var54.grb.tmp1 ')
          os.system('grib_set -s table2Version=128 '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var54.grb.tmp1 '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var54.grb.tmp2 ')
          os.system('grib_set -s indicatorOfParameter=54 '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var54.grb.tmp2 '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var54.grb ')


def transform_q2rh_LESM_files(namelist, SWL_info, member):
      
      """
      Transform specific humidity into relative humidity

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      member (int)   : member number

      """

      for case in SWL_info.keys():
          
          os.system('cdo merge '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var133.grb '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var130.grb '
                    +namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var54.grb '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var157.grb.tmp1 ')
          
          os.system('cdo expr,"var157=0.263*var54*var133/exp((17.67*(var130-273.15))/(var130-29.65))" '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var157.grb.tmp1 '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var157.grb.tmp2 ')
          os.system('grib_set -s table2Version=128 '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var157.grb.tmp2 '+namelist['output_dir']+'tmpDIR/'+case+'_m'+str(101+member)+'_var157.grb ')


def delta_ensemble(namelist, SWL_info):

      """
      Compute ensemble of PGW deltas compute as netcdf

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      
      """

      for case in list(list(SWL_info.keys())[1:]):
           for varname in namelist['varnames_sub']:

                if varname == 'var152': continue

                os.system('cdo -s -O -f nc copy -ensmean -del29feb '+namelist['output_dir']+'HCLIM_grib/PGW'+case[3]+'/m*/'
                          +varname+'_*.grb '+namelist['output_dir']+'HCLIM_grib/PGW'+case[3]+
                          '/ensemble/'+varname+'_PGW'+case[3]+'.nc')


def delta_timeinterp(namelist, SWL_info):

      """
      Modify time interpolation of ensemble deltas

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      
      """

      for cases in list(SWL_info.keys())[1:]:
           
           case = 'PGW'+cases[3]
           refyear = namelist[case]['refyear']

           for varname in namelist['varnames_sub']:

                if varname == 'var152': continue

                # Create a temporary file
                os.system('cdo -s setyear,'+str(refyear)+' '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp0.nc')
                
                os.system('cdo -s settime,12:00:00 '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp0.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp1.nc')

                # Select date
                os.system('cdo -s seldate,'+str(refyear)+'-01-01T12:00:00 '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp1.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/tmp_jan.nc')
                os.system('cdo -s seldate,'+str(refyear)+'-12-31T12:00:00 '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp1.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/tmp_dec.nc')
                
                # Change data 
                os.system('cdo -s setyear,'+str(refyear+1)+' '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/tmp_jan.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/jan.nc')
                os.system('cdo -s setyear,'+str(refyear-1)+' '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/tmp_dec.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/dec.nc')

                # merge files
                os.system('cdo -s mergetime '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/dec.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp1.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/jan.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp2.nc')

                # Interpolate times
                os.system('cdo -s inttime,'+str(refyear-1)+'-12-31,12:00:00,3hour '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp2.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp3.nc')

                os.system('cdo -s seldate,'+str(refyear)+'-01-01T00:00:00,'+str(refyear+1)+'-01-01T00:00:00 '+
                           namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp3.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'.nc')
               
                # Clean temporary files
                os.system('rm -rf '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/jan.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/dec.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/tmp_jan.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/tmp_dec.nc '+namelist['output_dir']+'HCLIM_grib/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp*.nc')


def harmonic_ac_analysis(ts):
    
    """
    Estimation of the harmonics according to formula 12.19 -
    12.23 on p. 264 in Storch & Zwiers
    
    Arguments:
        ts: a 1-d numpy array of a timeseries
    Returns:
        hcts: a reconstructed smoothed timeseries 
                (the more modes are summed the less smoothing)
        mean: the mean of the timeseries (needed for reconstruction)
    """

    # Substitue nan values with one    
    ts = np.nan_to_num(ts, nan=1)

    #calculate the mean of the timeseries (used for reconstruction)
    mean = ts.mean() 

    lt = len(ts) #how long is the timeseries?
    P = lt

    #initialize the output array. 
    #we will use at max 4 modes for reconstruction 
    #(for peformance reasons, it can be increased)
    hcts = np.zeros((4,lt))

    timevector=np.arange(1,lt+1,1)	#timesteps used in calculation	

    #a measure that is to check that the performed calculation 
    # is justified.
    q = math.floor(P/2.) 

    #create the reconstruction timeseries, mode by mode 
    #(starting at 1 until 5, if one wants more smoothing 
    #this number can be increased.)
    for i in range(1,4): 
        if i < q: #only if this is true the calculation is valid

            #these are the formulas from Storch & Zwiers
            bracket = 2.*math.pi*i/P*timevector
            a = 2./lt*(ts.dot(np.cos(bracket))) 
            #dot product (Skalarprodukt) for scalar number output!
            b = 2./lt*(ts.dot(np.sin(bracket))) 

            #calculate the reconstruction time series
            hcts[i-1,:] = a * np.cos(bracket) + b * np.sin(bracket) 

        else: #abort if the above condition is not fulfilled. In this case more programming is needed.
            sys.exit('Whooops that should not be the case for a yearly '+
            'timeseries! i (reconstruction grade) is larger than '+
            'the number of timeseries elements / 2.')

    smooths = sum(hcts[0:3,:]) + mean
    return smooths


def smooth(x, dim):

    """
    Apply smoothing of an annual timeseries 
    (typically daily resolution) using a spectral filter 
    (Bosshard et al. 2011).

    Parameters:
    x (float): variable to smooth

    Returns:
    float : smoothed variable
    """

    if dim == 2:
        return smooth_2d(x)
    elif dim == 3:
        return smooth_3d(x)

def smooth_2d(x):

    """
    Apply smoothing of an annual timeseries 
    (typically daily resolution) using a spectral filter 
    (Bosshard et al. 2011).

    Parameters:
    x (float): variable to smooth

    Returns:
    float : smoothed variable
    """

    n_y = x.shape[1]
    n_x = x.shape[2]


    for jj in range(n_y):
        for ii in range(n_x):

            x[:,jj,ii] = harmonic_ac_analysis(x[:,jj,ii])

    return x


def smooth_3d(x):

    """
    Apply smoothing of an annual timeseries 
    (typically daily resolution) using a spectral filter 
    (Bosshard et al. 2011).

    Parameters:
    x (float): variable to smooth

    Returns:
    float : smoothed variable
    """

    n_z = x.shape[1]
    n_y = x.shape[2]
    n_x = x.shape[3]


    for kk in range(n_z):
       for jj in range(n_y):
            for ii in range(n_x):

                x[:,kk,jj,ii] = harmonic_ac_analysis(x[:,kk,jj,ii])

    return x
 
def delta_smooth(namelist, SWL_info):

      """
      Modify time interpolation of ensemble deltas

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      
      """


      for case in list(SWL_info.keys())[1:]:

           for varname in namelist['varnames_sub']:

                print(varname)
                if varname == 'var152': continue

                filename = namelist['output_dir']+'HCLIM_grib/PGW'+case[3]+'/ensemble/'+varname+'_PGW'+case[3]+'.nc'

                # Apply filter
                nc = Dataset(filename,'r+')
                try:
                    var_nc = nc.variables[varname][:,:,:,:]
                    nc.variables[varname][:,:,:,:] = smooth(var_nc,3)

                except:
                    var_nc = nc.variables[varname][:,:,:]
                    nc.variables[varname][:,:,:] = smooth(var_nc,2)

                nc.close()


def get_varinfo():

    varinfo = {}
    varinfo['var34']  ={'varname':'delta_sst','long_name':'delta for sea surface temperature', 'units': 'K'}
    varinfo['var130'] ={'varname':'delta_t','long_name':'delta for T', 'units': 'K'}
    varinfo['var131'] ={'varname':'delta_u','long_name':'delta for zonal velocity', 'units': 'm/s'}
    varinfo['var132'] ={'varname':'delta_v','long_name':'delta for meridional velocity', 'units': 'm/s'}
    varinfo['var133'] ={'varname':'delta_q','long_name':'delta for specific humidity', 'units': 'kg/kg'}
    varinfo['var134'] ={'varname':'delta_sp','long_name':'delta for surface pressure', 'units': 'Pa'}
    varinfo['var157'] ={'varname':'delta_r','long_name':'delta for relative humidity', 'units': '%'}
    varinfo['var235'] ={'varname':'delta_skt','long_name':'delta for surface temperatures', 'units': 'K'}

    return varinfo


def delta_rewritenc(namelist, SWL_info):

      """
      Modify the metadata of the netcdf files

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      
      """

      vardict = get_varinfo()

      for case in list(SWL_info.keys())[1:]:

           for varname in namelist['varnames_sub']:
               

               if varname == 'var152': continue

               filein = namelist['output_dir']+'HCLIM_grib/PGW'+case[3]+'/ensemble/'+varname+'_PGW'+case[3]+'.nc'
               fileout = namelist['output_dir']+'HCLIM_grib/PGW'+case[3]+'/ensemble/'+vardict[varname]['varname']+'_PGW'+case[3]+'.nc'

               nvarname = vardict[varname]['varname']
            
               os.system('ncrename -v '+varname+','+nvarname+' '+filein+' '+fileout)
               os.system('ncatted -O -a long_name,'+nvarname+',o,c,"'+vardict[varname]['long_name']+'" -a units,'+nvarname+',o,c,"'+vardict[varname]['units']+'" '+fileout)

# -----------------------------------------------------------------------------

# 1 - Read namelist
namelist = read_namelist('../nam/namelist_deltas_HCLIM_grib.yaml')

# 2 - Read SWL info
SWL_info = read_SWL_info(namelist,['SWL09','SWL20','SWL30'])

# 3 - Create temporary work directory
create_workdir(namelist)

# 4 - Loop over all the members
for member in tqdm(range(namelist['num_members']), desc='Loop over members', leave=True):
  
   # if data already extracted skipped this step
   if not namelist['extract_data']: break

   # Create member directories
   create_workdir(namelist,member,['PGW2','PGW3'])
 
   # Extract data from LESM
   extract_LESM_files(namelist, SWL_info, member)

   # Separate grib files per variable
   separate_LESM_files_variables(namelist, SWL_info, member)

   # Create grib file for surface pressure
   transform_lnsp2sp_LESM_files(namelist, SWL_info, member)

   # Create grib file for pressure at all levels
   create_pres_LESM_files(namelist, SWL_info, member)

   # Create grib file for relative humidity
   transform_q2rh_LESM_files(namelist, SWL_info, member)

   # Compute deltas - substraction  
   delta_substract_LESM(namelist, SWL_info, member)

   # Archive delta information
   delta_archive(namelist, SWL_info, member)

   # Clean working directory
   clean_workdir(namelist, member)

# 5 - Compute ensemble of all deltas
delta_ensemble(namelist, SWL_info)

# 6 - Smooth deltas
if namelist['smooth']: delta_smooth(namelist, SWL_info)

# 7 - Adapt deltas -> interpolate in time
delta_timeinterp(namelist, SWL_info)

# 8 - Rewrite metadata in netcdf files
delta_rewritenc(namelist, SWL_info)

# 9 - Remove working directory
remove_workdir(namelist)

# ----------------------------------------------------------------------------------
