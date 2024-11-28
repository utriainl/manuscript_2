import numpy as np
import os
import sys
import yaml
import pandas as pd
import math
from netCDF4 import Dataset

from tqdm import tqdm
from yaml.loader import SafeLoader

# ------------------------------------------------------------------------------

def slp2sp(p0,th,hh):

    """
    Tranform sea level pressure to surface pressure

    Parameters:
    p0 (3D array) : Sea level pressure
    th (3D array) : Surface air pressure
    hh (2D array) : Altitude

    Returns:
    3D array : surface pressure

    """

    hh = np.repeat(hh[np.newaxis,:,:], p0.shape[0],axis=0)

    M = 0.02896968  # Molar mass of dry air (kg/mol)
    R = 8.314462618 # Universal gas constant (J/(mol K))
    g = 9.81        # gravitational acceleration (m/s2)

    return p0*np.exp(- g*M*hh/(R*th))


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
            os.system('mkdir -p '+namelist['output_dir']+'SURFEX_nc/'+case+'/m'+str(member+101))


def remove_workdir(namelist):

    """
    Remove working directory

    Parameters:
    namelist (dict): dictionary containing path information

    """

    os.system('rm -rf '+namelist['output_dir']+'tmpDIR')


def clean_workdir(namelist):

    """
    Clean working directory by removing all the temporary netCDF files

    Parameters:
    namelist (dict): dictionary containing path information
    
    """


    os.system('rm -rf '+namelist['output_dir']+'tmpDIR/*.nc')


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

    for case in SWL_info.keys():
        for varname in namelist['varnames_sub']+namelist['varnames_div']+namelist['varnames_other']:

                if varname == 'sp': continue

                yearref=SWL_info[case][member]
                yearlist=np.arange(yearref-yearshift,yearref+1)

                str1 = main_path.replace("CASECASE", str(namelist[case]['case']))
                str2 = str1.replace('EXPEXP', namelist[case]['exp'])

                if varname in namelist['varnames_other']:
                    str1 = str2.replace('day','Amon')
                    path_files = str1.replace('MMM',str(101+member))

                    forget_29feb = False
                else:
                    path_files = str2.replace('MMM',str(101+member))
                    forget_29feb = True

                for year in yearlist:
      
                   os.system('cdo -s -remapbil,'+namelist['grid_info']+' '+path_files+varname+'/gr/'+
                             namelist[case]['version']+'/*_'+str(year)+'*nc '+namelist['output_dir']+
                             'tmpDIR/'+varname+'_'+str(101+member)+'_'+case+'_'+str(year)+'.nc')

                   if is_year_leap(year) and forget_29feb:
                        
                        os.system('mv '+namelist['output_dir']+'tmpDIR/'+varname+'_'+
                                  str(101+member)+'_'+case+'_'+str(year)+'.nc '+
                                  namelist['output_dir']+'tmpDIR/temp.nc')

                        os.system('ncks -F -d time,1,59 -d time,61,366 '+
                                  namelist['output_dir']+'tmpDIR/temp.nc '+
                                  namelist['output_dir']+'tmpDIR/'+varname+
                                  '_'+str(101+member)+'_'+case+'_'+str(year)+'.nc')


def transform_slp2sp_LESM_files(namelist, SWL_info, member):
      
      """
      Transform sea level pressure file into surface pressure file

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      member (int)   : member number

      """

      for case in SWL_info.keys():

          yearref=SWL_info[case][member]
          yearlist=np.arange(yearref-5,yearref+1)

          for year in yearlist:

             os.system('cp -r '+namelist['output_dir']+'tmpDIR/psl_'+
                      str(101+member)+'_'+case+'_'+str(year)+'.nc '+namelist['output_dir']+'tmpDIR/sp_'+
                      str(101+member)+'_'+case+'_'+str(year)+'.nc')

             os.system('ncrename -v psl,sp '+namelist['output_dir']+'tmpDIR/sp_'+
                      str(101+member)+'_'+case+'_'+str(year)+'.nc > /dev/null 2>&1 ')

             sp_file = namelist['output_dir']+'tmpDIR/sp_'+str(101+member)+'_'+case+'_'+str(year)+'.nc'
             tas_file = namelist['output_dir']+'tmpDIR/tas_'+str(101+member)+'_'+case+'_'+str(year)+'.nc'
             orog_file = namelist['oro_dir']+namelist['oro_file']

             # Transform sea level pressure to surface pressure
             rewrite_sp_file(sp_file, tas_file, orog_file)   

             os.system('rm -rf '+namelist['output_dir']+'tmpDIR/psl_'+
                       str(101+member)+'_'+case+'_'+str(year)+'.nc')


def rewrite_sp_file(sp_file, tas_file, orog_file):

      """
      Rewrite surface pressure file, information and metadata

      Parameters:
      sp_file (str)   : path to the surface pressure file
      tas_file (str)  : path to the surace air temperature file
      orog_file (str) : path to the orography file

      """

      # Read orography
      nc = Dataset(orog_file,'r')
      hh = nc.variables['orog'][:,:]
      nc.close()

      # Read surface temperature
      nc = Dataset(tas_file,'r')
      tas = nc.variables['tas'][:,:,:]
      nc.close()

      # Rewrite surface pressure
      nc = Dataset(sp_file,'r+')
      sp = nc.variables['sp'][:,:,:]
      nc.variables['sp'][:,:,:] = slp2sp(sp,tas,hh) 
      nc.variables['sp'].long_name = "Surface pressure"
      nc.variables['sp'].standard_name = "surface_air_pressure"
      nc.variables['sp'].comment ="Surface pressure"
      nc.close()


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

           yearref=SWL_info['SWL09'][member]
           yearlist0=np.arange(yearref-5,yearref+1)

           yearref=SWL_info[case][member]
           yearlist1=np.arange(yearref-5,yearref+1)

           for varname in namelist['varnames_sub']+namelist['varnames_other']:
      
                if varname == 'psl': continue

                for year0,year1 in zip(yearlist0,yearlist1):
                       os.system('cdo -s sub '+namelist['output_dir']+'tmpDIR/'+varname+
                             '_'+str(101+member)+'_'+case+'_'+str(year1)+'.nc '+
                             namelist['output_dir']+'tmpDIR/'+varname+
                             '_'+str(101+member)+'_'+case0+'_'+str(year0)+'.nc '+
                             namelist['output_dir']+'tmpDIR/'+'/'+varname+'_'+str(101+member)+
                             '_'+case+'-'+case0+'_m'+str(yearref-year1)+'yr.nc ')


def delta_fraction_LESM(namelist, SWL_info, member):

      """
      Compute PGW deltas by dividing two SWL

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      member (int)   : member number

      """

      case0 = list(SWL_info.keys())[0]

      for case in list(list(SWL_info.keys())[1:]):

           yearref=SWL_info['SWL09'][member]
           yearlist0=np.arange(yearref-5,yearref+1)

           yearref=SWL_info[case][member]
           yearlist1=np.arange(yearref-5,yearref+1)

           for varname in namelist['varnames_div']:
               
               for year0,year1 in zip(yearlist0,yearlist1):

                       os.system('ncbo -v '+varname+' -y divide '+namelist['output_dir']+'tmpDIR/'+varname+
                             '_'+str(101+member)+'_'+case+'_'+str(year1)+'.nc '+
                             namelist['output_dir']+'tmpDIR/'+varname+
                             '_'+str(101+member)+'_'+case0+'_'+str(year0)+'.nc -O '+
                             namelist['output_dir']+'tmpDIR/'+varname+'_tmp.nc ')

                       os.system('ncap2 -s "where('+varname+'>1.e10) '+varname+'=0." '
                                 +namelist['output_dir']+'tmpDIR/'+varname+'_tmp.nc -O '
                                 +namelist['output_dir']+'tmpDIR/'+varname+'_'+str(101+member)+
                                 '_'+case+'-'+case0+'_m'+str(yearref-year1)+'yr.nc ')

                       os.system('rm -rf '+namelist['output_dir']+'tmpDIR/'+varname+'_tmp.nc ')


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
           for varname in namelist['varnames_sub']+namelist['varnames_div']+namelist['varnames_other']:

               if varname == 'psl': continue

               os.system('cdo -s -O mergetime '+namelist['output_dir']+'tmpDIR/'+varname+'*_'+case+'-'+
                         case0+'_m*.nc '+namelist['output_dir']+'SURFEX_nc/PGW'+case[3]+'/m'+
                         str(member+101)+'/'+varname+'_tmp.nc')

               os.system('cdo -s -O ydaymean '+namelist['output_dir']+'SURFEX_nc/PGW'+case[3]+'/m'+
                         str(member+101)+'/'+varname+'_tmp.nc '+namelist['output_dir']+'SURFEX_nc/PGW'+
                         case[3]+'/m'+str(member+101)+'/'+varname+'_'+str(101+member)+
                         '_'+case+'-'+case0+'.nc')

               os.system('rm -rf '+namelist['output_dir']+'SURFEX_nc/PGW'+case[3]+'/m'+
                         str(member+101)+'/'+varname+'_tmp.nc')

               
def delta_ensemble(namelist, SWL_info):

      """
      Compute ensemble of PGW deltas

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      
      """

      for case in list(list(SWL_info.keys())[1:]):
           for varname in namelist['varnames_sub']+namelist['varnames_div']+namelist['varnames_other']:

                if varname == 'psl': continue

                os.system('cdo -s -O ensmean '+namelist['output_dir']+'SURFEX_nc/PGW'+case[3]+'/m*/'
                          +varname+'_*.nc '+namelist['output_dir']+'SURFEX_nc/PGW'+case[3]+
                          '/ensemble/'+varname+'_PGW'+case[3]+'.nc')


def delta_timeinterp(namelist, SWL_info):

      """
      Modify time interpolation of ensemble deltas

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      
      """

      for case in list(SWL_info.keys())[1:]:
           
           refyear = namelist['PGW'+case[3]]['refyear']

           for varname in namelist['varnames_sub']:

                if varname == 'psl': continue 

                print(varname)

                # Create a temporary file
                os.system('cdo -s setyear,'+str(refyear)+' '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp1.nc')

                # Select date
                os.system('cdo -s seldate,'+str(refyear)+'-01-01T12:00:00 '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp1.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/tmp_jan.nc')
                os.system('cdo -s seldate,'+str(refyear)+'-12-31T12:00:00 '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp1.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/tmp_dec.nc')
                
                # Change data 
                os.system('cdo -s setyear,'+str(refyear+1)+' '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/tmp_jan.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/jan.nc')
                os.system('cdo -s setyear,'+str(refyear-1)+' '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/tmp_dec.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/dec.nc')

                # merge files
                os.system('cdo -s mergetime '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/dec.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp1.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/jan.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp2.nc')

                # Interpolate times
                os.system('cdo -s inttime,'+str(refyear-1)+'-12-31,12:00:00,3hour '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp2.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp3.nc')

                os.system('cdo -s seldate,'+str(refyear)+'-01-01T00:00:00,'+str(refyear+1)+'-01-01T00:00:00 '+
                           namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'_tmp3.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/'+varname+'_'+case+'.nc')
               
                # Clean temporary files
                os.system('rm -rf '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/jan.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/dec.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/tmp_jan.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
                          '/ensemble/tmp_dec.nc '+namelist['output_dir']+'SURFEX_nc/'+case+
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


def smooth(x):

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
 
def delta_smooth(namelist, SWL_info):

      """
      Modify time interpolation of ensemble deltas

      Parameters:
      namelist (dict): dictionary containing paths and orography files
      SWL_info (dict): SWL information
      
      """


      for case in list(SWL_info.keys())[1:]:

           for varname in namelist['varnames_sub']+namelist['varnames_div']:

                if varname == 'psl': continue

                filename = namelist['output_dir']+'SURFEX_nc/PGW'+case[3]+'/ensemble/'+varname+'_PGW'+case[3]+'.nc'

                # Apply filter
                nc = Dataset(filename,'r+')
                var_nc = nc.variables[varname][:,:,:]
                nc.variables[varname][:,:,:] = smooth(var_nc)
                nc.close()



# -----------------------------------------------------------------------------

# 1 - Read namelist
namelist = read_namelist('../nam/namelist_deltas_SURFEX_nc.yaml')

# 2 - Read SWL info
SWL_info = read_SWL_info(namelist,['SWL09','SWL20','SWL30'])

# 3 - Create temporary work directory
create_workdir(namelist)

# 4 - Loop over all the members
for member in tqdm(range(namelist['num_members']), desc='Loop over members', leave=True):
  
   # Create member directories
   create_workdir(namelist,member,['PGW2','PGW3'])
 
   # Extract data from LESM
   extract_LESM_files(namelist, SWL_info, member)

   # Transform slp to sp
   transform_slp2sp_LESM_files(namelist, SWL_info, member)   

   # Compute deltas - substraction  
   delta_substract_LESM(namelist, SWL_info, member)

   # Compute deltas - fraction
   delta_fraction_LESM(namelist, SWL_info, member)   

   # Archive delta information
   delta_archive(namelist, SWL_info, member)

   # Clean working directory
   clean_workdir(namelist)

# 5 - Compute ensemble of all deltas
delta_ensemble(namelist, SWL_info)

# 6 - Adapt deltas -> interpolate in time
delta_timeinterp(namelist, SWL_info)

# 7 - Smooth deltas
if namelist['smooth']: delta_smooth(namelist, SWL_info)

# 8 - Remove working directory
remove_workdir(namelist)

# ----------------------------------------------------------------------------------
