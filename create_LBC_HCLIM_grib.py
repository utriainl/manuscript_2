import pygrib
import yaml
import numpy as np
import os
import datetime

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

    with open(filename) as f:
        namelist = yaml.load(f, Loader=SafeLoader)
	
    return namelist


def num_days_in_month(year, month):
	
    daysinmonth = [29,31,28,31,30,31,30,31,31,30,31,30,31]

    if is_year_leap(year) and month==2:
          return daysinmonth[0]
    else:
          return daysinmonth[month]
    
def read_ERA5_ab(filename):
   
    #Load constant data for vertical levels: ERA5 a & b on half-levels, PGW pressure levels
    levm = np.loadtxt(filename)
    alevm = levm[:,0]
    blevm = levm[:,1]

    return alevm, blevm

def read_delta_ab(namelist):
    
    # Read sp deltas
    nc_file = Dataset(namelist['dir_deltas']+'/PGW2/ensemble/delta_sp_PGW2_3h.nc','r')
    am_delta = nc_file.variables['hyam'][:].data
    bm_delta = nc_file.variables['hybm'][:].data
    nc_file.close()

    return am_delta, bm_delta

def read_delta_var(namelist,pgw,tstep,varname):
    
    # Read sp deltas
    nc_file = Dataset(namelist['dir_deltas']+'/'+pgw+'/ensemble/delta_'+varname+'_'+pgw+'_3h.nc','r')
    if varname == 'sp':
      delta_var = nc_file.variables['delta_'+varname][tstep,0,:,:].data
    elif varname == 'skt' or varname == 'sst':
      delta_var = nc_file.variables['delta_'+varname][tstep,:,:].data
    else:
      delta_var = nc_file.variables['delta_'+varname][tstep,:,:,:].data

    nc_file.close()

    return delta_var

def read_ablevels(namelist,nx,ny):
   
    am_era5, bm_era5 = read_ERA5_ab(namelist['level_file'])

    am_delta, bm_delta = read_delta_ab(namelist)

    am_delta = np.expand_dims(am_delta,axis=(1,2))
    am_delta = np.tile(am_delta, (1,nx,ny))
    bm_delta = np.expand_dims(bm_delta,axis=(1,2))
    bm_delta = np.tile(bm_delta, (1,nx,ny))


    am_era5 = np.expand_dims(am_era5,axis=(1,2))
    am_era5 = np.tile(am_era5, (1,nx,ny))
    bm_era5 = np.expand_dims(bm_era5,axis=(1,2))
    bm_era5 = np.tile(bm_era5, (1,nx,ny))

    read_ab = False

    return am_era5, bm_era5, am_delta, bm_delta, read_ab
            
def read_delta_vinterp(namelist, pgw, tstep, p_era5, p_delta, case):
                  
    delta_var = read_delta_var(namelist, pgw, tstep, case)

    # Interpolate delta value
    delta_var_interp = np.zeros(p_era5.shape) * np.nan
    for jj in range(p_era5.shape[1]):
      for ii in range(p_era5.shape[2]):
        delta_var_interp[:,jj,ii] = np.interp(np.squeeze(p_era5[:,jj,ii]), p_delta[:,jj,ii], np.squeeze(delta_var[:,jj,ii]))

    return delta_var_interp


def transform_q2rh(hus, ta, pa):
    """
    Compute relative humidity from specific humidity.
    """
    hur = 0.263 * pa * hus *(np.exp(17.67*(ta - 273.15)/(ta-29.65)))**(-1)
    return(hur)


def transform_rh2q(hur, ta, pa):
    """
    Compute specific humidity from relative humidity.
    """
    hus = (hur  * np.exp(17.67 * (ta - 273.15)/(ta - 29.65))) / (0.263 * pa)
    return(hus)

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

def calculate_DOY(year, month, day):
    
    # Remove leap day, DOY same as 28th February
    if is_year_leap(year) and month ==2 and day==29:
       
      day_of_year = (datetime.datetime(1999, 2, 28) - datetime.datetime(1999, 1, 1)).days
    
    else:
   
      day_of_year = (datetime.datetime(1999, month, day) - datetime.datetime(1999, 1, 1)).days

    return day_of_year

def calculate_tstep(year, month, day, hh):
    
    doy = calculate_DOY(year,month,day)

    return int(8*doy + hh/3)


# ------------------------------------------------------------------------------

# 1 - Read namelist
namelist = read_namelist('../nam/namelist_create_LBC_HCLIM_grib.yaml')

read_ab = True

# 2 - Loop over dates
for yy in namelist['years']:
  for mm in [3]:#namelist['months']:
    for dd in [1]:#range(1,num_days_in_month(yy,mm)+1):

      for hh in [0]:#range(0,22,3): 

        tstep = calculate_tstep(yy,mm,dd,hh)
          
        pgw1fname = namelist['input_dir']+'/'+str(yy)+'/ERA5_NEUR_{0:04d}{1:02d}{2:02d}{3:02d}00+000H00M'.format(yy,mm,dd,hh)
        
        tempera5fname = namelist['template_dir']+namelist['template_file']

        if os.path.isfile(pgw1fname):

          # Open grib file
          grbs_era = pygrib.open(pgw1fname)
          grbs_tmp = pygrib.open(tempera5fname)

          grb_var = grbs_era.message(703)

          lnsp = grb_var.values[:,:]
          sp   = np.exp(lnsp)

          # Read a/b levels
          if read_ab: am_era5, bm_era5, am_delta, bm_delta, read_ab = read_ablevels(namelist, sp.shape[0],sp.shape[1])     

          for pgw in ['PGW2','PGW3']:

            delta_sp = read_delta_var(namelist,pgw,tstep,'sp')

            nsp = sp + delta_sp

            # ERA5 pressure levels
            p_era5 = am_era5 + bm_era5*nsp

            # Delta pressure levels
            p_delta = am_delta + bm_delta*nsp

            nera5fname = namelist['output_dir']+'/'+str(yy)+'/ERA5_NEUR_'+pgw+'_{0:04d}{1:02d}{2:02d}{3:02d}00+000H00M'.format(yy,mm,dd,hh)

            print(tstep, nera5fname)

            grbs_nera = open(nera5fname,'wb')

            #for igrb in tqdm(range(1,grbs_era.messages+1), desc='Variables'):
            for igrb in range(1,grbs_era.messages+1):

                grb_var = grbs_era.message(igrb)
                grb_tmp = grbs_tmp.message(igrb)

                if igrb == 3:

                  delta_sst = read_delta_var(namelist,pgw,tstep,'sst')

                  grb_tmp.values = grb_var.values + delta_sst

                elif igrb in [11, 13, 15, 16, 17]:

                  delta_ts = read_delta_var(namelist,pgw,tstep,'skt')

                  grb_tmp.values = grb_var.values + delta_ts

                elif igrb>=18 and igrb<=154:
                   
                  if igrb == 18: delta_var_interp = read_delta_vinterp(namelist, pgw, tstep, p_era5, p_delta, 't')

                  if igrb == 18:
                     tem = np.zeros((p_era5.shape))
                     ntem = np.zeros((p_era5.shape))

                  ilevel = grb_var.level
                  grb_tmp.values = grb_var.values + delta_var_interp[ilevel-1,:,:]

                  tem[ilevel-1,:,:] = grb_var.values
                  ntem[ilevel-1,:,:] = grb_tmp.values

                elif igrb>=155 and igrb<=291:
                   
                  if igrb == 155: delta_var_interp = read_delta_vinterp(namelist, pgw, tstep, p_era5, p_delta, 'r')

                  ilevel = grb_var.level

                  rh = transform_q2rh(grb_var.values, tem[ilevel-1,:,:], am_era5[ilevel-1,:,:]+bm_era5[ilevel-1,:,:]*sp)

                  nrh = rh + delta_var_interp[ilevel-1,:,:]

                  nq = transform_rh2q(nrh, ntem[ilevel-1,:,:],p_era5[ilevel-1,:,:])

                  grb_tmp.values = nq

                  if igrb == 291:
                     del tem 
                     del ntem

                elif igrb>=292 and igrb<=428:
                   
                  if igrb == 292: delta_var_interp = read_delta_vinterp(namelist, pgw, tstep, p_era5, p_delta, 'u')

                  ilevel = grb_var.level
                  grb_tmp.values = grb_var.values + delta_var_interp[ilevel-1,:,:]

                elif igrb>=429 and igrb<=565:
                   
                  if igrb == 429: delta_var_interp = read_delta_vinterp(namelist, pgw, tstep, p_era5, p_delta, 'v')

                  ilevel = grb_var.level
                  grb_tmp.values = grb_var.values + delta_var_interp[ilevel-1,:,:]

                elif igrb == 703:
                   
                  grb_tmp.values = np.log(nsp)

                else:

                  grb_tmp.values = grb_var.values

                grb_tmp.dataDate = yy*10000+mm*100+dd

                msg = grb_tmp.tostring()
                grbs_nera.write(msg)

            grbs_nera.close()  

          grbs_era.close()
          grbs_tmp.close()
