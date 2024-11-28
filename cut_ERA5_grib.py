import pygrib
import yaml
import numpy as np
import os

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
    

# ------------------------------------------------------------------------------

# 1 - Read namelist
namelist = read_namelist('../nam/namelist_cut_ERA5_grib.yaml')

# 2 - Loop over dates
for yy in namelist['years']:
  for mm in namelist['months']:
    for dd in range(1,num_days_in_month(yy,mm)+1):
      for hh in range(0,22,3): 

        era5fname = namelist['input_dir']+'ERA5_{0:04d}{1:02d}{2:02d}{3:02d}00+000H00M'.format(yy,mm,dd,hh)

        nera5fname = namelist['output_dir']+'ERA5_NEUR_{0:04d}{1:02d}{2:02d}{3:02d}00+000H00M'.format(yy,mm,dd,hh)

        tempera5fname = namelist['template_dir']+namelist['template_file']

        if os.path.isfile(era5fname):

          # Open grib file
          grbs_era = pygrib.open(era5fname)
          grbs_tmp = pygrib.open(tempera5fname)
          grbs_nera = open(nera5fname,'wb')

          for igrb in tqdm(range(1,grbs_era.messages+1), desc='Variables'):
               

            grb_var = grbs_era.message(igrb)

            new_var = np.empty((84,154))
            new_var[:,:17] = grb_var.values[60:144, -17:]
            new_var[:,17:] = grb_var.values[60:144, :137]

            grb_tmp = grbs_tmp.message(igrb) 
            grb_tmp.values = new_var
            grb_tmp.dataDate = yy*10000+mm*100+dd

            msg = grb_tmp.tostring()
            grbs_nera.write(msg)

          # Close files
          grbs_era.close()
          grbs_tmp.close()
          grbs_nera.close()
