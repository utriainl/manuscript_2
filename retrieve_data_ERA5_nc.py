import cdsapi
import numpy as np
import os.path
import os
import yaml
from yaml.loader import SafeLoader


# ------------------------------------------------------------------------------------------------------------------

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


def create_workdir(namelist):

	"""
	Create work directories 

	Parameters:
	namelist (dict): dictionary containing paths and names

	"""

	os.system('mkdir -p '+namelist['output_dir']+'tmpDIR')


def retrieve_nonaccumulate_data(year, filename):

	"""
	Retrieve non accumulate variables from ERA5

	Parameters:
        year(int)      : year of interest
	filename (str) : name of the output file

	"""

	c = cdsapi.Client()
	c.retrieve(
		'reanalysis-era5-single-levels',
		{
		'product_type': 'reanalysis',
		'variable': [
			'10m_u_component_of_neutral_wind', '10m_v_component_of_neutral_wind', '2m_dewpoint_temperature',
			'2m_temperature', 'geopotential', 'surface_pressure',
			],
		'year': str(year),
		'month': [
			'01', '02', '03',
			'04', '05', '06',
			'07', '08', '09',
			'10', '11', '12',
			],
		'day': [
			'01', '02', '03',
			'04', '05', '06',
			'07', '08', '09',
			'10', '11', '12',
			'13', '14', '15',
			'16', '17', '18',
			'19', '20', '21',
			'22', '23', '24',
			'25', '26', '27',
			'28', '29', '30',
			'31',
			],
		'time': [
			'00:00', '03:00', '06:00',
			'09:00', '12:00', '15:00',
			'18:00', '21:00',
			],
		'area': [
			75, -10, 45,
			45,
			],
		'format': 'netcdf',
		},
		filename)


def retrieve_accumulate_data(year, filename):

	"""
	Retrieve hourly accumulated variables from ERA5

	Parameters:
        year(int)      : year of interest
	filename (str) : name of the output file

	"""

	c = cdsapi.Client()
	c.retrieve(
		'reanalysis-era5-single-levels',
		{
		'product_type': 'reanalysis',
		'variable': [
    			'snowfall', 'surface_solar_radiation_downwards', 'surface_thermal_radiation_downwards',
    			'total_precipitation',
			],
		'year': str(year),
		'month': [
    			'01', '02', '03',
    			'04', '05', '06',
    			'07', '08', '09',
    			'10', '11', '12',
			],
		'day': [
			'01', '02', '03',
			'04', '05', '06',
			'07', '08', '09',
			'10', '11', '12',
			'13', '14', '15',
			'16', '17', '18',
			'19', '20', '21',
			'22', '23', '24',
			'25', '26', '27',
			'28', '29', '30',
			'31',
			],
		'time': [
			'00:00', '01:00', '02:00', '03:00', 
			'04:00', '05:00', '06:00', '07:00',
			'08:00', '09:00', '10:00', '11:00',
			'12:00', '13:00', '14:00', '15:00',
			'16:00', '17:00', '18:00', '19:00',
			'20:00', '21:00', '22:00', '23:00',
			],
		'area': [
			75, -10, 45,
			45,
		],
		'format': 'netcdf',
		},
		filename)


# ----------------------------------------------------------------------------------------

# 1 - Read namelist file
namelist = read_namelist('../nam/namelist_retrieve_data_ERA5.yaml')

# 2 - Retrieve non accumulated fields
outnc_nonaccu = namelist['output_dir']+'ERA5_NONACCU_'+str(namelist['year'])+'.nc'
if not os.path.isfile(outnc_nonaccu): 
	print('File '+outnc_nonaccu+' not found ... retrieving it')
	retrieve_nonaccumulate_data(namelist['year'],outnc_nonaccu)

# 3 - Retrieve accumulated fields
outnc_accu = namelist['output_dir']+'ERA5_ACCU_'+str(namelist['year'])+'.nc'
if not os.path.isfile(outnc_accu):
	print('File '+outnc_accu+' not found ... retrieving it')
	retrieve_accumulate_data(namelist['year'],outnc_accu)

# 4 - Remap to domain of interest
if namelist['remap']:

	os.system('cdo remapbil,'+namelist['grid_info']+' '+namelist['output_dir']+'ERA5_ACCU_'+str(namelist['year'])+'.nc '+namelist['output_dir']+'ERA5_ACCU_'+namelist['domain']+'_'+str(namelist['year'])+'.nc')
	os.system('cdo remapbil,'+namelist['grid_info']+' '+namelist['output_dir']+'ERA5_NONACCU_'+str(namelist['year'])+'.nc '+namelist['output_dir']+'ERA5_NONACCU_'+namelist['domain']+'_'+str(namelist['year'])+'.nc')    
	

