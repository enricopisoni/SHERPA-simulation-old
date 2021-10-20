'''
Created on Jun 23, 2015

Module 6 calculates for 1 cell the concentration change due to a 50 percent reductions in 
the snap sectors defined in the input file 'path_reduction_txt'. Emission are reduced in
each NUTS area in the input file 'path_area_cdf'
There are 2 outputs:
- a text file with all nuts codes and the DC/C/alpha (relative potential) as percent due to a reduction in that nuts area
- a map where each nuts has the value of the concentration change it causes in the target cell 

for compatibility the header is 'potency' in the output txt

@author: degraba
'''

# imports
from netCDF4 import Dataset
from numpy import lib, zeros, sum, power, array
from math import isnan
# path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_model_cdf_test,
import sys
from sherpa_globals import alpha_potency
from sherpa_auxiliaries import create_emission_reduction_dict, create_emission_dict, create_window, deltaNOx_to_deltaNO2
#EP 20210518
from sherpa_globals import sector_lst

# function that applies reductions per snap sector and precursor to the emission netcdf
def create_delta_emission(path_emission_cdf, precursor_lst, reduction_area_array, path_reduction_txt):
        
    # create a dictionary with reductions per precursor and macro sector
    emission_reduction_dict = create_emission_reduction_dict(path_reduction_txt)
    
    # open the emission netcdf
    emission_dict = create_emission_dict(path_emission_cdf, precursor_lst)
       
    # calculate a dictionary with the emission reductions per pollutant, macrosector and position
    delta_emission_dict = {}
    for precursor in precursor_lst:
        delta_emission_dict[precursor] = zeros(emission_dict[precursor].shape)
        # calculate the emission reduction
        # reductions are positive!
        # make the sum over all snap sectors
        #EP20210518
        for snap in range(1, sector_lst[-1]):
        #for snap in range(1, 13):
            delta_emission_dict[precursor][snap - 1, :, :] = emission_dict[precursor][snap - 1] * reduction_area_array * emission_reduction_dict[precursor][snap]
        
    # sum over all snap sectors
    for precursor in precursor_lst:
        delta_emission_dict[precursor] = sum(delta_emission_dict[precursor], axis=0)
              
    return delta_emission_dict

# function definition of source receptor model
def module6(path_emission_cdf, path_area_cdf, target_cell_lat, target_cell_lon, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf):
    
    # read the model netcdf
    # ---------------------
    rootgrp = Dataset(path_model_cdf, 'r')
    longitude_array = rootgrp.variables['lon'][0, :]
    latitude_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(longitude_array)  
    n_lat = len(latitude_array)

    #####
    #20180129 - EP - generalization to read 'radius of influence' variable, both written in matlab or python
    for nameatt in Dataset(path_model_cdf, 'r').ncattrs():
        if nameatt[0:6] == 'Radius':
            radiusofinfluence = nameatt
    inner_radius = int(getattr(Dataset(path_model_cdf, 'r'), radiusofinfluence))

    # inner_radius = int(getattr(rootgrp, 'Radius of influence'))
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    alpha = rootgrp.variables['alpha'][:, :, :]    
    omega = rootgrp.variables['omega'][:, :, :] 
    # put alpha and omega in a dictionary
    alpha_dict = {}
    omega_dict = {}
    for i in range(len(precursor_lst)):
        alpha_dict[precursor_lst[i]] = alpha[i, :, :]
        omega_dict[precursor_lst[i]] = omega[i, :, :]
        
    # close model netcdf
    rootgrp.close()
    
    # a better way should be found, e.g. variable inside the netcdf
    if (path_model_cdf.find('NO2eq') > -1):
        NOx = True
    else:
        NOx = False
    
    # open the area netcdf and get lat lon indexes of target cell
    #-------------------------------------------------------------
    rootgrp_nuts = Dataset(path_area_cdf, 'r')
    #20200414 EP modified for working with CAMS_EMEP
#    n_nuts = len(rootgrp_nuts.dimensions['nuts_id'])
    n_nuts = len(rootgrp_nuts.dimensions['NUTS'])
    
    #EP20210119
#    nuts_codes = rootgrp_nuts.variables['NUTS'].units.split('-')
    nuts_codes =getattr(rootgrp_nuts, 'NUTS').split('-')
    
#    nuts_codes_raw = rootgrp_nuts.variables['NUTS'][:]
#    nuts_codes = []
#    for i_code in range(len(nuts_codes_raw)):
#        code = ''
#        for letter in nuts_codes_raw[i_code]:
#            code = code + str(letter.decode('utf-8'))
#        nuts_codes.append(code)
    
    # convert latitude and longitude string in float
    target_cell_lat = float(target_cell_lat)
    target_cell_lon = float(target_cell_lon)
    
    # get row index of latitude and col index of longitude
    i_lat_target = 0
    lat_error = float('inf')
    for i in range(len(latitude_array)):
        lat_dist = abs(target_cell_lat - latitude_array[i])
        if lat_dist < lat_error:
            lat_error = lat_dist
            i_lat_target = i
    
    i_lon_target = 0
    lon_error = float('inf')
    for i in range(len(longitude_array)):
        lon_dist = abs(target_cell_lon - longitude_array[i])
        if lon_dist < lon_error:
            lon_error = lon_dist
            i_lon_target = i
    
    # read base concentrations and extract base case concentration in the target cell
    # -------------------------------------------------------------------------------
    rootgrp = Dataset(path_base_conc_cdf, 'r')
    if NOx == True:
        # in case of NO2, also get the NO2 base case concentration
        target_conc_basecase_no2 = array(rootgrp.variables['NO2'][i_lat_target, i_lon_target]) 
        target_conc_basecase_nox = rootgrp.variables['conc'][i_lat_target, i_lon_target]       # in case of NO2, conc == NOx
        # dictionary with the concentration change due to an emission reduction in a nuts, keys are nuts codes
        delta_conc_no2 = {} 
        delta_conc_nox = {}     # just in case we're doing NO2
    else:
        target_conc_basecase = rootgrp.variables['conc'][i_lat_target, i_lon_target]       # this is PM2.5 or PM10
        # dictionary with the concentration change due to an emission reduction in a nuts, keys are nuts codes
        delta_conc = {}
    # close model netcdf
    rootgrp.close()
    
    # array to store delta concentration spatially in reduction area
    DC_target_arrray = zeros((n_lat, n_lon))

    # make a window
    window = create_window(inner_radius)
    (n_lon_inner_win, n_lat_inner_win) = window.shape
            
    # loop over all nuts in 
    for nuts_id in range(n_nuts):
        # initialize delta_conc
        nuts_code = nuts_codes[nuts_id]
        if NOx == True:
            delta_conc_no2[nuts_code] = 0
            delta_conc_nox[nuts_code] = 0
        else:
            delta_conc[nuts_code] = 0
            
        # print the progress
        progress = float(nuts_id) / float(n_nuts) * 100
        sys.stdout.write('\r')
        sys.stdout.flush()
        sys.stdout.write('progress:%f\r' % progress)
        sys.stdout.flush()
    
        reduction_area_array = rootgrp_nuts.variables['AREA'][nuts_id,:,:] / 100.0
        
        # calculate the delta emissions, dictionary per pollutant a matrix of dimension n_lat x n_lon
        delta_emission_dict = create_delta_emission(path_emission_cdf, precursor_lst, reduction_area_array, path_reduction_txt)
       
        pad_delta_emission_dict = {}
        for precursor in precursor_lst:
            pad_delta_emission_dict[precursor] = lib.pad(delta_emission_dict[precursor], inner_radius, 'constant', constant_values=0)
        
        # apply source receptor relationships
        # -----------------------------------
        
        for precursor in precursor_lst:
            # get the model coefficients
            alpha_ij = alpha_dict[precursor][i_lat_target, i_lon_target]
            omega_ij = omega_dict[precursor][i_lat_target, i_lon_target]
            # if the model exists, use it
            if not(isnan(alpha_ij)):
                # select the emissions inside the weighting window for the precursor
                emissions_window = pad_delta_emission_dict[precursor][i_lat_target:(i_lat_target + n_lon_inner_win), i_lon_target:(i_lon_target + n_lat_inner_win)]
                # apply the weighting factors and sum over the whole window
                weighted_emissions_window = ((power(window, omega_ij)) * emissions_window).sum()
                # for NO2 the SR model gives delta NOx concentration
                if NOx == True:
                    delta_conc_nox[nuts_code] = delta_conc_nox[nuts_code] + alpha_ij * weighted_emissions_window
                # for PM2.5 and PM10 the delta concentration of the respective PM
                else:
                    delta_conc[nuts_code] = delta_conc[nuts_code] + alpha_ij * weighted_emissions_window
                
        # In the case of NOx, the NO2 concentration has to be calculated with the NO2 fraction correlation
        if NOx == True:
            delta_conc_no2[nuts_code] = deltaNOx_to_deltaNO2(delta_conc_nox[nuts_code], target_conc_basecase_nox, target_conc_basecase_no2)
            # there is a choice to be made, report delta NO2 or NOx in the spatial result
            # NOx
            DC_target_arrray = DC_target_arrray + delta_conc_nox[nuts_code] * reduction_area_array
            # or NO2
            # DC_target_arrray = DC_target_arrray + delta_conc_no2[nuts_code] * reduction_area_array
        else:
            # create an output map with in each nuts the DC in the target cell
            DC_target_arrray = DC_target_arrray + delta_conc[nuts_code] * reduction_area_array
        
    # close nuts cdf
    rootgrp_nuts.close()
    
    # sort nuts codes from delta_conc from high to low delta conc
    if NOx == True:
        sorted_nuts_codes_nox = sorted(delta_conc_nox, key=lambda i: delta_conc_nox[i], reverse=True)
        sorted_nuts_codes_no2 = sorted(delta_conc_no2, key=lambda i: delta_conc_no2[i], reverse=True)
    else:
        sorted_nuts_codes = sorted(delta_conc, key=lambda i: delta_conc[i], reverse=True) 
    
    # write the result to a netcdf file
    path_DC_target_cdf = path_result_cdf + 'radius_result.nc'
    rootgrp = Dataset(path_DC_target_cdf, 'w', format = 'NETCDF3_CLASSIC')
    rootgrp.createDimension('latitude', n_lat)
    rootgrp.createDimension('longitude', n_lon)
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    latitudes.units = "degrees_north"
    latitudes[:] = latitude_array
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    longitudes.units = "degrees_east"
    longitudes[:] = longitude_array
    area = rootgrp.createVariable('AREA', 'f4', ('latitude', 'longitude',))
    area[:] = DC_target_arrray
    rootgrp.close()
    
    # write a result file
    f_res = open(path_result_cdf + 'radius_result.txt', 'w')
    f_res.write('nuts_code\t%\n')
    if NOx == True:
        # there is a choice to be made, report NO2 or NOx relative potential
        # NOx
        for nuts_code in sorted_nuts_codes_nox:
            f_res.write('%s\t%e\n' % (nuts_code, delta_conc_nox[nuts_code] / target_conc_basecase_nox / (alpha_potency / 100) * 100)) # rel potential in percentage
        # or NO2
#         for nuts_code in sorted_nuts_codes_no2:
            # f_res.write('%s\t%e\n' % (nuts_code, delta_conc_no2[nuts_code] / target_conc_basecase_no2 / (alpha_potency / 100) * 100)) # rel potential in percentage
    # for PM2.5 or PM10
    else:
        for nuts_code in sorted_nuts_codes:
            f_res.write('%s\t%e\n' % (nuts_code, delta_conc[nuts_code] / target_conc_basecase / (alpha_potency / 100) * 100)) # rel potential in percentage
    f_res.close()
    
    
if __name__ == '__main__':
    
    # the testing is moved to test_module6.py    
     
    pass




