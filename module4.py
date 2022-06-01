'''
Created on Jul 16, 2015

Calculates potencies of:
    1 precursor
    1 macrosector
    in a selection of nuts 
for all nuts
The emission reduction to calculate potencies is 50% (alfa)

The calulated potencies are:
- Delta C / Delta E
- (Delta C / C) / (Delta E / E) = delta_C / (C*alfa)
- delta_C / alfa


@author: degraba
'''
# imports
from netCDF4 import Dataset
from module1 import module1
from time import time
from sherpa_globals import alpha_potency, path_result_cdf_test, path_reduction50all_txt_test
# path_emission_cdf_test, path_base_conc_cdf_test, path_model_cdf_test   
from sherpa_auxiliaries import create_emission_reduction_dict, write_progress_log, read_progress_log
from numpy import zeros
from os import remove

# function definition of source receptor model
def module4(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf, downscale_request, *progresslog):

    # get precursor list from model
    rootgrp = Dataset(path_model_cdf, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    rootgrp.close()
    
    # module 4 can be ran independently without progress log argument or...
    if progresslog:
        progresslog_filename = progresslog[0]
        progress_dict = read_progress_log(progresslog_filename)
    # called by another model which passes a progress log argument
    else:
        progresslog_filename = path_result_cdf + 'progresslog'
        write_progress_log(progresslog_filename, 0, 1)
        progress_dict = read_progress_log(progresslog_filename)
    
    # calculate the delta emissions, dictionary per pollutant a matrix of dimension n_lat x n_lon
    mod1_res = module1(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf, downscale_request, 
                       progresslog_filename)
    delta_conc = mod1_res['delta_conc']
    delta_emis_dict = mod1_res['delta_emis_dict']
    n_lat = mod1_res['n_lat']
    n_lon = mod1_res['n_lon']
    
    # remove progress log if he was made inside module 4, not if it was an external argument
    if not progresslog:
        remove(progresslog_filename)

    # only if only one precursor is reduced delta_emission is used for potency calculations
    # check this
    emission_reduction_dict = create_emission_reduction_dict(path_reduction_txt)
    number_reduced_precursors = 0
    for precursor in precursor_lst:
        sum_reductions = 0
        for snap in emission_reduction_dict[precursor].keys():
            sum_reductions += emission_reduction_dict[precursor][snap]
        if sum_reductions > 0:
            number_reduced_precursors += 1
            reduced_precursor = precursor
            delta_emis = delta_emis_dict[precursor]
            # print('precusor %s was reduced' % precursor)
    if number_reduced_precursors != 1:
        reduced_precursor = ''
    # print(sum(delta_emis))
    
    # read baseline concentrations
    rootgrp = Dataset(path_base_conc_cdf, 'r')
    conc = rootgrp.variables['conc'][:]

    #ENR 20180126 - fix for NO2. in case of NO2, load NO2 and not NOx variable as baseline
    if (path_model_cdf.find('NO2eq') > -1):
        conc = rootgrp.variables['NO2'][:]

    # close model netcdf
    rootgrp.close()
    
    # calculate different potencies
    #------------------------------
#     print('check delta conc')
#     print(sum(delta_conc))
    # Delta_C/alfa
    DC_alpha = delta_conc / (alpha_potency / 100.0)
    
    DC_C_alpha = delta_conc / conc / (alpha_potency / 100.0)
       
    if number_reduced_precursors == 1:
        ## select reduced precursor automatically !!!!!!!!!!!
        DC_DE = zeros((n_lat, n_lon))
        DC_DE = delta_conc / delta_emis_dict[reduced_precursor].sum() * 1000
#         for i in range(n_lat):
#             for j in range(n_lon):
#                 if delta_emis_dict[reduced_precursor][i, j] > 0:
                    # DC_DE[i, j] = delta_conc[i, j] / delta_emis_dict[reduced_precursor][i, j]
#                     DC_DE[i, j] = delta_conc[i, j] / DE

    # create a result netcdf 
    # -----------------------
    if progress_dict['start'] != -1:
        filename_result_cdf = path_result_cdf + 'potencies.nc'
        rootgrp = Dataset(filename_result_cdf, 'w', format='NETCDF3_CLASSIC')
        
        # create dimensions in the netcdf file
        rootgrp.createDimension('latitude', n_lat)
        rootgrp.createDimension('longitude', n_lon)
        latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
        latitudes.units = "degrees_north"
        longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
        longitudes.units = "degrees_east"
        latitudes[:] = mod1_res['latitude_array']
        longitudes[:] = mod1_res['longitude_array']
            
        # create potency variables data
        DC_alpha_var = rootgrp.createVariable('DC_alpha', 'f4', ('latitude', 'longitude',))
        DC_alpha_var[:] = DC_alpha
        DC_C_alpha_var = rootgrp.createVariable('DC_C_alpha', 'f4', ('latitude', 'longitude',))
        DC_C_alpha_var[:] = DC_C_alpha
        if number_reduced_precursors == 1:
            DC_DE_var = rootgrp.createVariable('DC_DE', 'f4', ('latitude', 'longitude',))
            DC_DE_var[:] = DC_DE
              
        delta_conc_var = rootgrp.createVariable('delta_conc_var', 'f4', ('latitude', 'longitude',))
        delta_conc_var[:] = delta_conc
        delta_emis_var = rootgrp.createVariable('delta_emis_var', 'f4', ('latitude', 'longitude',))
        delta_emis_var[:] = delta_emis
        
        rootgrp.close()
        
    # create a results object
    mod4_res = {}
    mod4_res['n_lat'] = n_lat
    mod4_res['n_lon'] = n_lon
    mod4_res['latitude_array'] = mod1_res['latitude_array']
    mod4_res['longitude_array'] = mod1_res['longitude_array']
    # add DC_alpha to the result dictionary
    mod4_res['DC_alpha'] = DC_alpha
    # add DC_C_alpha to the result dictionary
    mod4_res['DC_C_alpha'] = DC_C_alpha
    if number_reduced_precursors == 1:
        # add DC_DE to the result dictionary
        mod4_res['DC_DE'] = DC_DE
        
    return mod4_res

# test module 4    
if __name__ == '__main__':

    emission_1611_test = 'input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
    concentration_1611_test = 'input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
    model_1611_test = 'input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
    reduction_file = 'input/user_reduction_PM25.txt'

    # run module 4 without progress log
    start = time()
    module4(emission_1611_test, 'input/London_region.nc', reduction_file, concentration_1611_test, model_1611_test, path_result_cdf_test, downscale_request)
    stop = time()
    print('Module 4 calculation time = %f' % (stop - start))
    
    # run module 1 with progress log
#     proglog_filename = path_result_cdf_test + 'proglog_mod4'
#     write_progress_log(proglog_filename, 25, 2)
#     start = time()
#     module4(emission_1611_test, 'input/London_region.nc', path_reduction50all_txt_test, concentration_1611_test, model_1611_test, path_result_cdf_test, proglog_filename)
#     stop = time()
#     print('Module 4 run time: %s sec.' % (stop-start))
#     remove(proglog_filename)
      
    pass
    
    
    
