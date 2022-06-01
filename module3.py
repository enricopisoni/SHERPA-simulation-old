'''
Created on Aug 20, 2015

Module 3 executes module 4 for for:
- all macrosectros individually and all together for 1 or more precursor => module 3a
- all precursors individually and all together for 1 or more macrosectors => module 3b

The calculated potentials are:
- DC / alpha (absolute potential)
- DC / C / alfa (relative potential)

@author: degraba
'''
from netCDF4 import Dataset
from module4 import module4
from sherpa_auxiliaries import create_emission_reduction_dict, write_progress_log
from sherpa_globals import alpha_potency, sector_lst, \
    path_emission_cdf_test, path_area_cdf_test, path_reduction_mod3a1P_txt_test, path_reduction_mod3a2P_txt_test, \
    path_reduction_mod3b_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test
from time import time
from numpy import zeros
from os import remove

#--------------------------
# MODULE 3a
#--------------------------

# potency for all macrosectros individually and all together for one or more precursors
def module3a(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf, downscale_request):
          
    # get precursor list from model
    rootgrp = Dataset(path_model_cdf, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    rootgrp.close()

    # create a dictionary with reductions per precursor and macro sector
    emission_reduction_dict = create_emission_reduction_dict(path_reduction_txt)

# is this necessary???????????    
#     # check which precursor reductions are non-zero
#     for precursor in emission_reduction_dict.keys():
#         precursor_reduction = 0
#         for 
    
    # look up which precursor(s) is(are) reduced
    reduced_precursor_lst = []
    for precursor in emission_reduction_dict.keys():
        sum_over_snaps = 0
        for snap in emission_reduction_dict[precursor].keys():
            sum_over_snaps += emission_reduction_dict[precursor][snap]
        if sum_over_snaps > 0:
            reduced_precursor_lst.append(precursor)
    # print(reduced_precursor_lst)
            
    # make module 4 reduction input files for each snap sector
    f_red_mod_3 = open(path_reduction_txt, 'r')
    header = f_red_mod_3.readline()
    f_red_mod_3.close()
    
    # declare a dictonary to store results
    results = {}
    
    # progress counter of module 3a
    counter = 0.0
    
    # name of the file in which the progress of the calculation will be kept and transfered to sub processes
    progress_log_filename = path_result_cdf + 'proglogmod3.txt'
     
    for snap in sector_lst[0:-1]:
        
        # write progress log file
        start = float(counter) / (len(sector_lst[0:-1]) + 1) * 100
        divisor = len(sector_lst[0:-1]) + 1
        write_progress_log(progress_log_filename, start, divisor)
        
        # print(snap)
        # create the emission reduction file for the precusor and only one sector
        filename_mod4_reductions = path_result_cdf + 'mod4_reductions_snap_%s.txt' % (snap)
        f_red_mod_4_snap = open(filename_mod4_reductions, 'w')
        f_red_mod_4_snap.write(header)
        for precursor in precursor_lst:
            f_red_mod_4_snap.write(precursor)
            for snap2 in sector_lst[0:-1]:
                if precursor in reduced_precursor_lst and snap2 == snap:
                    f_red_mod_4_snap.write('\t' + str(alpha_potency))
                else:
                    f_red_mod_4_snap.write('\t0')
            f_red_mod_4_snap.write('\n')
        f_red_mod_4_snap.close()    

        # call module 4 with the newly created emission reduction file
        res_mod4_snap = module4(path_emission_cdf, path_area_cdf, filename_mod4_reductions, path_base_conc_cdf, path_model_cdf, path_result_cdf, 
                                downscale_request, progress_log_filename)

        # remove potencies output
        remove(path_result_cdf + 'potencies.nc')
        
        # update counter
        counter += 1
        # print(counter)
                
        # store the results for each individual snap sector
        results[snap] = res_mod4_snap
        
        # remove filename_mod4_reductions
        remove(filename_mod4_reductions)
   
    # write progress log file
    start = float(counter) / (len(sector_lst[0:-1]) + 1) * 100
    divisor = len(sector_lst[0:-1]) + 1
    write_progress_log(progress_log_filename, start, divisor)
            
    # execute module 4 with a reduction in all sectors together
    res_mod4_all = module4(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf, 
                           downscale_request, progress_log_filename)
    n_lat = res_mod4_all['n_lat'] 
    n_lon = res_mod4_all['n_lon']  
    n_nuts = len(sector_lst[0:-1])

    # remove potencies output
    remove(path_result_cdf + 'potencies.nc')
    
    # remove progress log file
    remove(progress_log_filename)

    # create results netcdf
    # -----------------------
    filename_result_cdf = path_result_cdf + 'potencies_overview_per_sector.nc'
    rootgrp = Dataset(filename_result_cdf, 'w', format='NETCDF3_CLASSIC')
    
    # create dimensions in the netcdf file
    rootgrp.createDimension('GNFRsector', n_nuts)
    rootgrp.createDimension('latitude', n_lat)
    rootgrp.createDimension('longitude', n_lon)
    GNFRsectors = rootgrp.createVariable('GNFRsector', 'f4', ('GNFRsector',))
    GNFRsectors[:] = sector_lst[0:-1]
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    latitudes.units = "degrees_north"
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    longitudes.units = "degrees_east"
    latitudes[:] = res_mod4_all['latitude_array']
    longitudes[:] = res_mod4_all['longitude_array']
    
    # create variables and initialize with zeros
    DC_alpha_snap_var = rootgrp.createVariable('DC_alpha_snap', 'f4', ('GNFRsector', 'latitude', 'longitude',))
    DC_alpha_snap_var.units = "ug/m3"
    DC_alpha_snap_var[:] = zeros((n_nuts, n_lat, n_lon))
    DC_C_alpha_snap_var = rootgrp.createVariable('DC_C_alpha_snap', 'f4', ('GNFRsector', 'latitude', 'longitude',))
    DC_C_alpha_snap_var.units = "%"
    DC_C_alpha_snap_var[:] = zeros((n_nuts, n_lat, n_lon))
#     DC_DE_snap_var = rootgrp.createVariable('DC_DE_snap', 'f4', ('GNFRsector', 'latitude', 'longitude',))
#     DC_DE_snap_var[:] = zeros((n_nuts, n_lat, n_lon))

    DC_alpha_all_var = rootgrp.createVariable('DC_alpha_all', 'f4', ('latitude', 'longitude',))
    DC_alpha_all_var.units = "ug/m3"
    DC_alpha_all_var[:] = res_mod4_all['DC_alpha']
    DC_C_alpha_all_var = rootgrp.createVariable('DC_C_alpha_all', 'f4', ('latitude', 'longitude',))
    DC_C_alpha_all_var.units = "%"
    DC_C_alpha_all_var[:] = res_mod4_all['DC_C_alpha'] * 100
    
    for snap in sector_lst[0:-1]:
        DC_alpha_snap_var[snap - 1, :, :] = results[snap]['DC_alpha']
        DC_C_alpha_snap_var[snap - 1, :, :] = results[snap]['DC_C_alpha'] * 100
    rootgrp.close()
    
#--------------------------
# MODULE 3b
#--------------------------

def module3b(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf, downscale_request):

    # get precursor list from model
    rootgrp = Dataset(path_model_cdf, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    rootgrp.close()

    # create a dictionary with reductions per precursor and macro sector
    emission_reduction_dict = create_emission_reduction_dict(path_reduction_txt)
    
    # look up which sector(s) that is(are) reduced
    reduced_sectors = []
    for sector in sector_lst[0:-1]:
        sum_over_precursors = 0
        for precursor in precursor_lst:
            sum_over_precursors += emission_reduction_dict[precursor][sector]
        if sum_over_precursors > 0:
            reduced_sectors.append(sector)
            
    # make module 4 reduction input files for each snap sector
    f_red_mod_3 = open(path_reduction_txt, 'r')
    header = f_red_mod_3.readline()
    f_red_mod_3.close()
    
    # declare a dictonary to store results
    results = {}

    # progress counter of module 3b
    counter = 0.0
    
    # name of the file in which the progress of the calculation will be kept and transfered to sub processes
    progress_log_filename = path_result_cdf + 'proglogmod3.txt'
    
    for precursor in precursor_lst:
        
        # write progress log file
        start = float(counter) / (len(precursor_lst) + 1) * 100
        divisor = len(precursor_lst) + 1
        write_progress_log(progress_log_filename, start, divisor)
        
        # create the emission reduction file for the precusor and only one sector
        filename_mod4_reductions = path_result_cdf + 'mod4_reductions_precursor_of_%s.txt' % (precursor)
        f_red_mod_4_snap = open(filename_mod4_reductions, 'w')
        f_red_mod_4_snap.write(header)
        for precursor2 in precursor_lst:
            f_red_mod_4_snap.write(precursor2)
            for sector in sector_lst:
                if sector in reduced_sectors and precursor2 == precursor:
                    f_red_mod_4_snap.write('\t' + str(alpha_potency))
                else:
                    f_red_mod_4_snap.write('\t0')
            f_red_mod_4_snap.write('\n')
        f_red_mod_4_snap.close()    

        # call module 4 with the newly created emission reduction file
        res_mod4_snap = module4(path_emission_cdf, path_area_cdf, filename_mod4_reductions, path_base_conc_cdf, path_model_cdf, path_result_cdf, 
                                downscale_request, progress_log_filename)
        
        # remove potencies output
        remove(path_result_cdf + 'potencies.nc')
        
        # update counter
        counter += 1
        # print(counter)

        # store the results for each individual precursor
        results[precursor] = res_mod4_snap
        
        # delete filename_mod4_reductions
        remove(filename_mod4_reductions)
        
    # write progress log file
    start = float(counter) / (len(precursor_lst) + 1) * 100
    divisor = len(precursor_lst) + 1
    write_progress_log(progress_log_filename, start, divisor)
     
    # execute module 4 with a reduction in all precursors
    res_mod4_all = module4(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf, 
                           downscale_request, progress_log_filename)
    n_lat = res_mod4_all['n_lat']  
    n_lon = res_mod4_all['n_lon']  
    
    # remove potencies output
    remove(path_result_cdf + 'potencies.nc')

    # remove progress log file
    remove(progress_log_filename)
    
    # create results netcdf
    # -----------------------
    filename_result_cdf = path_result_cdf + 'potencies_overview_per_precursor.nc'
    rootgrp = Dataset(filename_result_cdf, 'w', format='NETCDF3_CLASSIC')
    
    # create dimensions in the netcdf file
    rootgrp.createDimension('latitude', n_lat)
    rootgrp.createDimension('longitude', n_lon)
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    latitudes.units = "degrees_north"
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    longitudes.units = "degrees_east"
    latitudes[:] = res_mod4_all['latitude_array']
    longitudes[:] = res_mod4_all['longitude_array']
    
    # create variables and initialize with zeros
    for precursor in precursor_lst:
        DC_alpha_var = rootgrp.createVariable('DC_alpha_precursor_%s' % (precursor) , 'f4', ('latitude', 'longitude',))
        DC_alpha_var.units = "ug/m3"
        DC_alpha_var[:] = zeros((n_lat, n_lon))
        DC_C_alpha_var = rootgrp.createVariable('DC_C_alpha_precursor_%s' % (precursor), 'f4', ('latitude', 'longitude',))
        DC_C_alpha_var.units = "%"
        DC_C_alpha_var[:] = zeros((n_lat, n_lon))

        DC_alpha_var[:, :] = results[precursor]['DC_alpha']
        DC_C_alpha_var[:, :] = results[precursor]['DC_C_alpha'] * 100
    
    DC_alpha_all_var = rootgrp.createVariable('DC_alpha_all', 'f4', ('latitude', 'longitude',))
    DC_alpha_all_var.units = "ug/m3"
    DC_alpha_all_var[:] = res_mod4_all['DC_alpha']
    DC_C_alpha_all_var = rootgrp.createVariable('DC_C_alpha_all', 'f4', ('latitude', 'longitude',))
    DC_C_alpha_all_var.units = "%"
    DC_C_alpha_all_var[:] = res_mod4_all['DC_C_alpha'] * 100
    
    rootgrp.close()


# test function
if __name__ == '__main__':
    
    # test module 3a for one precursor
    start = time()
    module3a(path_emission_cdf_test, path_area_cdf_test, path_reduction_mod3a1P_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test,
             downscale_request)
    stop = time()
    print('Module 3a1P calculation time = %f seconds' % (stop - start))

    # test module 3a for multiple precursors
    start = time()
    module3a(path_emission_cdf_test, path_area_cdf_test, path_reduction_mod3a2P_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test,
             downscale_request)
    stop = time()
    print('Module 3a1P calculation time = %f seconds' % (stop - start))

    # test module 3b 
    start = time()
    module3b(path_emission_cdf_test, path_area_cdf_test, path_reduction_mod3b_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test, 
             downscale_request)
    stop = time()
    print('Module 3a1P calculation time = %f seconds' % (stop - start))
    
    # test module 3b for multiple snap sectors
#     fullmodel = 'input/fullFunction/SR_PM25_Y_fullFunction.nc'
#     start = time()
#     module3b(path_emission_cdf_test, path_area_cdf_test, path_reduction_mod3b_txt_test, path_base_conc_cdf_test, fullmodel, path_result_cdf_test)
#     stop = time()
#     print('Module 3b calculation time = %f' % (stop - start))

    pass