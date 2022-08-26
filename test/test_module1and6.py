'''
Created on 6 Apr 2017

this script runs some tests on module 6

@author: degraba
'''

from time import time
from os import rename, remove, getcwd
import sys
import re

if __name__ == '__main__':

    # Check the current working directory in order to import the main folder of the project
    # This will allow to correctly import the modules folder, otherwise the modules folder would not be visible
    if (getcwd().endswith('test')):
        sys.path.append("..")
    else:
        sys.path.append(".")
    from modules import module1, module6
    
    # run module 6
    # lastest model on 2017/04/04: O:/Integrated_assessment/SHERPA/20170322_v18_SrrResults_PotencyBased/
    model_path = 'O:/Integrated_assessment/SHERPA/20170322_v18_SrrResults_PotencyBased/'
    emission_folder = model_path + '1_base_emissions/'
    concentrations_folder = model_path + '2_base_concentrations/'
    model_folder = model_path + '3_source_receptors/'
    nuts_level = '1'
    
    for pollutant in ['PM25']: # ['PM10', 'PM25', 'NO2']:
        path_emission_cdf = emission_folder + 'BC_emi_' + pollutant + '_Y.nc'
        reduction_multi_area = 'input/EMI_RED_ATLAS_NUTS' + nuts_level + '.nc'
        reduction_snap = 'input/user_reduction_all50.txt'
        if pollutant == 'NO2':
            path_base_conc_cdf = concentrations_folder + 'BC_conc_NO2_NO2eq_Y_mgm3.nc'
            path_model_cdf = model_folder + 'SR_NO2eq_Y_20170322_potencyBased.nc'
        else:
            path_base_conc_cdf = concentrations_folder + 'BC_conc_' + pollutant + '_Y.nc'
            path_model_cdf = model_folder + 'SR_' + pollutant + '_Y_20170322_potencyBased.nc' #SR_PM10_Y_20170322_potencyBased
            
        output_path = 'output/test/'
        # London coordinates
        tag = pollutant + '_London'
        target_cell_lat = 51.51    
        target_cell_lon = -0.13
         
        # run module 1 with progress log
        start = time()
        # module6(path_emission_cdf, reduction_multi_area, target_cell_lat, target_cell_lon, reduction_snap, path_base_conc_cdf, path_model_cdf, output_path)
        # print(DC)
        stop = time()
        print('Module 6 run time: %s sec.' % (stop-start))
        
        # rename output
#         rename(output_path + 'radius_result.nc', output_path + 'radius_result_' + tag + '.nc')
        new_radius_result_name = output_path + 'radius_result_' + tag + '.txt'
#         rename(output_path + 'radius_result.txt', new_radius_result_name)
        
        # consistency check with module1
        # ------------------------------
        
        # read module6 output
        f_res = open(new_radius_result_name)
        radius_result_dict = {'nuts_codes': [], 'potential': []}
        line = f_res.readline().rstrip()        # header
        line = f_res.readline().rstrip()        # first line
        while len(line) > 0:
            [nuts_code, potential] = re.split('\t', line)
            radius_result_dict['nuts_codes'].append(nuts_code.rstrip())
            radius_result_dict['potential'].append(float(potential))
            line = f_res.readline().rstrip()
        f_res.close()
        
        # run module 1 for the 3 most influential nuts
        for i_nuts in range(0, 3):
            nuts_code = radius_result_dict['nuts_codes'][i_nuts]
            print(nuts_code)
            reduction_area = 'D:/workspace/sherpa/trunk/input/EMI_RED_ATLAS_NUTS' + nuts_level + '/' + nuts_code + '.nc'
            start = time()
            module1(path_emission_cdf, reduction_area, reduction_snap, path_base_conc_cdf, path_model_cdf, output_path)
            stop = time()
            print('Module 6 run time: %s sec.' % (stop-start))
            # rename output delta_concentration
            rename(output_path + 'delta_concentration.nc', output_path + 'delta_concentration_' + tag + '_from_' + nuts_code + '.nc')
            # delete delta emissions
            remove(output_path + 'delta_emission.nc')
    
    pass