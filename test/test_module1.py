'''
Created on 5 Apr 2017

This module runs some tests for module 1

@author: degraba
'''
from time import time
from os import rename, getcwd
import sys

if __name__ == '__main__':


    # Check the current working directory in order to import the main folder of the project
    # This will allow to correctly import the modules folder, otherwise the modules folder would not be visible
    if (getcwd().endswith('test')):
        sys.path.append("..")
    else:
        sys.path.append(".")
    from modules import module1
    # lastest model on 2017/04/04: O:/Integrated_assessment/SHERPA/20170322_v18_SrrResults_PotencyBased/
    model_path = 'O:/Integrated_assessment/SHERPA/20170322_v18_SrrResults_PotencyBased/'
    emission_folder = model_path + '1_base_emissions/'
    concentrations_folder = model_path + '2_base_concentrations/'
    model_folder = model_path + '3_source_receptors/'
    
    pollutant = 'NO2'
    path_emission_cdf = emission_folder + 'BC_emi_' + pollutant + '_Y.nc'
    reduction_area = 'input/London_region.nc'
    # reduction_area = 'input/EMI_RED_ATLAS_NUTS1/UKI.nc'
    reduction_snap = 'input/user_reduction_snap7.txt'
    if pollutant == 'NO2':
        path_base_conc_cdf = concentrations_folder + 'BC_conc_NO2_NO2eq_Y_mgm3.nc'
        path_model_cdf = model_folder + 'SR_NO2eq_Y_20170322_potencyBased.nc'
    else:
        path_base_conc_cdf = concentrations_folder + 'BC_conc_' + pollutant + '.nc'
        path_model_cdf = model_folder + 'SR_' + pollutant + '_Y_20170322_potencyBased.nc'
        
    output_path = 'output/test/'
    tag = pollutant + '_London_region'
 
    # run module 1 with progress log
    proglog_filename = output_path + 'proglog'
    # write_progress_log(proglog_filename, 25, 2)
    start = time()
    module1(path_emission_cdf, reduction_area, reduction_snap, path_base_conc_cdf, path_model_cdf, output_path) 
    
    stop = time()
    print('Module 1 run time: %s sec.' % (stop-start))
    
    # rename output
    rename(output_path + 'delta_concentration.nc', output_path + 'delta_concentration_' + tag + '.nc')
    
    
    pass