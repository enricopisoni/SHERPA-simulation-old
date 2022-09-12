import sys
from os import getcwd
from time import time

if __name__ == '__main__':
    if (getcwd().endswith('test')):
        sys.path.append("..")
    else:
        sys.path.append(".")
    from modules import module4
    from globals import path_result_cdf_test, downscale_request
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