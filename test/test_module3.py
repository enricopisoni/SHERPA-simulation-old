import sys
from time import time
from os import getcwd

if __name__ == '__main__':
    if (getcwd().endswith('test')):
        sys.path.append("..")
    else:
        sys.path.append(".")
    from modules import module3a, module3b
    from globals import path_emission_cdf_test, path_area_cdf_test, path_reduction_mod3a1P_txt_test, path_reduction_mod3a2P_txt_test, \
        path_reduction_mod3b_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test, downscale_request
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