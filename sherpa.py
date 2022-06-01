# 20220530

# for downscaling of NO2
# - run module 1 for NO2
# - run module from Bertrand to downscale to 0.01 deg

# for downscaling of PM2.5 or PM10:
# - run module 1 for PM with downscale_request=0 (reduce all precursors) - get total delta PM
# - run module 1 for PM with downscale_request=1 (as before, but reduce only PPM) - deta delta PPM
# - run module from Bertrand to downscale to 0.01 deg (please note that the module uses PM and PPM for its computation...only PPM is downscaled)

#note that now automatic names are created, so that results are not overwritten

###############################################################################

'''
Created on Jun 23, 2015

simulates the main tool developped by teraria and call the source-receptor model

The sherpa tool can be ran in different modes

1) Module 1: scen_nuts
----------------------
In this module the concentrations are calculated for a given emission reductions scenario
Inputs: - baseline emissions (per precursor and cell),
        - a netcdf defining the area where emission reduction will be applied
        - emission reductions (per precursor and macrosector)
        - model coefficients (per pollutant, precursor, cell)
        - path where the results will be written
        - eventually a progress log to be able to report on the progress of the whole calculation when module 1
        is called by another moduel
        
output: - netcdf with concentration changes per pollutant and cell
        - delta emission netcdf with emission changes per precursor and cell
       
The calculation is optimized using a flat weight over the whole domain. This allows to update only the scale
factor of this flat weight. The bell shape central weighting factors have to be recalculated for each cell.

2) Module 2: Scen_Atlas
-----------------------
This module calculates the concentration reduction in each cell of a nuts sector due to emission reductions in that same nuts sector.
Input: - baseline emissions (per precursor and cell)
       - a netcdf defining the areas where emission reductions will be applied
       - txt file with reductions (per snap and precursor) to be applied in each nuts
       - model parameters
       - result path
Output: - netcdf with concentration reductions, one layer per NUTS

3) Module 3: source apportionment
---------------------------------
Module 3 executes module 4 for for:
- all macrosectros individually and all together for 1 or more precursor => module 3a
- all precursors individually and all together for 1 or more macrosectors => module 3b

The calculated potencies are:
- DC / C
- DC / C / alfa

4) Module 4: potencies
-----------------------
Calculate potencies for an emission reduction in a given area. 

5) Module 5: potency atlas
--------------------------

6) Module 6: radius calculation
-------------------------------

7) Module 7: source a. 
    
-------------------------------

8) Moulde 8: health impact assessment

9) Module 9: aggregation - dipslay aggregated values in graphical interface
    
A target cell is selected. For this cell the concentration changes is calculated due to a defined emission reduction
applied in each NUTS individually
Input: -

Output: - text file with NUTS code and corresponding concentration change in the target cell.
        - netcdf with each NUTS with the vallue of its concentration change caused in the target cell.


@author: degraba
'''

from module1 import module1
#from module2 import module2
from module3 import module3a, module3b
from module4 import module4
#from module5 import module5
from module6 import module6
#from module7_SA import module7
from module8_healthia import module8_healthia
from module9_aggregation import  module9_aggregation
#20211020 - EP
from module10 import module10
from sherpa_globals import path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, \
    path_model_cdf_test, path_result_cdf_test, path_nuts0_cdf_test, path_nuts1_cdf_test, path_nuts2_cdf_test, \
    path_base_conc_cdf_test, path_healthbl_test, path_config_json_test, \
    fua_intersect_dir, dbf_dir, target_list, nuts_intersect_dir, path_natural_dir_test, aggr_zones, \
    path_logo_test, aggrinp_txt, sector_lst, downscale_request
from sherpa_auxiliaries import is_number
from sys import argv
import os.path
from time import time

#20220527, last commigg before starting the downscaling work


if __name__ == '__main__':
    
    if len(argv) == 1:
        # no command arguments are provided, sherpa is ran in test mode with fixed input arguments
        
        # run module 1 with test inputs
        start = time()
        module1(path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test,
                downscale_request) #downscale_request=0 nothing change, downscale_request=1 you only reduce PPM emissions
        stop = time()
        print('Module 1 run time: %s sec.' % (stop-start))

###############################################################################        
        #20200318 - module not used anymore
        # run module 2 with test inputs
#        start = time()
#        module2(path_emission_cdf_test, path_nuts0_cdf_test, path_reduction_txt_test, path_model_cdf_test, path_result_cdf_test)
#        stop = time()
#        print('Module 2 run time: %s sec.' % (stop-start))
###############################################################################
        
        # run module 3a test inputs with test inputs
        start = time()
        # module3a(path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test,
                    # downscale_request)
        stop = time()
        print('Module 3a calculation time = %f' % (stop - start))

        # run module 3b with test inputs
        start = time()
        # module3b(path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test,
                    # downscale_request)
        stop = time()
        print('Module 3b calculation time = %f' % (stop - start))

        # run module 4  with test inputs
        start = time()
        # module4(path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test, 
                # downscale_request)
        stop = time()
        print('Module 4 run time: %s sec.' % (stop-start))
        
        
###############################################################################        
        #20200318 - module not used anymore        
#        # run module 5 test inputs
#        start = time()
#        module5(path_emission_cdf_test, path_nuts0_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test)
#        stop = time()
#        print('Module 5 calculation time = %f' % (stop - start))
###############################################################################
        
        # run module 6 test inputs
        start = time()
        # paris
        target_cell_lat = 48.85     # 51.51
        target_cell_lon = 2.35  # 9.19      #-0.13
        # module6(path_emission_cdf_test, path_nuts1_cdf_test, target_cell_lat, target_cell_lon, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test)
        stop = time()
        print('Module 6 calculation time = %f' % (stop - start))

###############################################################################        
        #20220531 - not used anymore
        # # run module 7 test inputs    
        # start = time()
        # # module7(path_emission_cdf_test, path_base_conc_cdf_test, path_natural_dir_test,
        # #         path_model_cdf_test, fua_intersect_dir, nuts_intersect_dir, dbf_dir, 
        # #         target_list, path_result_cdf_test, path_logo_test, aggr_zones, 'PM25')
        # stop = time()
        # print('Module 7 calculation time = %f' % (stop - start))
###############################################################################

        # run module 8 test inputs        
        start = time()
        # module8_healthia(path_healthbl_test, path_result_cdf_test,
        #                 path_config_json_test, path_base_conc_cdf_test)
        stop = time()
        print('Module 8 calculation time = %f' % (stop - start))
        
        # run module 9 test inputs        
        start = time()
        # module9_aggregation('../input_EDGAR_EMEP/config/aggrinp.json')        
        stop = time()
        print('Module 9 calculation time = %f' % (stop - start))
        
        #EP 20211020
        # run module 10 test inputs        
        start = time()
        #Input file from Alessandro
        path_input_txt='./inputBB/costs_input.csv'
        #path_input_txt='./inputBB/costs_input_nognfr13.csv'
        #Static database of costs
        path_datab_txt='./inputBB/SHERPA_unitc_DB_CAO2_2015.csv'
        #Output file to be managed by Alessandro 
        path_costd_txt='./outputBB/costdetails.txt'
        # module10(path_input_txt,path_datab_txt,path_costd_txt)
        
        stop = time()
        print('Module 10 calculation time = %f' % (stop - start))
    
    else:
        # check which module has to be ran
        module = int(argv[1])

        # check if all files (emissions, emission reduction, model) exist
        for i_input in range(2, len(argv)):    
            if not(os.path.exists(argv[i_input])) and not(is_number(argv[i_input])):
                    print('WARNING: if %s is a file it does not exist!' % argv[i_input])
        
        # ---------#
        # module 1 #
        # ---------#
        if module == 1:
            path_emission_cdf = argv[2]     
            path_area_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            downscale_request = argv[8]
            
            # run module 1
            print('running module1')
            module1(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf, downscale_request)

            
        # ---------#
        # module 2 #
        # ---------#
        # elif module == 2:
        #     path_emission_cdf = argv[2]     
        #     path_nuts_cdf = argv[3]
        #     path_reduction_txt = argv[4]
        #     path_model_cdf = argv[5]
        #     path_result_cdf = argv[6]
            
        #     module2(path_emission_cdf, path_nuts_cdf, path_reduction_txt, path_model_cdf, path_result_cdf)

        # ---------#
        # module 3a #
        # ---------#
        elif module == 31:
            path_emission_cdf = argv[2]     
            path_area_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            module3a(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)

        # ----------#
        # module 3b #
        # ----------#
        elif module == 32:
            path_emission_cdf = argv[2]     
            path_area_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            module3b(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)
            
        # ---------#
        # module 4 #
        # ---------#
        elif module == 4:
            path_emission_cdf = argv[2]     
            path_area_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            module4(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)
            
        # ---------#
        # module 5 #
        # ---------#
        elif module == 5:
            path_emission_cdf = argv[2]     
            path_nuts_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            module5(path_emission_cdf, path_nuts_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)

        # ---------#
        # module 6 #
        # ---------#
        elif module == 6:
            path_emission_cdf = argv[2]     
            path_nuts_cdf = argv[3]
            cell_lat = argv[4]
            cell_lon = argv[5]
            path_reduction_txt = argv[6]
            path_base_conc_cdf = argv[7]
            path_model_cdf = argv[8]
            path_result_cdf = argv[9]
            
            module6(path_emission_cdf, path_nuts_cdf, cell_lat, cell_lon, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)

        # ---------#
        # module 7 #
        # ---------#

        
        elif module == 7:

            path_emission_cdf = argv[2]     
            path_base_conc_cdf = argv[3]
            path_natural_dir = argv[4]
            path_model_cdf = argv[5]
            fua_intersect_dir = argv[6]
            nuts_intersect_dir = argv[7]
            dbf_dir = argv[8]
            target_list = argv[9]
            path_result_cdf = argv[10]
            path_logo = argv[11]
            aggr_zones_in = argv[12]
            pollutant = argv[13]
            
            module7(path_emission_cdf, path_base_conc_cdf, path_natural_dir,
                path_model_cdf, fua_intersect_dir, nuts_intersect_dir, dbf_dir, 
                target_list, path_result_cdf, path_logo, aggr_zones_in, pollutant)

        # ---------#
        # module 8 #
        # ---------#
        
        elif module == 8:
            path_healthbl = argv[2]     
            path_result_cdf = argv[3]
            path_config_json = argv[4]
            if argv[5]:
                path_base_conc_cdf = argv[5]           
            module8_healthia(path_healthbl, path_result_cdf, path_config_json, path_base_conc_cdf)
       
        # ---------#
        # module 9 #
        # ---------#
        
        elif module == 9:
            path_aggrinp_txt = argv[2]     
            module9_aggregation(path_aggrinp_txt)
            
        elif module == 10:
            path_input_txt = argv[2]     
            path_datab_txt = argv[3]
            path_costd_txt = argv[4]

            module10(path_input_txt,path_datab_txt,path_costd_txt)      
     
     
        else:
            print('unknown module %d' % module)
        

    pass



