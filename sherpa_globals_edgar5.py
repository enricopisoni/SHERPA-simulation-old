'''
Created on Jul 14, 2015

define some global variables

@author: degraba
'''

# variabels for testing
# absolute emission per cell and macro sector
path_emission_cdf_test = '../input_EDGAR_EMEP/createBCemiGNFR/output/BC_emissions_pm25_no2_o3.nc'
#path_emission_cdf_test = '../input_EDGAR_EMEP/createBCemiGNFR/output/BC_emissions_pm10.nc'

# netcdf with cells where reductions have to be applied (value between 0 and 1)
# path_area_cdf_test = 'input_EDGAR_EMEP/EMI_RED_NUTS3_ITALY.nc'
path_area_cdf_test = '../input_EDGAR_EMEP/createRedArea/output/emiRedOn_Germany.nc'#London_emepCams_0_100_FLIP.nc'
# reductions per precursor and macro sector
path_reduction_txt_test = '../input_EDGAR_EMEP/createRedText/user_reduction_GNFR_all.txt'#user_reduction_GNFR_all.txt'
path_reduction50all_txt_test = '../input_EDGAR_EMEP/user_reduction_GNFR_50p.txt'
# reductions per precursor and macro sector for module 3a and 3b
path_reduction_mod3a1P_txt_test = '../input_EDGAR_EMEP/user_reduction_GNFR_50p.txt'
path_reduction_mod3a2P_txt_test = '../input_EDGAR_EMEP/user_reduction_GNFR_50p.txt'
path_reduction_mod3b_txt_test = '../input_EDGAR_EMEP/user_reduction_GNFR_50p.txt'

# netcdf with model parameters per cell
#path_model_cdf_test = '../input_EDGAR_EMEP/SRR/SR_SURF_ug_PM25_rh50.nc' 
#path_model_cdf_test = '../input_EDGAR_EMEP/SRR/SR_SURF_ug_PM10_rh50.nc'
#path_model_cdf_test = '../input_EDGAR_EMEP/SRR/SR_SURF_ug_NO2.nc'  
#path_model_cdf_test = '../input_EDGAR_EMEP/SRR/SR_SURF_ug_NO.nc'  
#path_model_cdf_test = '../input_EDGAR_EMEP/SRR/SR_SURF_ppb_O3.nc'  
path_model_cdf_test = '../input_EDGAR_EMEP/SRR/SR_SURF_MAXO3.nc'  

# folder where output will be put
path_result_cdf_test = '../output_EDGAR_EMEP/'

# progress log is used when module 1 is called by another module
path_nuts0_cdf_test = '../input_EDGAR_EMEP/createRedArea/output/EMI_RED_ATLAS_NUTS_Lv0_EdgarEmep.nc'
path_nuts1_cdf_test = '../input_EDGAR_EMEP/createRedArea/output/EMI_RED_ATLAS_NUTS_Lv1_EdgarEmep.nc'
path_nuts2_cdf_test = '../input_EDGAR_EMEP/createRedArea/output/EMI_RED_ATLAS_NUTS_Lv2_EdgarEmep.nc'
path_nuts3_cdf_test = '../input_EDGAR_EMEP/createRedArea/output/EMI_RED_ATLAS_NUTS_Lv3_EdgarEmep.nc'

#path_base_conc_cdf_test = '../input_EDGAR_EMEP/createBCConc/output/BCconc_SURF_ug_PM25_rh50.nc'
#path_base_conc_cdf_test = '../input_EDGAR_EMEP/createBCConc/output/BCconc_SURF_ug_PM10_rh50.nc'
#path_base_conc_cdf_test = '../input_EDGAR_EMEP/createBCConc/output/BCconc_SURF_ug_NO2.nc'
#path_base_conc_cdf_test = '../input_EDGAR_EMEP/createBCConc/output/BCconc_SURF_ug_NO.nc'
#path_base_conc_cdf_test = '../input_EDGAR_EMEP/createBCConc/output/BCconc_SURF_ppb_O3.nc'
path_base_conc_cdf_test = '../input_EDGAR_EMEP/createBCConc/output/BCconc_SURF_MAXO3.nc'

#for health
path_healthbl_test = '../input_EDGAR_EMEP/createImpactFileForSherpa/input/impacts/healthbl_nc.nc'
path_config_json_test = '../input_EDGAR_EMEP/config/sharedvariables.json'

fua_intersect_dir = '../input_EDGAR_EMEP/createGridIntersect/output_mod7/fua/'
nuts_intersect_dir = '../input_EDGAR_EMEP/createGridIntersect/output_mod7/nuts/'
dbf_dir = '../input_EDGAR_EMEP/createGridIntersect/output_mod7/'
target_list = '../input_EDGAR_EMEP/createGridIntersect/output_mod7/AM_targets.txt'
path_natural_dir_test = '../input_EDGAR_EMEP/createDustSalt/output/'
aggr_zones='fua'
path_logo_test=''
aggrinp_txt=''
# list of precursors
# order important, it's the order in the alpha and omega arrays
# precursor_lst = ['NOx', 'NMVOC', 'NH3', 'PM25', 'SOx']  

# order important, it's the order in the alpha and omega arrays
sector_lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]   #  

# fixed reduction percentage for potency calculation
alpha_potency = float(50)
