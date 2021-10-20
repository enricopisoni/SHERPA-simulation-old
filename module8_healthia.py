# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:32:28 2017
    
    Health impacts of PM2.5 - WHO-Europe method (HRAPIE reccomendations)
    to calculate mortality according to the Risk Rates
   
    Corresonding YLLs or days of life lost are calculated considering the 
    distribution of mortality and population by age and by country, 
    where data is not available in the ICD-10 format YLLs are calculated 
    considering the average value for the countries that 
    are availble.
    
    ASSUMPTION: baseline values (i.e. mortality and average years of life 
    loss) are averaged in border cells between countries.
    
    NB: A positive delta means a reduction!
    
    INPUT: 
        - path_healthbl: baseline values for the impact calculation \
          produced by precompute_healthia.py
        - path_config_json_test: configuration file that the SHERPA interface 
          will also use (if it is not there default values are taken)
        - path_base_conc_cdf_test: optional argument (it is not needed when 
          using the module from the interface)
    OUTPUT: 
        - healthia.nc 
        
    - Bibliography:
    
    [1] Estimating Local Mortality Burdens associated with Particulate Air 
    Pollution 2014 Public Health England 
    
    [2] World Health Organization Europe, 2013. Health risks of air pollution
    in Europe - HRAPIE project - Recommendations for concentration–response
    functions for cost–benefit analysis of particulate matter, ozone and
    nitrogen dioxide, Copenhagen Ø, Denmark.

    [3] Holland, M., 2014. Cost-benefit Analysis of Final Policy Scenarios
    for the EU Clean Air Package Version. Version 2

    [4] World Health Organization Europe, 2017. AirQ+: software tool for health
    risk assessment of air pollution.
     
    [5] Holland, M., 2014. Implementation of the HRAPIE Recommendations for 
    European Air Pollution CBA work. EMRC. \
      
    [6] Data for baseline population 
    
        ICD codes: ICD-10: A00-B99,C00-D48,D50-D89,E00-E88,F01-F99,G00-G98,
        H00-H59,H60-H93,I00-I99,J00-J98,K00-K92,L00-L98,M00-M99,N00-N98,
        O00-O99,P00-P96,Q00-Q99,R00-R99
        Age: '30 - 85 +'
        Sex: Both
        http://data.euro.who.int/dmdb/ [Accessed December 13, 2016].

    @author: peduzem
    """


from netCDF4 import Dataset  
import numpy as np
import os as os
import json

#
#from sherpa_globals import (path_result_cdf_test,
#                            path_healthbl_test, path_config_json_test, path_base_conc_cdf_test,
#                            )

def health_impact(pop30plus, pm25_conc, ar_drate, ar_lyl, approx='l'):

    """
    Function that caclulates the health impact
    
    INPUT : 
        - pop30plus = array with the distribution of the population over 30 \
                      years of age
        - pm25_conc = array with antrhopogenic concentration of PM2.5 \
                      (total)
        - ar_drate = array with the distribution of baseline death rate \
                     (from all cause mortality)
        - ar_lyl = array with the average years of life lost per death \
                    over 30 years of age
        - approx = 'e' for exponential and 'l' for linear       
    
    OUTPUT :
        - mort = array with mortality (from lower bound to upper bound)
        - dll = array with the days of life lost per year \
                (from lower bound to upper bound)
        - dll_spec = array with the days of life lost per person per year \
                     (from lower bound to upper bound)
          
          
    @author: peduzem
    """
# -----------------------------------------------------------------------------
    # create empty arrays to store results
    mort = np.zeros(np.shape(pop30plus))
    dll = np.zeros(np.shape(pop30plus))
    
# -----------------------------------------------------------------------------
    # CONCENTRATION RESPONSE FUNCTION:
    # From [2] Table 1 
    # Estimate of mortality, all-cause (natural) age 30+ years
    # PM2.5 Annual mean 
    # RR = 1.062 (1.04-1.083) 95% CI per 10 microg/m3 
    lrr = 1.04
    mrr = 1.062
    hrr = 1.083
    
    # From [4] Beta: 
    lbeta = 0.003922071315328133  # lower bound
    mbeta = 0.006015392281974714  # average value
    hbeta = 0.007973496801885352  # higher bound
    # RR = e^(beta x)=(e^(beta*10))^(x/10) = 1.062^(x/10) CI = 1.04, 1.083
    # (we obtain the same values reported in [2])

    if approx == 'l':
    # Linear approximation: f(x-x0) = f(x0)+ f'(x0)(x-x0)
    #                       RR = 1 + e^(beta*10)^(x0/10)*ln(e^(beta*10))*(x-x0)/10
    #                       RR = 1 + ln(e^(beta*10))*(x-x0)/10 = 1 + coef*x/10
    # AF = (RR -1)/RR = coef*x/10 / (1+ coef*x/10)
    
    # linear approximation
        lcoef = np.log(lrr)
        mcoef = np.log(mrr)
        hcoef = np.log(hrr)

        coef = [lcoef, mcoef, hcoef]
        pt = len(coef)
        mort = mort + (np.where(np.isnan(pm25_conc), 0, (
                     [(coef[i]*pm25_conc/10)/(1+coef[i]*pm25_conc/10)*pop30plus*ar_drate
                     for i in range(len(coef))]))) 
    
    # exponential approximation
    elif approx == 'e':
    # AF = (RR -1)/RR = e^bx -1 / (e^bx) = 1 -e^(-bx) 
        beta = [lbeta, mbeta, hbeta]
        pt = len(beta)
        mort = mort + (np.where(
                    np.isnan(pm25_conc), 0, (
                            [(1-(np.exp(-beta[i]*pm25_conc))) *
                             pop30plus * ar_drate
                             for i in range(len(beta))]))) 

# -----------------------------------------------------------------------------
    # ESTIMATE OF THE YLL (Not in the Guidelines!)
    # days of life lost per year 
    dll = dll + (np.where(np.isnan(pm25_conc), 0,
                             [mort[i] * ar_lyl * 365
                             for i in range(pt)])) 
    # days of life lost per person per year 
    dll_spec = [np.divide(dll[i], pop30plus, out=np.zeros_like(dll[i]), where=pop30plus!=0) for i in range(pt)] 
    
# -----------------------------------------------------------------------------
    # return results    
    return mort, dll, dll_spec



def module8_healthia(path_healthbl, path_result_cdf, path_config_json, *path_base_conc_cdf):
    """
    Main functin that calculates the health impacts given the paths: 
    input: 
        - path_base_conc_cdf_test = base case concentration \
          optional input argument if value_conc is not in the results
        - path_dust_conc_cdf_test = path of baseline dust concentration 
        - path_salt_conc_cdf_test = path of baseline salt concentration 
        - path_healthbl = path where results are stored (health baseline)
        - path_result_cdf_test: path of the delta concentrations
           (output of module1) 
    @author: peduzem
    """
    # value of concentration from dust and salt (i.e. natural)
    fh_pm25_natural = Dataset(path_healthbl, mode='r')
    pm25_natural = fh_pm25_natural.variables['conc'][:]
    fh_pm25_natural.close()

    # delta concentration from model resutls
    path_conc_nc = path_result_cdf + 'delta_concentration.nc'
    fh_deltapm25 = Dataset(path_conc_nc, mode='r')
    d_pm25_conc = fh_deltapm25.variables['delta_concentration'][:]
#       pm25_delta = fh_deltapm25.variables['conc'][:]
    fh_deltapm25.close()
    
    # SHERPA interface produces the scenario nc file..     
    path_value_nc = path_result_cdf + 'value_conc.nc'
    # if it is not present the scenario concentration has to be calculated 
    # from the base concentration
    if not os.path.exists(path_value_nc):
#        if path_base_conc_cdf[0]:
        if path_base_conc_cdf[0]:
            print('Calculating scenario value from base case concentration')
            fh_pm25_base = Dataset(path_base_conc_cdf[0], mode='r')
            pm25_base = fh_pm25_base.variables['conc'][:]
            fh_pm25_base.close()
            pm25_conc = pm25_base - d_pm25_conc                      
        else: 
            print('Error')
    else: 
    # if the scenario value is in the results
        fh_pm25_conc = Dataset(path_value_nc, mode='r')
        pm25_conc = fh_pm25_conc.variables['conc'][:]
        fh_pm25_conc.close()
    
    # Anthropogenic concentration: scenario values minus natural concentration 
    # -- there are different views on this. At the moment natural background is 
    # substracted to the concentration. See for example: 
    # -- [1]	H. Fintan, A. Hunt, H. Cowie, M. Holland, B. Miller, S. Pye, and
    # -- P. Watkiss, “Service Contract for Carrying out Cost-Benefit Analysis
    # -- of Air Quality Related Issues, in particular in the Clean Air for 
    # -- Europe (CAFE) Programme,” 2005.
    # -- [2]   The ALPHA Benefit Assessment Model, EC4MACS 2013

    sce_pm25_conc = pm25_conc - pm25_natural  

    # get baseline data from nc file
    fh = Dataset(path_healthbl, mode='r')
    pop30plus = fh.variables['ppl30+'][:]
    fh.close()
    fh = Dataset(path_healthbl, mode='r')
    ar_drate = fh.variables['deathsppl30+'][:]
    fh.close()
    fh = Dataset(path_healthbl, mode='r')
    ar_lyl = fh.variables['lyl30+'][:]
    fh.close()
    
 # -----------------------------------------------------------------------------   
    # calculate impacts   
    sce_mort, sce_dll, sce_dll_spec = health_impact(pop30plus, sce_pm25_conc,
                                                    ar_drate, ar_lyl, approx='l')
    bc_pm25 = sce_pm25_conc + d_pm25_conc 
    
    # this could be improved (the valuse are always the same)
    bc_mort, bc_dll, bc_dll_spec = health_impact(pop30plus, bc_pm25,
                                                    ar_drate, ar_lyl, approx='l')
    delta_mort = bc_mort - sce_mort
    delta_dll = bc_dll - sce_dll
    delta_dll_spec = np.array(bc_dll_spec) - np.array(sce_dll_spec)

# -----------------------------------------------------------------------------
    # Saving results:        
    # default dictionary to save results:           
    dflt_dict = {
	"d_mort":{
		"impact": "Mortality",
		"data": "Delta",
		"ci":["d_mort_lb", "d_mort", "d_mort_ub"],

		"long_description":["delta mortality lower bound", "delta mortality", "delta mortality upper bound"],
		"aggregation":"sum",
		"units":"people/year"},
	"v_mort":{
		"impact": "Mortality",
		"data": "Value",
		"ci":["v_mort_lb", "v_mort", "v_mort_ub"],

		"long_description":["mortality lower bound", "mortality", "mortality upper bound"],
		"aggregation":"sum",
		"units":"people/year"},
	"d_dll":
	{
		"impact": "Days of life loss",
		"data": "Delta",
		"ci":["d_dll_lb", "d_dll", "d_dll_ub"],

		"long_description":["delta days of life loss lower bound", "delta days of life loss", "delta days of life loss upper bound"],
		"aggregation":"sum",
		"units":"dll/year"},
	"v_dll":
	{
		"impact": "Days of life loss",
		"data": "Value",
		"ci":["v_dll_lb", "v_dll", "v_dll_ub"],

		"long_description":["days of life loss lower bound", "days of life loss", "days of life loss upper bound"],
		"aggregation":"sum",
		"units":"dll/year"},
	"d_dll_pp":
	{
		"impact": "Days of life loss per person",
		"data": "Delta",
		"ci":["d_dll_pp_lb", "d_dll_pp", "d_dll_pp_ub"],

		"long_description":["delta days of life loss per person lower bound", "delta days of life loss per person", "delta days of life loss per person upper bound"],
		"aggregation":"population weighted average",
		"units":"dll/(person year)"},
	"v_dll_pp":
	{
		"impact": "Days of life loss per person",
		"data": "Value",
		"ci":["v_dll_pp_lb", "v_dll_pp", "v_dll_pp_ub"],

		"long_description":["days of life loss per person lower bound", "days of life loss per person", "days of life loss per person upper bound"],
		"aggregation":"population weighted average",
		"units":"dll/(person year)"}
    }

    # If it exists we use the json config file which is used also by the 
    # SHERPA interface
    if os.path.exists(path_config_json):    
        print('Using stored json file')
        json_file = open(path_config_json)
        json_str = json_file.read()
        cfg_dct = json.loads(json_str)
    else:
        print('Not using stored json file')
        cfg_dct = dflt_dict
    
    # Generation of results files: 
    outfile=path_result_cdf + 'healthimp.nc'
    if os.path.exists(outfile):
        os.remove(outfile)   
    for key in cfg_dct.keys():
        if key == 'd_mort':
            for it in enumerate(cfg_dct[key]['ci']):  
                write_nc(delta_mort[it[0]], outfile, it[1], cfg_dct[key]['units'], path_healthbl,
                     addnutsid=False, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'v_mort': 
            for it in enumerate(cfg_dct[key]['ci']):
                write_nc(sce_mort[it[0]], outfile, it[1], cfg_dct[key]['units'], path_healthbl,
                     addnutsid=False, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'd_dll':
            for it in enumerate(cfg_dct[key]['ci']):
                write_nc(delta_dll[it[0]], outfile, it[1], cfg_dct[key]['units'], path_healthbl,
                     addnutsid=False, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'v_dll': 
            for it in enumerate(cfg_dct[key]['ci']):  
                write_nc(sce_dll[it[0]], outfile, it[1], cfg_dct[key]['units'], path_healthbl,
                     addnutsid=False, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'd_dll_pp':
            for it in enumerate(cfg_dct[key]['ci']):  
                write_nc(delta_dll_spec[it[0]], outfile, it[1], cfg_dct[key]['units'], path_healthbl,
                     addnutsid=False, l_name=cfg_dct[key]['long_description'][it[0]])
        if key == 'v_dll_pp': 
            for it in enumerate(cfg_dct[key]['ci']):  
                write_nc(sce_dll_spec[it[0]], outfile, it[1], cfg_dct[key]['units'], path_healthbl,
                     addnutsid=False, l_name=cfg_dct[key]['long_description'][it[0]])
    
## SUPPORT FUNCTIONS (IDEALLY IN THE AUXIALIARIES FILE)

def write_nc(array, path_nc, name_var, unit_var, path_healthbl, addnutsid=False, l_name=None):
    ''' Function to write an array in a netcdf file,
        if the file already exist it is going to write in append mode, 
        otherwise in write mode. 
        input:
            - array: data to write
            - path_nc: path of netcdf file
            - name_var: name for data in array
            - unit_var: units for data in array
            - path_healthbl: ncd file template (for lon and lat arrays)
            - addnutsid: if True the layer nuts_id is added so that the
                nectcdf file is consistent with the ones provided
                by terraria
        ouput: 
            - nc file
    @author: peduzem
    '''
    rootgrp = Dataset(path_healthbl, 'r')
    lon_array = rootgrp.variables['longitude'][:]
    lat_array = rootgrp.variables['latitude'][:]
    rootgrp.close()

    if not os.path.exists(path_nc):
        mode = 'w' 
        fh=Dataset(path_nc, mode=mode, format='NETCDF3_CLASSIC') 
        fh.createDimension('latitude', len(lat_array))
        fh.createDimension('longitude', len(lon_array))
        latitude = fh.createVariable('latitude', 'f4', ('latitude',))
        longitude = fh.createVariable('longitude', 'f4', ('longitude',)) 
        if addnutsid is True:
#        fh.createDimension('z', 10)
            fh.createDimension('nuts_id', 1)
            var = fh.createVariable(name_var, 'f8',
                                    ('nuts_id', 'latitude', 'longitude',))
            nutsid = fh.createVariable('NUTS', 'i4', ('nuts_id',))
            longitude[:] = lon_array
            latitude[:] = lat_array
            nutsid[0] = 1
            var[0, :] = array
        elif addnutsid is False:
            longitude[:] = lon_array
            latitude[:] = lat_array
            var = fh.createVariable(name_var, 'f8', ('latitude', 'longitude'))
            var[:] = array          
    else:
        mode = 'a'
        fh=Dataset(path_nc, mode=mode, format='NETCDF3_CLASSIC')
        if addnutsid is True:
            var = fh.createVariable(name_var, 'f8',
                                    ('nuts_id', 'latitude', 'longitude',))
            var[0, :] = array
        elif addnutsid is False:
            var = fh.createVariable(name_var, 'f8', ('latitude', 'longitude'))
            var[:] = array

    fh.variables[name_var].units = unit_var
    if l_name is not None:
            fh.variables[name_var].long_name =l_name   
    fh.close()  


if __name__ == '__main__':
    
#    module8_healthia(path_healthbl_test, path_result_cdf_test, path_config_json_test, path_base_conc_cdf_test)

    pass
