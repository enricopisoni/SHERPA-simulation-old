'''
Created on Jun 23, 2015

In this module the concentrations are calculated for a given emission reductions scenario
Inputs: - baseline emissions (per precursor and cell),
        - a netcdf defining the area where emission reduction will be applied
        - emission reductions (per precursor and macrosector)
        - model coefficients (per pollutant, precursor, cell)
        - path where the results will be written
        - eventually a progress log to be able to report on the progress of the whole calculation when module 1
        is called by another moduel
        - baseline concentrations. The are needed for NO2. First the NOx reduction is calculated and then a correlation
        is applied that predicts the NO2 fraction. Finally the delta NO2 is calculated
        
output: - netcdf with concentration changes per pollutant and cell
        - delta emission netcdf with emission changes per precursor and cell
       
The calculation is optimized using a flat weight over the whole domain. This allows to update only the scale
factor of this flat weight. The bell shape central weighting factors have to be recalculated for each cell.
        
@author: degraba
Enrico agrees with this nice explanation of module 1
'''

# imports
from math import isnan
from multiprocessing import shared_memory
import sys
from time import time
from netCDF4 import Dataset
from numpy import lib, zeros, sum, power, sqrt
from numpy.ma import is_masked
import multiprocessing as mp
from utils import (create_emission_reduction_dict,
                   create_emission_dict, create_window, read_progress_log,
                   deltaNOx_to_deltaNO2)
# EP 20210518
from globals import sector_lst
from tqdm import tqdm
import pickle

# Window class that returns aggregated weighting windows for a given omega


class OmegaPowerWindows:
    def __init__(self, hires_window_size):
        self.hires_window_size = hires_window_size
        self.hires_window_radius = (hires_window_size - 1) / 2
        self.hires_inverse_distance = zeros(
            (self.hires_window_size, self.hires_window_size))
        for iw in range(self.hires_window_size):
            for jw in range(self.hires_window_size):
                cell_dist = sqrt((float(iw - self.hires_window_radius))
                                 ** 2 + (float(jw - self.hires_window_radius)) ** 2)
                self.hires_inverse_distance[iw, jw] = 1 / (1 + cell_dist)

        # dictionary to store previously calculated windows elevated to omega
        self.hires_omega_windows = {}

    def getOmegaPowerWindow(self, omega):
        if omega in self.hires_omega_windows.keys():
            # for this omega the weighting window has been calculated
            pass
        else:
            self.hires_omega_windows[omega] = power(
                self.hires_inverse_distance, omega)

        return self.hires_omega_windows[omega]

# function that applies reductions per snap sector and precursor to the emission netcdf


def create_delta_emission(path_emission_cdf, precursor_lst, path_area_cdf,
                          path_reduction_txt, path_result_cdf,
                          write_netcdf_output, pollName, downscale_request):
    """
    Function that applies reductions per snap sector and precursor to the
    emission netcdf.
    """
    # create a dictionary with reductions per precursor and macro sector
    emission_reduction_dict = create_emission_reduction_dict(
        path_reduction_txt)

    # 20220530 if you do downscaling, only reduce PPM
    if downscale_request == 1:
        emission_reduction_dict['NOx'] = {
            1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0}
        emission_reduction_dict['NMVOC'] = {
            1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0}
        emission_reduction_dict['NH3'] = {
            1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0}
        emission_reduction_dict['SOx'] = {
            1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0}

    # open the emission netcdf
    emission_dict = create_emission_dict(path_emission_cdf, precursor_lst)

    # open the area netcdf
    rootgrp = Dataset(path_area_cdf, 'r')
    reduction_area = rootgrp.variables['AREA'][:] / 100.0
    rootgrp.close()

    # calculate a dictionary with the emission reductions per pollutant, macrosector and position
    delta_emission_dict = {}
    for precursor in precursor_lst:
        delta_emission_dict[precursor] = zeros(emission_dict[precursor].shape)
        # calculate the emission reduction
        # reductions are positive!
        # make the sum over all snap sectors
        # for snap in range(1, 13):
        # EP20210518
        for snap in range(1, sector_lst[-1]):
            delta_emission_dict[precursor][snap - 1, :, :] = emission_dict[precursor][snap -
                                                                                      1] * reduction_area * emission_reduction_dict[precursor][snap]

    # before summing over all snap sectors write the delta emissions per precursor and snap to a netcdf
    # create an output netcdf with delta emissions
    # --------------------------------------------
    if write_netcdf_output == True:
        filename_delta_emission_cdf = path_result_cdf + \
            'DCemis_emepV434_camsV42_' + pollName + '.nc'

        # change name of emission file in case of downscaling
        if downscale_request == 1:
            filename_delta_emission_cdf = path_result_cdf + \
                'DCemis_emepV434_camsV42_' + \
                pollName[0:8] + 'P' + pollName[8:] + '_.nc'

        rootgrp = Dataset(filename_delta_emission_cdf,
                          'w', format='NETCDF3_CLASSIC')

        # create dimensions in the netcdf file
        rootgrp.createDimension('latitude', len(emission_dict['lat_array']))
        rootgrp.createDimension('longitude', len(emission_dict['lon_array']))
        rootgrp.createDimension('GNFRsector', len(emission_dict['GNFRsector']))
        latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
        latitudes.units = "degrees_north"
        latitudes[:] = emission_dict['lat_array']
        longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
        longitudes.units = "degrees_east"
        longitudes[:] = emission_dict['lon_array']
        GNFRsector = rootgrp.createVariable(
            'GNFRsector', 'f4', ('GNFRsector',))
        GNFRsector[:] = emission_dict['GNFRsector']

        # 20220413, units for the emission file
        unitsForEmis = Dataset(path_emission_cdf, 'r').variables['NOx'].units

        # create delta emission data
        for precursor in precursor_lst:
            delta_emission_precursor = rootgrp.createVariable(
                precursor, 'f4', ('GNFRsector', 'latitude', 'longitude',))
            delta_emission_precursor.units = unitsForEmis  # "Mg/km2"
            delta_emission_precursor[:] = delta_emission_dict[precursor]

        rootgrp.close()

    # sum over all snap sectors
    for precursor in precursor_lst:
        delta_emission_dict[precursor] = sum(
            delta_emission_dict[precursor], axis=0)

    return delta_emission_dict


# function definition of source receptor model
def module1(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf, downscale_request, *progresslog):
    global shared_dictionary

    pollName = path_model_cdf.split('SR_')[1].split('.nc')[0]

    # check if a progess log file was passed as argument
    if progresslog:
        progress_dict = read_progress_log(progresslog[0])
        write_netcdf_output = False
    else:
        progress_dict = {'start': 0.0, 'divisor': 1.0}
        write_netcdf_output = True

    # M read the model netcdf
    # ---------------------
    rootgrp = Dataset(path_model_cdf, 'r')
    longitude_array = rootgrp.variables['lon'][0, :]
    latitude_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(longitude_array)  # len(rootgrp.dimensions['longitude'])
    n_lat = len(latitude_array)  # len(rootgrp.dimensions['latitude'])

    #####
    # 20180129 - EP - generalization to read 'radius of influence' variable, both written in matlab or python
    for nameatt in Dataset(path_model_cdf, 'r').ncattrs():
        if nameatt[0:6] == 'Radius':
            radiusofinfluence = nameatt
    inner_radius = int(
        getattr(Dataset(path_model_cdf, 'r'), radiusofinfluence))
    #inner_radius = int(getattr(rootgrp, 'Radius of influence'))
    #####

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

    # calculate the delta emissions, dictionary per pollutant a matrix of dimension n_lat x n_lon
    delta_emission_dict = create_delta_emission(path_emission_cdf, precursor_lst, path_area_cdf, path_reduction_txt, path_result_cdf, write_netcdf_output,
                                                pollName, downscale_request)

    # make a window
    window = create_window(inner_radius)
    (n_lon_inner_win, n_lat_inner_win) = window.shape

    pad_delta_emission_dict = {}
    for precursor in precursor_lst:
        pad_delta_emission_dict[precursor] = lib.pad(
            delta_emission_dict[precursor], inner_radius, 'constant', constant_values=0)

    # apply source receptor relationships
    # -----------------------------------
    last_progress_print = time()
    delta_conc = zeros((n_lat, n_lon)) * float('nan')
    cell_counter = 0
    n_cell = n_lat * n_lon

    # initialize a OmegaPowerWindows class
    win_pow_omega = OmegaPowerWindows(2 * inner_radius + 1)

    # loop over all cells of the domain
    shared_dictionary = {
        'alpha_dict': alpha_dict,
        'omega_dict': omega_dict,
        'pad_delta_emission_dict': pad_delta_emission_dict,
        'n_lon_inner_win': n_lon_inner_win,
        'n_lat_inner_win': n_lat_inner_win,
        'win_pow_omega': win_pow_omega,
        'n_lon': n_lon,
        'precursor_lst': precursor_lst
    }
    f = open('dict', 'wb')
    bytes_dict = pickle.dumps(shared_dictionary)
    sm = shared_memory.SharedMemory(create=True, size=len(bytes_dict))
    sm.buf[0:] = bytes_dict
    f.close()
    pool = mp.Pool(initializer=init, initargs=(sm.name,), processes=6)
    res = list(tqdm(pool.imap(work, range(n_lat)),
               total=n_lat, desc="progress:", bar_format='{desc}{percentage:3.6f}'))
    for i in range(n_lat):
        delta_conc[i] = res[i]
    pool.close()
    pool.join()
    sm.close()
    sm.unlink()
    # In the case of NO2 the variable 'delta_conc' contains the NOx concentrations as NO2-equivalent.
    # NO2 concentration and concentration difference are calculated applying an empiric formula
    # check if the pollutant is NO2, if so NO2 has to be calculated from NOx results w/ function 'deltaNOx_to_deltaNO2'
    if (path_model_cdf.find('NO2eq') > -1):
        rootgrp = Dataset(path_base_conc_cdf, 'r')
        base_conc_nox = rootgrp.variables['conc'][:]
        base_conc_no2 = rootgrp.variables['NO2'][:]
        rootgrp.close()
        delta_conc = deltaNOx_to_deltaNO2(
            delta_conc, base_conc_nox, base_conc_no2)

    # 20210202EP, read units to save correct units as output, for delta_conc
    rootgrp = Dataset(path_base_conc_cdf, 'r')
    units_for_output = rootgrp.variables['conc'].units

    # create a result netcdf
    # -----------------------
    if write_netcdf_output == True:
        filename_result_cdf = path_result_cdf + \
            'DCconc_emepV434_camsV42_' + pollName + '.nc'

        # 20220530 change name of emission file in case of downscaling
        if downscale_request == 1:
            filename_result_cdf = path_result_cdf + 'DCconc_emepV434_camsV42_' + \
                pollName[0:8] + 'P' + pollName[8:] + '_.nc'

        rootgrp = Dataset(filename_result_cdf, 'w', format='NETCDF3_CLASSIC')

        # create dimensions in the netcdf file
        rootgrp.createDimension('latitude', n_lat)
        rootgrp.createDimension('longitude', n_lon)
        latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
        latitudes.units = "degrees_north"
        longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
        longitudes.units = "degrees_east"
        latitudes[:] = latitude_array
        longitudes[:] = longitude_array

        # create delta concentration data
        delta_conc_pol = rootgrp.createVariable(
            'delta_conc', 'f4', ('latitude', 'longitude',))
        # u'\u03BC'+'g/m'+ u'\u00B3' # 'ug/m3' '\u03bcg/m\u00B3'
        delta_conc_pol.units = units_for_output
        delta_conc_pol[:] = delta_conc

        rootgrp.close()

    # create a results object
    mod1_res = {}
    mod1_res['delta_conc'] = delta_conc
    mod1_res['delta_emis_dict'] = delta_emission_dict
    mod1_res['n_lat'] = n_lat
    mod1_res['n_lon'] = n_lon
    mod1_res['latitude_array'] = latitude_array
    mod1_res['longitude_array'] = longitude_array

    return mod1_res


def init(smm):
    global shared_dictionary
    sm = shared_memory.SharedMemory(name=smm)
    shared_dictionary = pickle.loads(sm.buf)
    sm.close()


def work(ie):
    d = shared_dictionary
    res = zeros(d['n_lon']) * float('nan')
    for je in range(d['n_lon']):
        for precursor in d['precursor_lst']:
            alpha_ij = d['alpha_dict'][precursor][ie, je]
            omega_ij = d['omega_dict'][precursor][ie, je]
            if not(is_masked(alpha_ij)):
                # if the model is available remove NaN value
                if isnan(res[je]):
                    res[je] = 0
                # select the window of emissions around the target cell
                emissions_centre = d['pad_delta_emission_dict'][precursor][ie:(
                    ie + d['n_lon_inner_win']), je:(je + d['n_lat_inner_win'])]
                # apply the weights to the emissions and sum them over the whole window
                weighted_emissions_centre = (
                    d['win_pow_omega'].getOmegaPowerWindow(omega_ij) * emissions_centre).sum()
                # sum the contribution of the precursor
                res[je] = res[je] + alpha_ij * weighted_emissions_centre
    return res


if __name__ == '__main__':

    # testing is know done in a separate script
    pass
