'''
Created on Jul 14, 2015

auxiliary functions for SHERPA
@author: degraba
'''

from netCDF4 import Dataset
import numpy as np
from numpy import zeros, sqrt, array
import pandas as pd
import sys
#EP 20210518
from sherpa_globals import sector_lst

# create a dictionary with emission reductions per precursor and nuts
#--------------------------------------------------------------------

# file format:
# POLL    MS1    MS2    MS3    MS4    MS5    MS6    MS7    MS8    MS9    MS10
# NOx    0    0    0    0    0    100    0    0    0    0
# NMVOC    0    0    0    0    0    100    0    0    0    0
# NH3    0    0    0    0    0    100    0    0    0    0
# PM25    0    0    0    0    0    100    0    0    0    0
# SOx    0    0    0    0    0    100    0    0    0    0

def create_emission_reduction_dict(path_reduction_txt):
    # read emission reductions per precursor and macro sector
    f = open(path_reduction_txt, 'r')
    emission_reduction_dict = {}
    f.readline()
    while True:
        line = f.readline().rstrip()
        if len(line) == 0:
            break
        value_lst = line.split('\t')
        # sub dictionary per precursor
        precursor = value_lst[0]
        emission_reduction_dict[precursor] = {}
        #EP20210518
        for snap in range(1, sector_lst[-1]):
#        for snap in range(1, 13):
            emission_reduction_dict[precursor][snap] = float(value_lst[snap]) / 100.0
    f.close()

    return emission_reduction_dict

# function making a list of the precursors that are reduced
def create_reduced_precursor_lst(emission_reduction_dict):
    reduced_precursor_lst = []
    for precursor in emission_reduction_dict.keys():
        sum_reductions = 0
        for snap in emission_reduction_dict[precursor].keys():
            sum_reductions += emission_reduction_dict[precursor][snap]
        if sum_reductions > 0:
            reduced_precursor_lst.append(precursor)
    return reduced_precursor_lst


# create a dictionary with emissions per precursor, macrosector and postion (lat, lon)
#-------------------------------------------------------------------------------------
def create_emission_dict(path_emission_cdf, precursor_lst):
    # open the emission netcdf
    rootgrp = Dataset(path_emission_cdf, 'r')

    emission_dict = {}
    for precursor in precursor_lst:
        emission_dict[precursor] = rootgrp.variables[precursor][:, :, :]

    # get snap, longitude and latitude arrays from emission file

    #snap_array = range(1, 13)
    #EP20210518
    snap_array = range(1, sector_lst[-1])
    lon_array = rootgrp.variables['longitude'][:]
    lat_array = rootgrp.variables['latitude'][:]
    emission_dict['GNFRsector'] = snap_array
    emission_dict['lon_array'] = lon_array
    emission_dict['lat_array'] = lat_array

    # close the emission file
    rootgrp.close()

    return emission_dict


# make a window with cell distances to the central cell
# -----------------------------------------------------
def create_window(radius):
    # the window contains the distance between each cell and the centre of the window
    # the distance is expressed in cells
    n_lon_win = 2 * radius + 1
    n_lat_win = 2 * radius + 1

    window = zeros((n_lat_win, n_lon_win))
    i_centre = radius
    j_centre = radius
    for iw in range(n_lon_win):
        for jw in range(n_lon_win):
            cell_dist = sqrt((float(iw - i_centre)) ** 2 + (float(jw - j_centre)) ** 2)
            window[iw, jw] = 1 / (1 + cell_dist)

    return window

# convert to progress log file to a dictionary
def read_progress_log(progresslog):
    progress_dict = {}
    f_prog = open(progresslog, 'r')
    line = f_prog.readline().rstrip()
    [start, divisor] = line.split('\t')
    progress_dict['start'] = float(start)
    progress_dict['divisor'] = float(divisor)
    progress_dict['netcdf_output'] = False

    return progress_dict

# write progress log file
def write_progress_log(progress_log_filename, start, divisor):
    # write progress log file
    f_prog = open(progress_log_filename, 'w')
    f_prog.write('%f\t%f' % (start, divisor))
    f_prog.close()

# define a function that applies the NO2 fraction correlation
def fno2_corr(nox_array):
    # this is the correlation used in GRAL
    fno2_array = 30 / (nox_array + 35) + 0.18
    return fno2_array


# function that converts a delta_conc(NOx) into a delta_conc(NO2)
def deltaNOx_to_deltaNO2(delta_conc_nox, base_conc_nox, base_conc_no2):
    # read the NO2 concentration (should be inside the NO2eq/NOx file)
    base_fno2 = base_conc_no2 / base_conc_nox
    base_fno2_rel_error = fno2_corr(base_conc_nox) / base_fno2

    # calculate NO2 fraction and the absolute NO2 concentration
    # delta_conc = -(scen_conc_nox - base_conc_nox)
    scen_conc_nox = base_conc_nox - delta_conc_nox
    # the NO2 fraction given by the correlation has to be corrected with the NO2 fraction of each cell
    # from the baseline scenario. Otherwise delta_NO2's will be created when the NO2 fraction is different for the
    # correlation and the basecase model results.
    scen_fno2 = array(fno2_corr(scen_conc_nox) / base_fno2_rel_error)
    # correlation can lead to NO2 fractions above 1, avoid this
    scen_fno2[scen_fno2 > 1] = 1
    scen_conc_no2 = scen_conc_nox * scen_fno2
    # recalculate delta_conc
    delta_conc_no2 = base_conc_no2 - scen_conc_no2

#     # add diagnostic variables in the case of NOx
#     delta_conc_nox_var = rootgrp.createVariable('delta_conc_nox', 'f4', ('latitude', 'longitude',))
#     delta_conc_nox_var.units = 'ug/m3'
#     delta_conc_nox_var[:] = delta_conc_nox

    return delta_conc_no2

def read_nc(nc_file):
    '''
    NAME
        Reads SHERPA ncdf file with Python
    PURPOSE
        To read matrix data and put them in a multindexed dataframe
    PROGRAMMER(S)
        Denise Pernigotti
    REVISION HISTORY
        20/02/2018 Emanuela Peduzzi    
    REFERENCES
    
    '''
    nc_data = Dataset(nc_file, 'r')
    nc_dims = [dim for dim in nc_data.dimensions]
    nc_vars = [var for var in nc_data.variables]
    #sometimes the latitude is written just with lat as in model data
    latname=list(filter(lambda x: x in nc_vars, ['Lat','lat','latitude']))[0]
    #latname=list(filter (lambda x: 'lat' in x, nc_vars))[0]
    lats = nc_data.variables[latname][:]
    lonname=list(filter(lambda x: x in nc_vars, ['Lon','lon','longitude']))[0]
    #lonname=list(filter (lambda x: 'lon' in x, nc_vars))[0]
    lons = nc_data.variables[lonname][:]
    #if there are three dimensional arrays
    if len(nc_dims)==3:
        ncz=str(list(set(nc_dims)-set(['latitude','longitude']))[0])
        nz=range(len(nc_data.dimensions[ncz]))
        if ncz=='pollutant':
            strpoll=nc_data.Order_Pollutant
            nznames=strpoll.split(', ')
        else:
            nznames=[ncz +"{:02d}".format(x+1) for x in nz]
            #nznames=[ncz + s for s in map(str,range(1,len(nc_data.dimensions[ncz])+1))]
    #create an index with lat and lon
    #latrep=map(str, np.repeat(lats,len(lons)))
    #lonrep=map(str, np.tile(lons,len(lats)))
    #trasform variables arrays in vectors
    #allvar={'lat_lon':map(lambda (x,y): x+'_'+y, zip(latrep, lonrep))}
    #create lat and lon info
    if len(lats.shape)==2 and len(lons.shape)==2:
        nrow=lats.shape[0]
        ncol=lats.shape[1]
        lon=lons.ravel()
        lat=lats.ravel()
    else:
        nrow=len(lats)
        ncol=len(lons)
        lon=np.tile(lons,nrow)
        lat=np.repeat(lats,ncol)

    y=np.repeat(range(1, nrow+1),ncol)
    x=np.tile(range(1, ncol+1),nrow)
    row=list(map(str,y))
    col=list(map(str,x))
    index_grid=list(map(lambda x: '_'.join(x),list(zip(col,row))))

    allvar={}
    allvar['coord']=pd.DataFrame(lon,columns=['lon'])
    allvar['coord']['lat']=lat
    allvar['coord']['x']=x
    allvar['coord']['y']=y
    allvar['coord'].index=index_grid
    nc_vars.remove(latname)
    nc_vars.remove(lonname)
    # added by EPE to deal with the delta_emissions file created byt 
    # the GUI which has an extra variable 
    # @todo this condition can be removed in the future if the GUI 
    # is fixed (20180219)
    if 'GNFRsector' in nc_vars: 
        nc_vars.remove('GNFRsector')
        
    for var in nc_vars:
        varnc=nc_data.variables[var][:]
        if len(nc_dims)==3:
            allvar[var]=pd.concat(map(lambda sn : pd.Series(varnc[sn].ravel()),nz),axis=1)
            allvar[var].columns=nznames
        else:
            allvar[var]=pd.DataFrame(varnc.ravel())
            allvar[var].columns=[var]
        allvar[var].index=index_grid
        #allvarnc[var]=allvarnc[var].transpose()
        #index_var = pd.MultiIndex.from_tuples(zip(np.repeat(var,len(nz)),nznames), names=['vars', ncz])
        #allvar[var].columns=index_var
    nc_data.close()
    reform = {(outerKey, innerKey): values for outerKey, innerDict in allvar.items() for innerKey, values in innerDict.items()}
    df=pd.DataFrame(reform)
    return df.transpose()

def read_nuts_area(filenuts, calcall=False, nullnut=False, nutsall=None):
    '''
    NAME
        Import info on grid points attribution to nuts or specific area type from ascii file
    PURPOSE
        Import info on grid points attribution to nuts or specific area type from ascii file/s.
        If the file is single then it must contain the column 'Area [km2]' relative to % of the area in the finest nut,
        this datum will be set to each nut but it will then aggregated for larger nuts when nutsarea will be calculated
        If the files are two, then each nut will have its own % area for each grid point, then the data will be merged here
    PROGRAMMER(S)
        Denise Pernigotti
    REVISION HISTORY
        WARNING by EPE: this function does not work andymore as originally 
        intended because the gridintersect structure has changed
        for example there is no nullnut  which has the same name for all 
        cells! the 'rect' par has not been tested.
        EPE: added reading population. 
    REFERENCES
    
    '''   
#    filenuts='D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/selection/grid_intersect'
#    filenuts='D:/sherpa.git/Sherpa/input/selection/gridnew/fua/grid_intersect'
    nuts_info_all={}
    if(filenuts != 'rect'):
        nuts_def= filenuts +'.txt'
        nuts_info = pd.read_csv(nuts_def,delimiter="\t") # EPE added dtype because of warning message
#        pd.to_numeric(nuts_info['POPULATION'])
        nuts_info=nuts_info.dropna(axis=1,how='all')
        nutsnames=list(nuts_info.columns[~nuts_info.columns.isin(['POP','COL','ROW','AREA_km2','LAT','LON','CENTROID_X', 'CENTROID_Y', 'PERCENTAGE', 'POPULATION'])])
        #optional 'nut' comprising all grid points
        if calcall :
        #nutsnames.insert(0, 'ALL')
            nutsnames.insert(0, 'ALL_'+nutsnames[0])
            nuts_info[nutsnames[0]]=nutsnames[0]
        nuts_info['grid']=['_'.join(str(i) for i in z) for z in zip(nuts_info['COL'],nuts_info['ROW'])]
        if 'AREA_km2' in nuts_info.columns:
            nuts_area=pd.concat(map(lambda p: nuts_info['AREA_km2'],nutsnames),axis=1)
            #nuts_area.index=nuts_info['grid']
            nuts_area.columns=nutsnames
        if 'POPULATION' in nuts_info.columns:
            nuts_pop=pd.concat(map(lambda p: nuts_info['POPULATION'],nutsnames),axis=1)
            #nuts_area.index=nuts_info['grid']
            nuts_pop.columns=nutsnames
           #nuts_info=nuts_info[nutsnames]
        else:
            sys.exit("missing infos on grid cells area per nut")

        #aggregate data for each nut, create a dictionary
        nut_info_nut={}
        nut_info={}
#        nut_info_area={}
#        nut_info_pop={}
        nullnut_dct={'NUTS_Lv0':'0', # at this level there is no background
                     'NUTS_Lv1':'1', 
                     'NUTS_Lv2':'11',
                     'NUTS_Lv3':'111'}
        for nut in nutsnames:
            #create a multindex
            index = pd.MultiIndex.from_tuples(list(zip(nuts_info[nut],nuts_info['grid'])), names=['nutname','grid'])
            nut_info_area=pd.Series(list(nuts_area[nut]), index=index)
            nut_info_area=nut_info_area.to_frame(name='area')
            nut_info_pop=pd.Series(list(nuts_pop[nut]), index=index)
            nut_info_pop=nut_info_pop.to_frame(name='pop')       
            nut_info=pd.concat([nut_info_area, nut_info_pop], axis=1)

            #aggregate data on these nuts if necessary
            nut_info_nut=nut_info.groupby(level=[0,1]).sum()
            #find total area
            grid_area_tot=nut_info_nut.groupby(level=['grid']).sum()
            nut_info_nut['parea']=nut_info_nut['area']/grid_area_tot['area']
            nut_info_nut.loc[nut_info_nut['area']==0,'parea']=0.
            #eventually remove the fillng code
            if nullnut is True:
                for nutkey in set(nut_info_nut.index.get_level_values(0)):
                    if nutkey[2:]==nullnut_dct[nut]:
#                        print(nutkey)
                        nut_info_nut=nut_info_nut.drop(nutkey, level='nutname')
            nuts_info_all[nut]=nut_info_nut.reindex(columns =['area', 'parea', 'pop'])
       
    else:
        nuts_rect=nutsall
        nuts_rect.index=nuts_rect.index.droplevel(level=0)
        grid_inrect=nuts_rect.index
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lon']>=rect_coord['ll']['lon']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lon']<=rect_coord['ur']['lon']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lat']>=rect_coord['ll']['lat']]
        grid_inrect=grid_inrect[coordinates.loc[grid_inrect,'lat']<=rect_coord['ur']['lat']]
        nuts_rect=nuts_rect.loc[list(grid_inrect)]
        nuts_rect['nutname'] = 'rect'
        nuts_rect.set_index('nutname', append=True, inplace=True)
        nuts_info_all['rect']=nuts_rect.swaplevel(i=-2, j=-1, axis=0)
   
    return nuts_info_all



if __name__ == '__main__':

    # check the window function
    radius = 200
    testwindow = create_window(radius)
    window_file = open('C:/temp/source_recptor_window.txt', 'w')
    for i in range(2 * radius + 1):
        for j in range(2 * radius + 1):
            window_file.write('%e\t' % testwindow[i,j])
        window_file.write('\n')
    window_file.close()
    pass


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
