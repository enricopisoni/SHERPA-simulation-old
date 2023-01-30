# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:32:03 2018

This module is used to do all the aggregations for the GUI, as a postcompute. 
it produces equivalent files to the old fortran code (but emissions are in Mg, 
and not densities)

NB: the values for mortality are treated as a special case. 

STILL MISSING: average over threshold 
@author: peduzem
"""
import pandas as pd
#import numpy as np
from time import time
import json
import ast 

from sherpa_auxiliaries import read_nuts_area 
from sherpa_auxiliaries import read_nc
#EP 20210518
from sherpa_globals import sector_lst


def area_calc(rw, lat_area_dct):
    ''' function to look up the area of the cell from the index used in the 
    dataframe. This workoround is because in the grid intersect the area of the 
    cells by the water refers to the fraction on land and not the whole cell!
    '''
    rw_cell=int(rw.split('_')[1])                     
    area_cell=lat_area_dct[rw_cell]
    return area_cell
                    
def module9_aggregation(aggrinp_txt):
    
    '''Function that aggregates results by area (for both nuts level and fuas)
    
    Inputs: 
        -aggrinp_txt: txt file with the information:
        {
        "delta": {
                "path": "D:/programs/sherpa/app/data/temp/delta_concentration.nc",
                "var": "delta_concentration",
                "aggregation": "('avg','pop')"},
        "bc": {
             "path": "D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/base_concentrations/BC_conc_PM25_Y.nc",
             "var": "conc",
             "aggregation": "('avg','pop')"},
     	"grid-intersect":"D:/programs/sherpa/app/data/input/models/chimere_7km_nuts/selection/grid_intersect",
    	"output-dir":"D:/programs/sherpa/app/data/temp/"
        }
    
    Warning for health impact: when pointing to base case it is actually pointing 
    to the scenario results. This should be improved harmonized. 
    
    Outputs:
        -txt files that contain the value
        'NUTS0', 'NUTS1','NUTS2', 'NUTS3' 
        
    '''

    json_file = open(aggrinp_txt)
    json_str = json_file.read()
    dct = json.loads(json_str)
    grd_int_txt = dct['grid-intersect']
    out_path =  dct['output-dir']
    areas_path=dct['areas-dir']+'areas.csv'
    start = time()    # read the grid intersect (fua or nuts) 
    lat_area_dct= pd.read_csv(areas_path, header=None, index_col=0, skiprows=2).to_dict()[1]
    area=read_nuts_area(grd_int_txt, calcall=True)
    # levels (they are colled the same in fua and nuts)
    nuts_lvs= ['NUTS_Lv0', 'NUTS_Lv1','NUTS_Lv2', 'NUTS_Lv3']

    if sector_lst[-1] ==14 :
        dct_ms={'GNFR1': 'GNFR01',
         'GNFR10': 'GNFR10',
         'GNFR11': 'GNFR11',
         'GNFR12': 'GNFR12',
         'GNFR13': 'GNFR13',
         'GNFR2': 'GNFR02',
         'GNFR3': 'GNFR03',
         'GNFR4': 'GNFR04',
         'GNFR5': 'GNFR05',
         'GNFR6': 'GNFR06',
         'GNFR7': 'GNFR07',
         'GNFR8': 'GNFR08',
         'GNFR9': 'GNFR09',
         'ALL': 'GNFRALL'}
    elif sector_lst[-1] == 13 :
        dct_ms={'GNFR1': 'GNFR01',
         'GNFR10': 'GNFR10',
         'GNFR11': 'GNFR11',
         'GNFR12': 'GNFR12',
         'GNFR2': 'GNFR02',
         'GNFR3': 'GNFR03',
         'GNFR4': 'GNFR04',
         'GNFR5': 'GNFR05',
         'GNFR6': 'GNFR06',
         'GNFR7': 'GNFR07',
         'GNFR8': 'GNFR08',
         'GNFR9': 'GNFR09',
         'ALL': 'GNFRALL'}
    
    dct_units={'conc': '[\u03bcg/m\u00B3]', '[delta_concentration]':'[\u03bcg/m\u00B3]',
               'v_dll_pp':"[dll/(person year)]",'d_dll_pp':'[dll/(person year)]', 
               'v_mort': '[people/year]', 'd_mort': '[people/year]'}
    
    for nuts_lv in nuts_lvs:
        # prepare pd dataframe for results (index are the names of the 
        # geographical units in each level and keys are the delta and the base
        # case values)
        res = pd.DataFrame(index=area[nuts_lv]['area'].index.levels[0],
                           columns=['bc','delta'])
        # for the delta and base case value
        for key in ['bc','delta']:          
            # read the corresonding nc
            nc = read_nc(dct[key]['path'])   
            if dct[key]['aggregation']=="'sume'":
#                # @todo to be removed when path is updated from GUI
#                if key =='delta':
#                    pathdelta=dct[key]['path'][:-3]+'_all.nc'
#                    nc =  read_nc(pathdelta)                 
                t=ast.literal_eval(dct[key]['var'])    
                ms=dct_ms[t[1]] #'Nsnaps11'
                tpl_new=(t[0], ms)
                nct=pd.DataFrame(columns=[(tpl_new)])
                # assign to the dataframe the tansposed matrix of
                 # the corresponding value 
                nct[tpl_new]=nc.loc[tpl_new].transpose()
#                    print(nct.sum())

            else: 
            # drop level if necessary (e.g. when reading delta_concentration)
                tpl = dct[key]['var']
                if len(tpl)!=nc.index.nlevels:
                       nc = nc.reset_index(level=0, drop=True)
    
                # create a pd dataframe which has as columns the delta or base 
                # case value
                nct=pd.DataFrame(columns=[tpl])
                # assign to the dataframe the tansposed matrix of
                # the correscponding value 
                nct[tpl]=nc.loc[tpl].transpose()

            arealist = area[nuts_lv]['area'].index.levels[0]    
            aggr=ast.literal_eval(dct[key]['aggregation'])
            if len(aggr)==2:
                opt1=aggr[0]
                opt2=aggr[1]
            else:
                opt1=aggr
                opt2=None
                 
            for areait in arealist:
            # create a pd dataframe with the areas of each cell 
            # in the geographical unit
                # if the aggregation mode is average (e.g. concentration levels)
                if opt1 =='avg':
                    # area averaged values 
                    df_areas = pd.DataFrame(area[nuts_lv][opt2].loc[areait])
                    df_areas['mult'] = df_areas[opt2].multiply(nct.reindex(df_areas.index)[tpl])
                    areaxvar=df_areas[df_areas['mult'].notnull()]['mult'].sum()                                      
                    areatot=df_areas[df_areas['mult'].notnull()][opt2].sum()
                    # if the aggregation mode is sum (e.g. mortality)
                elif opt1 =='sum':
                    # sum the values in each cell considering the fraction of the 
                    # cell belonging to a geographical entity
                    df_areas = pd.DataFrame(area[nuts_lv]['parea'].loc[areait])
                    df_areas['mult'] = df_areas['parea'].multiply(nct.reindex(df_areas.index)[tpl])
                    areaxvar=df_areas[df_areas['mult'].notnull()]['mult'].sum()
                    areatot = 100 #20210201EP, change as percentage in the grid_intersect file are between 0 and 100, not between 0 and 1
                elif opt1 =='sume':
                    # sum the values in each cell considering the area of the 
                    # cell belonging to a geographical entity. 
                    df_areas = pd.DataFrame(area[nuts_lv][['area','parea']].loc[areait]) 
                   
                    df_areas['area_cell']=[area_calc(rw, lat_area_dct) 
                                           for rw in df_areas.index]

                    df_areas['mult'] = df_areas['parea'].multiply(df_areas['area_cell']).multiply(nct.reindex(df_areas.index)[tpl_new])
                    
                    
                    areaxvar=df_areas[df_areas['mult'].notnull()]['mult'].sum()
                    areatot = 1 #20220606, there was a bug
                # this if statement is to take care of areas in the 
                # shape file which are outside the domain 
                if areatot is not 0: 
                    value=areaxvar/areatot
                else:
                    value= float('NaN')
                res[key].loc[areait]=value           
        try:    
            units = dct_units[dct['bc']['var']]
        except: #KeyError
            units='[Mg]'
            
        print('Saving results')

#            print(dct['bc']['var'])
        if dct['bc']['var'] == 'v_mort' or dct['bc']['var'] == 'v_dll' or dct['bc']['var']=='v_dll_pp':
#                print('Quick and dirty fix of bug - see comments')
            # I (EPE) made a mistake - when reading values to aggregate for the interface, the 
            # label 'bc' actually refers to value... (the scenario), therefore 
            # I need to substitue only for this case the values before saving results. 
            # likle this I do not need to change anything in the interface. 
            # this should be fixed better in the future. 
            res['value']=res['bc']
            res['bc']=res['value']+res['delta']
        else: 
            res['value']=res['bc']-res['delta']

        res['per']=(res['delta'])/res['bc']*100
        
        res[['value', 'delta', 'per']].rename(
            columns={'value':'value'+units, 'delta':'delta'+units, 'per':'per[%]'}).to_csv(
                    out_path+nuts_lv[0:4]+nuts_lv[-1]+'.txt', header=True, index=True, sep='\t', na_rep='NaN', mode='w', 
                    encoding='utf-8')  
    end = time()
    print('Calculation time  for aggregation', end-start)
    

    
if __name__ == '__main__': 
    
    pass  
