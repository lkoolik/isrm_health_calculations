#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Population Data Object

@author: libbykoolik
last modified: 2022-08-10
"""

# Import Libraries
import pandas as pd
import geopandas as gpd
import logging
import numpy as np
import matplotlib.pyplot as plt
import pyarrow
from scipy.io import netcdf_file as nf
import os
from os import path
import sys

#%% Define the Population Object
class population:
    '''
    Defines a new object for storing and manipulating concentration data.
    
    INPUTS:
        - file_path: the file path of the raw population data
        - load_file: a Boolean indicating whether or not the file should be loaded 
        - verbose: a Boolean indicating whether or not detailed logging statements 
          should be printed
          
    CALCULATES:
        - valid_file: a Boolean indicating whether or not the file provided is valid
        - geometry: geospatial information associated with the emissions input
        - pop_all: complete, detailed population data from the source
        - pop_geo: a geodataframe with population IDs and spatial information
        - crs: the inherent coordinate reference system associated with the emissions input
        - pop_exp: a geodataframe containing the population information with associated 
          spatial information, summarized across age bins
        - pop_hia: a geodataframe containing the population information with associated
          spatial information, broken out by age bin
          
    EXTERNAL FUNCTIONS:
        - allocate_population: reallocates population into new geometry using a 
          spatial intersect
        
    '''
    def __init__(self, file_path, load_file=True, verbose=False):
        ''' Initializes the Population object'''        
        logging.info('\n << Reading Population Census Data for EJ Exposure Calculations >>')
        
        # Gather meta data
        self.file_path = file_path
        self.file_type = file_path.split('.')[-1].lower()
        self.load_file = load_file
        self.verbose = verbose
        verboseprint = logging.info if self.verbose else lambda *a, **k:None # for logging
        verboseprint('- Creating a new Census Tract population object from {}'.format(self.file_path))
        
        # Initialize population object by reading in the feather file
        self.valid_file = self.check_path()
        
        if not self.valid_file:
            logging.info('\n << ERROR: The filepath provided is not correct. Please correct and retry. >>')
            sys.exit()
        
        # Read in the data
        if self.load_file == True and self.valid_file:
            verboseprint('- Attempting to load the population data. This step may take some time.')
            self.pop_all, self.pop_geo, self.crs = self.load_population()
            self.pop_exp, self.pop_hia = self.split_pop_objs()
            verboseprint('- Population successfully loaded.')        
            
    def __str__(self):
        return '< Population object for year '+str(self.year)+ '>'

    def __repr__(self):
        return '< Emissions object created from '+self.file_path + '>'

    def check_path(self):
        ''' Checks if file exists at the path specified '''
        # Use the os library to check the path and the file
        path_exists = path.exists(self.file_path)
        file_exists = path.isfile(self.file_path)
        return path_exists and file_exists
    
    def load_population(self):
        ''' Loads the population file, depending on the extension ''' 
        # Based on the file extension, run different load functions
        if self.file_type == 'shp':
            pop_all = self.load_shp()
        
        if self.file_type == 'feather':
            pop_all = self.load_feather()
            
        # Create a variable that is just geometry and IDs
        pop_geo = pop_all[['POP_ID','geometry']].copy().drop_duplicates()
        pop_crs = pop_geo.crs
        
        return pop_all, pop_geo, pop_crs
    
    def load_shp(self):
        ''' Loads population data from a shapefile. '''
        # Shapefiles are read using geopandas
        pop_all = gpd.read_file(self.file_path)
        
        return pop_all
                
    def load_feather(self):
        ''' Loads population data from a feather file. '''
        # Feather file is read using geopandas
        pop_all = gpd.read_feather(self.file_path)
        
        return pop_all

    def split_pop_objs(self):
        ''' '''
        # Create a copy of the population data to avoid overwriting
        pop_tmp = self.pop_all.copy()
        
        ## Create the exposure calculation population object
        # For the exposure calculations, we do not need the age bins
        pop_exp = pop_tmp[['POP_ID', 'YEAR', 'TOTAL', 'ASIAN', 'BLACK', 'HISLA', 
                           'INDIG', 'PACIS', 'WHITE', 'OTHER']].copy()
        
        # Sum across POP_ID and YEAR
        pop_exp = pop_exp.groupby(['POP_ID','YEAR'])[['TOTAL', 'ASIAN', 'BLACK', 
                                                      'HISLA', 'INDIG', 'PACIS', 
                                                      'WHITE', 'OTHER']].sum().reset_index()
        
        # Add geometry back in
        pop_exp = pd.merge(self.pop_geo, pop_exp, on='POP_ID')
        
        ## Create the health impact calculation population object
        # TBD
        pop_hia=...
        
        return pop_exp, pop_hia

    def project_pop(self, pop_obj, new_crs):
        ''' Projects the population data into a new crs '''
        pop_obj_prj = pop_obj.to_crs(new_crs)
    
        return pop_obj_prj
    
    def allocate_population(self, new_geometry, new_geometry_ID):
        ''' Reallocates the population into the new geometry using a spatial intersect '''
        if self.verbose:
            logging.info('- Allocating population from Census tracts to ISRM grid cells.')
        
        # Confirm that the coordinate reference systems match
        #assert pop_tmp.crs == new_geometry.crs, 'Coordinate reference system does not match. Population cannot be reallocated'
        if self.crs == new_geometry.crs:
            pop_tmp = self.pop_exp.copy(deep=True)
        else:
            pop_tmp = self.project_pop(self.pop_exp, new_geometry.crs)
        
        # Add the land area as a feature of this dataframe
        pop_tmp['AREA_M2'] = pop_tmp.geometry.area/(1000.*1000.)
        
        
        # Create intersect object
        intersect = gpd.overlay(pop_tmp, new_geometry, how='intersection')
        pop_totalarea = intersect.groupby('POP_ID').sum()['AREA_M2'].to_dict()
    
        # Add a total area and area fraction to the intersect object
        intersect['AREA_POP_TOTAL'] = intersect['POP_ID'].map(pop_totalarea)
        intersect['AREA_FRAC'] = intersect['AREA_M2'] / intersect['AREA_POP_TOTAL']

        # Define the racial/ethnic groups and estimate the intersection population
        cols = ['TOTAL', 'ASIAN', 'BLACK', 'HISLA', 'INDIG', 'PACIS', 'WHITE','OTHER']
        for c in cols:
            intersect[c] = intersect[c] * intersect['AREA_FRAC']
            
        # Sum across new geometry grid cells
        new_alloc_pop = intersect.groupby([new_geometry_ID])[cols].sum().reset_index()
        
        # Merge back into the new geometry using the new_geometry_ID
        new_alloc_pop = new_geometry.merge(new_alloc_pop, how='left',
                                           left_on=new_geometry_ID,
                                           right_on=new_geometry_ID)

        # Fill the missing cells with zero population
        new_alloc_pop[cols] = new_alloc_pop[cols].fillna(0)
        
        # Confirm that the population slipt was close
        old_pop_total = pop_tmp[cols].sum()
        new_pop_total = pop_tmp[cols].sum()
        
        for c in cols:
            assert np.isclose(old_pop_total[c], new_pop_total[c])
        
        # Print statement
        if self.verbose:
            logging.info('- Census tract population data successfully re-allocated to the ISRM grid.')
        
        return new_alloc_pop
    
test_pop = '/Users/libbykoolik/Documents/Research/OEHHA Project/scripts/isrm_health_calculations/data/ca2010.feather'
test = population(test_pop)