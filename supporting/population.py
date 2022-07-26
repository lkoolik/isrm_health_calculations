#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Population Data Object

@author: libbykoolik
last modified: 2022-07-19
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
        - tbd
        
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
            self.geometry, self.pop_data, self.crs, self.pop_gdf = self.load_population()
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
        if self.file_type == 'feather':
            geometry, pop_data, crs, pop_gdf = self.load_feather()
        
        return geometry, pop_data, crs, pop_gdf
    
                
    def load_feather(self):
        ''' 
        Loads population data from a feather file.
        
        Requirements:
            - TBD
        
        '''
        # Feather file is read using geopandas
        pop_gdf = gpd.read_feather(self.file_path)
        pop_gdf['POP_ID'] = 'POP_'+pop_gdf.index.astype(str)
        pop_gdf.columns = pop_gdf.columns.str.upper()
        pop_gdf.rename(columns={'GEOMETRY':'geometry'}, inplace=True)
        
        # Split off geometry from emissions data
        geometry = pop_gdf[['POP_ID','geometry']].copy()
        pop_data = pd.DataFrame(pop_gdf.drop(columns='geometry'))
        
        # Separately save the coordinate reference system
        crs = pop_gdf.crs
        
        return geometry, pop_data, crs, pop_gdf

    
    def project_pop(self, new_crs):
        ''' Projects the population data into a new crs '''
        pop_gdf_prj = self.pop_gdf.to_crs(new_crs)
    
        return pop_gdf_prj
    
    def allocate_population(self, new_geometry, new_geometry_ID):
        ''' Reallocates the population into the new geometry using a spatial intersect '''
        if self.verbose:
            logging.info('- Allocating population from Census tracts to ISRM grid cells.')
        
        # Confirm that the coordinate reference systems match
        #assert pop_tmp.crs == new_geometry.crs, 'Coordinate reference system does not match. Population cannot be reallocated'
        if self.crs == new_geometry.crs:
            pop_tmp = self.pop_gdf.copy(deep=True)
        else:
            pop_tmp = self.project_pop(new_geometry.crs)
        
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