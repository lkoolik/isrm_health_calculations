#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ISRM Data Object

@author: libbykoolik
last modified: 2022-03-23
"""

# Import Libraries
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file as nf
import os
from os import path
import sys

#%% Define the ISRM Object
class isrm:
    '''
    Defines a new object for storing and manipulating ISRM data.
    
    INPUTS:
        - file_path: file path to the ISRM grid NetCDF file (underlying data)
        - geo_file_path: file path to the shapefile of geographic information
        - load_file: set to True to import emissions, otherwise will just run checks
        - verbose: enable for more detailed outputs 
        
    '''
    def __init__(self, file_path, geo_file_path, load_file=True, verbose=False):
        ''' Initializes the ISRM object'''        
        # Initialize paths and check that they are valid
        self.file_path = file_path
        self.geo_file_path = geo_file_path
        self.valid_file, self.valid_geo_file = self.check_path()
        
        # Grab other meta-parameters
        self.load_file = load_file
        self.verbose = verbose
        verboseprint = print if self.verbose else lambda *a, **k:None # for logging
        verboseprint('\nCreating a new ISRM object from {}'.format(self.file_path))
        
        # If the files do not exist, quit before opening
        if not self.valid_file:
            print('\n<< ERROR: The filepath provided for the ISRM netCDF is not correct. Please correct and retry. >>')
            sys.exit()
        elif not self.valid_geo_file:
            print('\n<< ERROR: The filepath provided for the ISRM boundaries is not correct. Please correct and retry. >>')
            sys.exit()
        else:
            verboseprint('- Filepaths and files found. Proceeding to import ISRM data.')
        
        # Read ISRM data and geographic information
        if self.valid_file == True and self.load_file == True and self.valid_geo_file == True:
            verboseprint('- Beginning to import ISRM netCDF data. This step may take some time.')
            self.PM25, self.NH3, self.NOX, self.VOC, self.SOX = self.load_isrm()
            verboseprint('- ISRM netCDF data imported. Five pollutant variables created')
            verboseprint('- Beginning to import ISRM geographic data. This step may take some time.')
            self.geodata = self.load_geodata()
            verboseprint('- ISRM geographic data imported.')
            self.crs = self.geodata.crs
            self.geometry = self.geodata['geometry']
            self.ISRM_ID = self.geodata['ISRM_ID']
    
    def __str__(self):
        return 'ISRM object created from '+self.file_path

    def __repr__(self):
        return '< ISRM object created from '+self.file_path+'>'

    
    def check_path(self):
        ''' Checks if file exists at the path specified '''
        # Use the os library to check the path and the file
        # First, check ISRM netCDF exists
        path_exists = path.exists(self.file_path)
        file_exists = path.isfile(self.file_path)
        
        # Second, check ISRM geodata exists
        geo_path_exists = path.exists(self.geo_file_path)
        geo_file_exists = path.isfile(self.geo_file_path)
            
        return (path_exists and file_exists, geo_path_exists and geo_file_exists)
    
    def load_isrm(self):
        ''' Loads ISRM File '''
        # Open the NetCDF file 
        isrm_nf = nf(self.file_path, mode='r', mmap=False)
        
        # Define the Pollutant Layers and Loop through to Slice into Variables
        pollutant_layers = ['PrimaryPM25', 'pNH4', 'pNO3', 'pSO4', 'SOA']
        pollutants = []
        
        for p in pollutant_layers:
            pollutants.append(self.split_pollutant(isrm_nf, p))    
        
        # Close the NetCDF File before exiting
        isrm_nf.close()
        
        return pollutants
    
    def split_pollutant(self, isrm_nf, pol_name):
        ''' Grabs relevant pollutant layer data '''
        pol_data = isrm_nf.variables[pol_name].data.copy()
        
        return pol_data
    
    def load_geodata(self):
        ''' Loads shapefile into geopandas dataframe '''
        isrm_gdf = gpd.read_file(self.geo_file_path)
        isrm_gdf.columns = ['ISRM_ID', 'geometry']
        
        return isrm_gdf
    
    def get_pollutant_layer(self, pol_name):
        ''' Returns pollutant layer '''
        # Define a pollutant dictionary for convenience
        pollutant_dict = {'PM25':self.PM25,
                         'NH3':self.NH3,
                         'VOC':self.VOC,
                         'NOX':self.NOX,
                         'SOX':self.SOX}
        
        # Confirm pol_name is valid
        assert pol_name in pollutant_dict.keys()
        
        # Return pollutant layer
        return pollutant_dict[pol_name]
    
    def map_isrm(self):
        ''' Creates map of ISRM grid  '''
        # Note to build this out further at some point in the future, works for now
        fig, ax = plt.subplots(1,1)
        self.geodata.plot(ax = ax, edgecolor='black', facecolor='none')
        ax.set_title('ISRM Grid')
        fig.tight_layout()
        return fig