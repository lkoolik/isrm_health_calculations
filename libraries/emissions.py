#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Emissions Data Object

@author: libbykoolik
last modified: 2022-03-22
"""

# Import Libraries
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import os
from os import path
import sys

#%% Define the Emissions Object
class emissions:
    '''
    Defines a new object for storing and manipulating emissions data.
    
    INPUTS:
        - file_path: File path to an emissions file (currently allows .shp only)
        - units: Units of the emissions file (default is ug/s)
        - name: plain English tag for emissions, otherwise will use the filename
        - details_to_keep: any additional aggregation field (e.g., FUEL_TYPE)
        - load_file: set to True to import emissions, otherwise will just run checks
        - verbose: enable for more detailed outputs
        
    '''
    def __init__(self, file_path, units='ug/s', name='', details_to_keep=[], load_file=True, verbose=False):
        ''' Initializes the emissions object'''        
        # Initialize path and check that it is valid
        self.file_path = file_path
        self.file_type = file_path.split('.')[-1].lower()
        self.valid_file = self.check_path()

        # Define additional Meta Variables
        self.emissions_name = self.get_name(name)
        self.details_to_keep = details_to_keep
        self.verbose = verbose
        verboseprint = print if self.verbose else lambda *a, **k:None # for logging
        verboseprint('Creating a new emissions object from {}'.format(self.file_path))
        
        # If the file does not exist, quit before opening
        if not self.valid_file:
            print('\n<< ERROR: The filepath provided is not correct. Please correct and retry. >>')
            sys.exit()
        else:
            verboseprint('- Filepath and file found. Proceeding to import emissions data.')
        
        # Confirm Units are Valid
        self.units = units 
        self.valid_units = self.check_units()
        if not self.valid_units:
            print('\n<< ERROR: The units provided ({}) are not valid. Please convert emissions to ug/s (micrograms per second) and try again. >>'.format(self.units))
            sys.exit()
        else:
            verboseprint('- Units of provided emissions ({}) are valid.'.format(self.units))
            
        # If load_file is True, import the emissions data
        if load_file == True and self.valid_file:
            verboseprint('- Attempting to load the emissions data. This step may take some time.')
            self.geometry, self.emissions_data, self.crs = self.load_emissions()
            verboseprint('- Emissions successfully loaded.')
            self.emissions_data.columns = map(str.upper, self.emissions_data.columns)
            self.valid_emissions = self.check_emissions()
            
            if self.valid_emissions:
                verboseprint('- Emissions formatting has been checked and confirmed to be valid.')
                verboseprint('- Beginning data cleaning and processing.')
                # Should add a statement about details_to_keep and the aggfunc
                self.emissions_data_clean = self.clean_up(self.details_to_keep)
                self.PM25 = self.split_pollutants(self.emissions_data_clean, 'PM25', self.details_to_keep)
                self.NH3 = self.split_pollutants(self.emissions_data_clean, 'NH3', self.details_to_keep)
                self.NOX = self.split_pollutants(self.emissions_data_clean, 'NOX', self.details_to_keep)
                self.VOC = self.split_pollutants(self.emissions_data_clean, 'VOC', self.details_to_keep)
                self.SOX = self.split_pollutants(self.emissions_data_clean, 'SOX', self.details_to_keep)
    
    def __str__(self):
        return 'Emissions object created from '+self.file_path

    def __repr__(self):
        return '< Emissions object created from '+self.file_path


    def get_file_path(self):
        ''' Returns the file path '''
        print('Emissions file path: ', self.file_path)
        return self.file_path
    
    def get_name(self, name):
        ''' Gets emissions file name from the file_path unless provided '''
        if name == '':
            emissions_name = path.basename(self.file_path).split('.')[0]
        else:
            emissions_name = name
        return emissions_name
    
    def get_unit_conversions(self):
        ''' Hard-coded dictionary of unit conversions '''
        # The purpose of this function is to create a singular place where units have
        # to be updated if future iterations add more units.
        
        # mass_units and time_units are dictionaries with keys = current units and values = conversion factor to ug or s
        mass_units = {'ug':1.0, 'g':10.0**6, 'lb':453592370.0, 'ton':907184740000.0, 
                      'mt':1000000000000.0, 'kg':10.0**9} # micrograms per unit
        
        time_units = {'s': 1.0, 'min':60.0, 'hr':3600.0, 'day':86400.0, 'yr':31536000.0}
        
        return mass_units, time_units
    
    def check_path(self):
        ''' Checks if file exists at the path specified '''
        # Use the os library to check the path and the file
        path_exists = path.exists(self.file_path)
        file_exists = path.isfile(self.file_path)
        return path_exists and file_exists

    def check_units(self):
        ''' Checks if units are compatible with program '''
        units = self.units.lower() # force lower for simplicity
        
        # Make this somewhat flexible for more possibilities
        unit_mass, unit_time = units.split('/')
        mass_units, time_units = self.get_unit_conversions()
        
        # Check if units are in the dictionaries
        valid_units = (unit_mass in mass_units.keys()) and (unit_time in time_units.keys()) 
        
        return valid_units
    
    def load_emissions(self):
        ''' Loads the emissions file, depending on the extension ''' 
        # Based on the file extension, run different load functions
        if self.file_type == 'shp':
            geometry, emissions_data, crs = self.load_shp()
        
        return geometry, emissions_data, crs

                
    def load_shp(self):
        ''' 
        Loads emissions data from a shapefile.
        
        Requirements:
            - Emissions data must have the columns I_CELL and J_CELL that should be uniquely 
            indexed to spatial information
        
        '''
        # Shapefiles are read using geopandas
        emissions_gdf = gpd.read_file(self.file_path)
        
        # Split off geometry from emissions data
        geometry = emissions_gdf[['I_CELL','J_CELL','geometry']]\
            .copy().drop_duplicates().reset_index(drop=True)
        emissions_data = pd.DataFrame(emissions_gdf.drop(columns='geometry'))
        
        # Separately save the coordinate reference system
        crs = emissions_gdf.crs
        
        return geometry, emissions_data, crs

    
    def check_emissions(self):
        ''' Runs a number of checks to make sure that emissions data is valid  '''
        ## (TEST 1) Check for I_CELL and J_CELL
        has_indices = ('I_CELL' in self.emissions_data.columns) & ('J_CELL' in self.emissions_data.columns)
        
        ## (TEST 2) Check that columns are correct, or replace them if incorrect
        # Define the set of correct pollutants
        correct_pollutants = ['PM25', 'NH3', 'VOC', 'NOX', 'SOX']
        
        # Define a few dummy variables for storing information
        missing = []
        
        # Iterate through the list of correct pollutants to find and address issues
        for pol in correct_pollutants:
            if pol in self.emissions_data.columns: # True if correct pollutant is in the emissions_data
                pass # No action necessary
            else: # Enter the missing pollutant decision tree if missing
                choice=self.map_pollutant_name(pol) # Returns the missing pollutant if there is no suitable replacement
                if choice != pol:
                    self.emissions_data.rename(columns={choice:pol}, inplace=True) # Updates the emissions data
                else:
                    missing.append(pol) # Do not spit out error and exit yet
                    
        # After iterating through all of the pollutants, check which couldn't be identified or mapped
        if len(missing)>0:
            print('\n<< ERROR: required pollutants are missing from data. Please check data and try again. >>')
            print('- The tool requires that these five pollutants be included as columns in emissions data:', ', '.join(correct_pollutants))
            print('- The tool was not able to find or map:', ', '.join(missing))
            sys.exit()
            
        else: # All five pollutants are accounted for
            has_five_pollutants = True

        return has_indices and has_five_pollutants
    
    def map_pollutant_name(self, missing_pol):
        ''' If a pollutant name is not found in the emissions data, tries to map it before quitting '''
        # Below, there are a number of potential erroneous column names for each pollutant
        # These are ordered for most appropriate
        possible_wrong_names_dict = {'PM25':['PM2.5', 'PM2_5'],
                                     'NH3':['NH_3'],
                                     'VOC':['ROG', 'TOG', 'OG'],
                                     'NOX':['NO_X', 'NO2', 'NO_2'],
                                     'SOX':['SO_X', 'SO2', 'SO_2']}
        
        # Cut the dictionary into just the relevant ones for missing_pol
        possible_wrong_names = possible_wrong_names_dict[missing_pol]
        first_choice = 'na'
        
        # Iterate through the list of potential wrong names, stop if one is found
        for choice in possible_wrong_names:
            if choice in self.emissions_data.columns:
                print('\n<< Warning: {} was not found in the emissions data. The program will continue using {} as a replacement. >> \n'.format(missing_pol, choice))
                first_choice = choice
                break
        
        # If nothing is found, return just the name of the missing pollutant
        if first_choice == 'na':
            return missing_pol
        
        # Otherwise, return the first choice
        else:
            return first_choice
    
    def clean_up(self, details_to_keep, func=np.sum):
        ''' Simplifies emissions file by reducing unnecessary details '''
        # Start by deep-copying to avoid any accidental overwriting
        emissions_data_tmp = self.emissions_data.copy(deep=True)
        
        # Use pandas processing to perform a groupby
        groupby_features = ['I_CELL', 'J_CELL']+details_to_keep
        emissions_data_clean = emissions_data_tmp.groupby(groupby_features).agg(func)
        if self.verbose:
            print('- Emissions have been reduced to contain the {} of emissions for each {}'.format(func.__name__, ', '.join(groupby_features)))
        
        # Clean up indices
        emissions_data_clean = emissions_data_clean.reset_index()
        
        # Scale emissions to ug/s
        scaling_factor = self.convert_units()
        emissions_data_clean[['PM25', 'NH3', 'VOC', 'NOX', 'SOX']] *= scaling_factor
        if scaling_factor != 1.0 and self.verbose:
            print('- Scaled emissions data by a factor of {:e} to convert from {} to ug/s.'.format(scaling_factor, self.units))
        
        # Limit columns
        emissions_data_clean = emissions_data_clean[groupby_features+['PM25', 'NH3', 'VOC', 'NOX', 'SOX']]
        
        return emissions_data_clean
    
    def convert_units(self):
        ''' Based on provided units, converts all emissions columns to ug/s '''
        # Import the units and the unit conversions
        units = self.units.lower() # force lower for simplicity
        unit_mass, unit_time = units.split('/')
        mass_units, time_units = self.get_unit_conversions()
        
        # Determine the scaling factor based on the provided units
        scaling_factor = mass_units[unit_mass]/time_units[unit_time]
        
        return scaling_factor
        
    
    def split_pollutants(self, emissions_object, pollutant, details_to_keep):
        ''' Creates separate objects for a pollutant '''
        # Start by deep-copying to avoid any accidental overwriting
        emissions_data_tmp = emissions_object.copy(deep=True)
        
        # Pull only the I_CELL, J_CELL, pollutant emissions, and other details
        groupby_features = ['I_CELL','J_CELL',pollutant]+details_to_keep
        pollutant_emissions = emissions_data_tmp[groupby_features].copy().drop_duplicates().reset_index()
        
        # Add geometry back in
        pollutant_emissions = self.geometry.merge(pollutant_emissions, left_on=['I_CELL','J_CELL'], right_on=['I_CELL','J_CELL'])
        
        # Rename emissions column to 'Emissions_ug/s' for easier functionality later
        pollutant_emissions.rename(columns={pollutant:'EMISSIONS_UG/S'}, inplace=True)
        pollutant_emissions.drop(columns=['index'], inplace=True)
        
        if self.verbose:
            print('- Successfully created emissions object for {} for {}.'.format(self.emissions_name, pollutant))
        
        return pollutant_emissions
    
    def visualize_emissions(self, emissions_object, pollutant_name=''):
        ''' Creates map of emissions using simple chloropleth '''
        # Note to build this out further at some point in the future, works for now
        fig, ax = plt.subplots(1,1)
        emissions_object.plot(column='EMISSIONS_UG/S',
                              figsize=(20,10),
                              legend=True,
                              legend_kwds={'label':r'Emissions ($\mu$g/s)'},
                              ax = ax)
        ax.set_title(pollutant_name+' Emissions')
        fig.tight_layout()
        return fig

    # def reallocate_emissions(self, emissions_object, new_geometry):
    #     ''' 
    #     Reallocates single-pollutant emissions spatially based on new_geometry.
        
    #     INPUTS:
    #         - emissions_object: assumed to be a geopandas dataframe with columns
    #                             ['I_CELL', 'J_CELL', 'geometry', 'EMISSIONS_UG/S']
    #         - new_geometry: the new geometry to reallocate to (geopandas dataframe)
    #     '''
    #     # Create a copy of the emissions object and assign temporary IDs
    #     emissions = emissions_object.copy(deep=True)
        
    #     # Assign temporary IDs to emissions and new_geometry
    #     emissions['EMISSIONS_ID'] = emissions.apply(lambda x: str(x['I_CELL'])+'_'+str(x['J_CELL']), axis=1)
    #     new_geometry['NEW_GEO_ID'] = new_geometry.apply(lambda x: 'NEW_GEO_ID_'+str(x.index))
        
    #     # Project this new emissions object to the geometry of the new_geometry object
    #     crs_use = new_geometry.crs
    #     emissions = emissions.to_crs(crs_use)
        
    #     # Get total area of emissions cell (should be ~1 km2)
    #     emissions['AREA_KM2'] = emissions.geometry.area/(1000*1000)
        
    #     # Create intersect between emis and grid
    #     intersect = gpd.overlay(emissions, new_geometry, how='intersection')
    #     emissions_totalarea = intersect.groupby('EMISSIONS_ID').sum()['AREA_KM2'].to_dict()
        
    #     # Add a Total Area and Area Fraction 
    #     intersect['AREA_TOTAL_KM2'] = intersect['EMISSIONS_ID'].map(emissions_totalarea)
    #     intersect['AREA_FRAC_KM2'] = intersect['AREA_KM2'] / intersect['AREA_TOTAL_KM2']
        
    #     # Allocate emissions based on area fraction
    #     intersect['EMISSIONS_UG/S_ALLOC'] = intersect['AREA_FRAC_KM2'] * intersect['EMISSIONS_UG/S']  
            
        
    #     # Sum over new geometry shape
    #     reallocated_emissions = intersect.groupby('NEW_GEO_ID')[['EMISSIONS_UG/S_ALLOC']].sum().reset_index()
    #     reallocated_emissions = new_geometry[['NEW_GEO_ID','geometry']].merge(reallocated_emissions,
    #                                                       left_on='NEW_GEO_ID',
    #                                                       right_on='NEW_GEO_ID')
        
    #     # Rename the column
    #     reallocated_emissions.rename(columns={'EMISSIONS_UG/S_ALLOC':'EMISSIONS_UG/S'}, inplace=True)
        
    #     return reallocated_emissions

#file_path = '/Users/libbykoolik/Documents/Research/OEHHA Project/data/emissions (dl 2022-03-07)/2000_annual_oehha_MPOv10.shp'
file_path = '/Users/libbykoolik/Documents/Research/OEHHA Project/data/test_emissions_data/test_data.shp'
tmp = emissions(file_path, load_file=True, verbose=True, units='lb/yr')
tmp.visualize_emissions(tmp.NH3, 'NH3')