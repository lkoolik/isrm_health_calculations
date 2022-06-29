#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Health Impact Function Meta Data Object

@author: libbykoolik
last modified: 2022-06-29
"""

# Import Libraries
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import os
from os import path
import sys
sys.path.append('/Users/libbykoolik/Documents/Research/OEHHA Project/scripts/isrm_health_calculations/supporting')

#%% Define the Health Data Object
class health_data:
    '''
    Defines a new object for storing and manipulating health impact input data.
    
    INPUTS:
        - filepath_dict: a dictionary with the filepaths of each input feather.
          Should have the following keys:
              (1) POPULATION
              (2) INCIDENCE
              (3) POPULATION-INCIDENCE
        - verbose: a Boolean enabling more detailed output statements
        
    '''
    def __init__(self, filepath_dict, verbose):
        ''' Initializes the Health Input object'''        
        # Get object metadata
        self.filepath_dict = filepath_dict
        self.verbose = verbose
        verboseprint = print if self.verbose else lambda *a, **k:None # for logging
        verboseprint('\nDownloading the input data for calculating excess mortality.')

        # Add source information (hard-coded for now)
        self.population_source = 'US 2010 Census County Level from BenMAP CE'
        self.incidence_source = 'All-Cause Mortality Incidence (2010) from BenMAP CE'
        
        # Initialize object by loading the health data
        self.population, self.incidence, self.pop_inc = self.load_data()
        verboseprint('- Population data source {} imported.'.format(self.population_source))
        verboseprint('- Incidence data source {} imported.'.format(self.incidence_source))
        
            
    def __str__(self):
        return 'Health impact function input population from '+self.population_source+' and incidence from '+self.incidence_source

    def __repr__(self):
        return '< Health impact object.>'
    
    def load_data(self):
        ''' Loads population and incidence data from feather files. '''
        data_dict = {}
        
        # Get filepaths from the filepath dictionary
        population = gpd.read_feather(self.filepath_dict['POPULATION'])
        incidence = gpd.read_feather(self.filepath_dict['INCIDENCE'])
        pop_inc = gpd.read_feather(self.filepath_dict['POPULATION-INCIDENCE'])            
        
        return population, incidence, pop_inc
    