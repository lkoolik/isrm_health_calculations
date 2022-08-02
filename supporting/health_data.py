#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Health Impact Function Meta Data Object

@author: libbykoolik
last modified: 2022-08-01
"""

# Import Libraries
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import logging
import os
from os import path
import sys

#%% Define the Health Data Object
class health_data:
    '''
    Defines a new object for storing and manipulating health impact input data.
    
    INPUTS:
        - filepath_dict: a dictionary with the filepaths of each input feather.
          Should have the following keys:
              (1) POPULATION
              (2) INCIDENCE
        - verbose: a Boolean enabling more detailed output statements
        - race_stratified: Boolean indicating whether race-stratified incidence
          rates should be used
          
    CALCULATES:
        - population: a geodataframe containing the raw population data from BenMAP
        - incidence: a geodataframe containing the raw incidence data from BenMAP
        - pop_inc: a geodataframe containing the combined population and incidence 
          data based on the requested geographies
        
    '''
    def __init__(self, filepath_dict, verbose, race_stratified):
        ''' Initializes the Health Input object'''   
        logging.info('\n << Loading BenMAP Health Inputs >>')
        # Get object metadata
        self.filepath_dict = filepath_dict
        self.verbose = verbose
        verboseprint = logging.info if self.verbose else lambda *a, **k:None # for logging
        verboseprint('- Downloading the input data for calculating excess mortality.')
        self.race_stratified = race_stratified

        # Add source information (hard-coded for now)
        self.population_source = 'US 2010 Census County Level from BenMAP CE reallocated to the ISRM grid'
        if not self.race_stratified:
            self.incidence_source = 'Mortality Incidence (2010) from BenMAP CE'
        else:
            self.incidence_source = 'Race-Strafified Mortality Incidence (2007-16) from BenMAP CE'
        
        # Initialize object by loading the health data
        self.population, self.incidence = self.load_data()
        verboseprint('- Population data (source: {}) successfully imported.'.format(self.population_source))
        verboseprint('- Incidence data (source: {}) successfully imported.'.format(self.incidence_source))
        
        # Combine the population and incidence data into one dataframe
        verboseprint('- Combining population and incidence data for health calculations. This step may take some time.')
        self.pop_inc, self.tmp = self.combine_pop_inc()
        verboseprint('- Population and incidence data succesfully combined and ready for health calculations.')
            
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
        
        return population, incidence
    
   
    def update_pop(self, population):
        ''' Performs a few population dataset updates before combining '''
        # Un-pivot the population data to have separate columns for RACE and POPULATION
        population = population.melt(id_vars=['ROW','START_AGE', 'END_AGE','geometry'], 
                                     value_vars=['ASIAN','BLACK','HISLA','INDIG',
                                                 'WHITE','TOTAL'], 
                                     var_name='RACE', value_name='POPULATION', 
                                     ignore_index=False)
        
        population.rename(columns={'ROW':'ISRM_ID'}, inplace=True)
        
        return population
    
    def update_inc(self, incidence):
        ''' Performs a few incidence dataset updates before combining '''
        
        inc_geo = incidence[['NAME','STATE_NAME', 'geometry']].drop_duplicates()
        
        # Pivot the incidence data around the endpoints to make one column for each
        incidence = incidence.pivot(index=['NAME', 'STATE_NAME', 'RACE', 'START_AGE', 'END_AGE'],
                                    columns='ENDPOINT', values='INCIDENCE').reset_index()
        
        # Clean up columns
        incidence.columns = incidence.columns.str.upper() # Capitalize
        incidence = incidence[['NAME', 'STATE_NAME', 'RACE', 'START_AGE', 'END_AGE', 
                                'ALL CAUSE','ISCHEMIC HEART DISEASE','LUNG CANCER']]
        
        incidence = pd.merge(inc_geo, incidence, left_on=['NAME','STATE_NAME'],
                             right_on=['NAME','STATE_NAME'])
        
        return incidence

    def get_incidence_lookup(self, incidence, name, race, start_age, endpoint):
        ''' Creates a small incidence lookup table based on the name and age ranges '''
        inc_small = incidence[(incidence['NAME']==name)]#&(incidence['RACE']==race)]
        inc_small = inc_small[(inc_small['END_AGE']>=start_age)&(inc_small['START_AGE']<=start_age)]
        
        return inc_small[endpoint].values[0]
    
    def get_incidence_pop(self, name, race, start_age, end_age, endpoint, incidence_lookup):
        ''' Helper function for returning incidence given name, race, age range, and endpoint '''
        n = end_age - start_age + 1
        
        snip = incidence_lookup[incidence_lookup['NAME']==name]
        snip = snip[(snip['START_AGE']>=start_age)&(snip['START_AGE']<=end_age)]
        
        inc = snip[endpoint].sum()/n
        
        return inc
    
    def make_incidence_lookup(self, incidence, keys, lookups):
        ''' Creates a dictionary for faster lookup table processing '''
        # Create a temporary dataframe for building the lookup table
        lookup_table_builder = keys.copy()
        
        # Add the fields required for building the lookup table
        lookup_table_builder['START_AGE'] = lookup_table_builder['KEY'].str.split('_').str[0].astype(int)
        lookup_table_builder['END_AGE'] = lookup_table_builder['KEY'].str.split('_').str[1].astype(int)
        lookup_table_builder['NAME'] = lookup_table_builder['KEY'].str.split('_').str[2]
        lookup_table_builder['RACE'] = lookup_table_builder['KEY'].str.split('_').str[3]
        
        # Use the get_incidence_pop function for estimating incidence for each key
        lookup_table_builder['ALL CAUSE'] = lookup_table_builder.apply(lambda x: self.get_incidence_pop(x['NAME'], x['RACE'], x['START_AGE'], x['END_AGE'], 'ALL CAUSE', lookups), axis=1)
        lookup_table_builder['ISCHEMIC HEART DISEASE'] = lookup_table_builder.apply(lambda x: self.get_incidence_pop(x['NAME'], x['RACE'], x['START_AGE'], x['END_AGE'], 'ISCHEMIC HEART DISEASE', lookups), axis=1)
        lookup_table_builder['LUNG CANCER'] = lookup_table_builder.apply(lambda x: self.get_incidence_pop(x['NAME'], x['RACE'], x['START_AGE'], x['END_AGE'], 'LUNG CANCER', lookups), axis=1)

        # Convert into a dictionary
        lookup_table_builder['INCIDENCE'] = lookup_table_builder[['ALL CAUSE','ISCHEMIC HEART DISEASE', 'LUNG CANCER']].to_numpy().tolist()
        lookup_dict = lookup_table_builder.set_index('KEY')['INCIDENCE'].to_dict()
        
        return lookup_dict
    
    def incidence_by_age(self, incidence, population):
        ''' Create a smaller incidence lookup tables for merging '''
        # Get an array of the start ages that increment by 1 only
        start_ages = pd.DataFrame({'START_AGE':range(population['END_AGE'].max())})
        
        # Get a dataframe of the 'NAME' and 'RACE' combination *** TO ADD RACE
        incidence_names = incidence[['NAME']].drop_duplicates()
        
        # Add a dummy column for performing an outer join
        start_ages['key'] = 'tmp'
        incidence_names['key'] = 'tmp'        
        
        # Create a lookup table
        incidence_lookup = pd.merge(incidence_names, start_ages, on='key', how='outer')
        incidence_lookup['RACE'] = 'ALL'
        
        # Add each endpoint to the lookup table
        incidence_lookup['ALL CAUSE'] = incidence_lookup.apply(lambda x: self.get_incidence_lookup(incidence,
                                                                                              x['NAME'],
                                                                                              x['RACE'],
                                                                                              x['START_AGE'],
                                                                                              'ALL CAUSE'),axis=1)

        incidence_lookup['ISCHEMIC HEART DISEASE'] = incidence_lookup.apply(lambda x: self.get_incidence_lookup(incidence,
                                                                                              x['NAME'],
                                                                                              x['RACE'],
                                                                                              x['START_AGE'],
                                                                                              'ISCHEMIC HEART DISEASE'),axis=1)

        incidence_lookup['LUNG CANCER'] = incidence_lookup.apply(lambda x: self.get_incidence_lookup(incidence,
                                                                                              x['NAME'],
                                                                                              x['RACE'],
                                                                                              x['START_AGE'],
                                                                                              'LUNG CANCER'),axis=1)
        
        return incidence_lookup
    
    def combine_pop_inc(self):
        ''' Combines the population and incidence data into one dataset '''
        # Make copies of the population and incidence datasets
        population, incidence = self.population.copy(), self.incidence.copy()
        
        # Update both dataframes
        population = self.update_pop(population)
        incidence = self.update_inc(incidence)
        
        # Re-project incidence onto population
        incidence = incidence.to_crs(population.crs)
        
        ## Create intersect of just the geometries
        # Grab just the incidence geometry
        inc_geo = incidence[['NAME','geometry']].drop_duplicates()
        
        # Grab the population geometry and assign an ID for simplicity
        # ** THIS SHOULD BE UPDATED TO ISRM_ID
        pop_geo = population[['ISRM_ID','geometry']].drop_duplicates()

        # Perform the intersect
        pop_inc = gpd.overlay(pop_geo, inc_geo, how='intersection') 
        
        # # Merge in the population data on the ID field
        pop_inc = pd.merge(pop_inc, population[['ISRM_ID','START_AGE', 'END_AGE',
                                                'RACE', 'POPULATION']],
                            on='ISRM_ID')
        
        # Trim data to only include population age groups over 30 for faster processing
        pop_inc = pop_inc[pop_inc['START_AGE']>= 30]
        
        # Set up lookup keys
        pop_inc['KEY'] = pop_inc['START_AGE'].astype(str) + '_' + pop_inc['END_AGE'].astype(str) + '_' + pop_inc['NAME'].astype(str) + '_' + pop_inc['RACE'].astype(str)
        
        # Create a smaller incidence lookup tables for merging
        keys = pop_inc[['KEY']].drop_duplicates()
        lookups = self.incidence_by_age(incidence, population)
        lookup_dict = self.make_incidence_lookup(incidence, keys, lookups)
        
        # Map values from the lookup_dict
        pop_inc['INCIDENCE'] = pop_inc['KEY'].map(lookup_dict)
        pop_inc['ALL CAUSE INC'] = pop_inc['INCIDENCE'].str[0]
        pop_inc['ISCHEMIC HEART DISEASE INC'] = pop_inc['INCIDENCE'].str[1]
        pop_inc['LUNG CANCER INC'] = pop_inc['INCIDENCE'].str[2]
        
        # Clean up
        pop_inc = pop_inc[['ISRM_ID', 'NAME', 'RACE', 'POPULATION','ALL CAUSE INC', 
                           'ISCHEMIC HEART DISEASE INC','LUNG CANCER INC','geometry']]
        
        return pop_inc, lookup_dict
