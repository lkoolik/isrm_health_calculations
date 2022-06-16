#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Control File Reading Object

@author: libbykoolik
last modified: 2022-06-07
"""

# Import Libraries
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import os
from os import path
import sys
import re
import difflib
            
class control_file:
    '''
    Defines a new object for reading the control file and storing information.
    
    INPUTS:
        - file_path: File path to the control file 
        
    '''
    def __init__(self, file_path):
        ''' Initializes the control file object'''        
        # Initialize path and check that it is valid
        self.file_path = file_path
        self.valid_file = self.check_path()
        
        # Hardcode the current keywords for simplicity
        self.keywords = ['BATCH_NAME', 'RUN_NAME','EMISSIONS_FILENAME',
                         'EMISSIONS_UNITS','CHECK_INPUTS','VERBOSE',
                         'REGION_OF_INTEREST','REGION_CATEGORY','OUTPUT_RESOLUTION']
        self.blanks_okay = [True, True, False, 
                            False, True, True,
                            True, True, True]
        
        # Run basic checks on control file
        if self.valid_file:
            self.valid_structure, self.no_incorrect_blanks = self.check_control_file()
        else:
            self.valid_structure, self.no_incorrect_blanks = 'NA'
            
        # If checks are good, import values
        if self.valid_structure and self.no_incorrect_blanks and self.valid_file:
            self.batch_name, self.run_name, self.emissions_path, self.emissions_units, self.check, self.verbose, self.region_of_interest, self.region_category, self.output_resolution = self.get_all_inputs()
            self.valid_inputs = self.check_inputs()
            if self.valid_inputs:
                print('\n<< Control file was successfully imported and inputs are correct >>')
                self.ready = True
            else:
                print('\n<< Control file was successfully imported but inputs are not correct >>')
                self.ready = False
        else:
            if not self.valid_structure:
                print('\n* Control file did not have the correct structure.')
                print('* Please confirm that the control file has the following keywords exactly once each:')
                print('   * '+'\n   * '.join(self.keywords))
            if not self.no_incorrect_blanks:
                print('* Some keywords were left blank incorrectly in the control file.')
                print('* Only the following keywords can be left blank:')
                print('   * '+'\n   * '.join(pd.Series(self.keywords)[self.blanks_okay].tolist()))
            if not self.valid_file:
                print('* The control file path is not valid.')
            
            print('\n<< Control File was not correct. Please revise errors and try again. >>\n')
            self.ready = False
            
    def check_path(self, file=''):
        ''' Checks if file exists at the path specified '''
        # Use the os library to check the path and the file
        if file == '':
            path_exists = path.exists(self.file_path)
            file_exists = path.isfile(self.file_path)
        else:
            path_exists = path.exists(file)
            file_exists = path.isfile(file)
        return path_exists and file_exists
    
    def get_input_value(self, keyword, upper=False):
        ''' Gets the input for the given keyword '''
        
        # Iterate through each line of the file to find the keyword
        for line in open(self.file_path):
            re_k = '- '+keyword+':' # Grabs exact formatting
            if re_k in line:
                line_val = line.split(':')[1].strip('\n').strip(' ')
            
        if upper: # Should be uppercased
            line_val = line_val.upper()
            
        return line_val
    
    
    def check_control_file(self):
        ''' Runs a number of checks to make sure that control file is valid  '''
        ## (TEST 1) Check for correct keywords
        # Create a dictionary that logs whether or not each keyword exists and 
        # how many times
        check_dict = {}
        
        # Loop through the text file and capture keyword usage
        for line in open(self.file_path):
            for k in self.keywords:
                re_k = '- '+k+':' # Grabs exact formatting
                if re_k in line and k not in check_dict.keys(): # Found keyword, doesn't exist in dictionary
                    check_dict[k] = 1
                elif k in line and k in check_dict.keys(): # Found keyword, already exists in dictionary
                    check_dict[k] += 1
                else:
                    pass # not found, keep looping
                    
        # confirm that keywords all present
        all_keywords = set(self.keywords) == set(check_dict.keys())
        
        # confirm that all values are 1
        correct_count = sum(check_dict.values()) == len(self.keywords)
        
        # Combine for test #1 output
        valid_structure = all_keywords & correct_count
        
        ## (TEST 2) Check for blank inputs
        incorrect_blanks = 0 # holder for incorrect blanks
        for line in open(self.file_path):
            for k in zip(self.keywords, self.blanks_okay):
                if k[1]: # If blanks are okay, ignore
                    pass
                else:
                    line_val = self.get_input_value(k[0])
                    if line_val == '': # Blanks will report as ''
                        incorrect_blanks += 1 # Add to holder
        no_incorrect_blanks = incorrect_blanks == 0

        return valid_structure, no_incorrect_blanks
    
    
    def get_all_inputs(self):
        ''' Once it passes the basic control file checks, import the values '''
        mapper = {'Y':True, 'N':False}
        
        # Get each input individually
        batch_name = self.get_input_value('BATCH_NAME')
        run_name = self.get_input_value('RUN_NAME')
        emissions_path = self.get_input_value('EMISSIONS_FILENAME')
        emissions_units = self.get_input_value('EMISSIONS_UNITS')
        region_of_interest = self.get_input_value('REGION_OF_INTEREST', upper=True)
        region_category = self.get_input_value('REGION_CATEGORY', upper=True)
        output_resolution = self.get_input_value('OUTPUT_RESOLUTION', upper=True)
        
        # For CHECK_INPUTS and VERBOSE, assume something if no value is given
        check = self.get_input_value('CHECK_INPUTS')
        if check == '':
            print('* No value provided for the CHECK_INPUTS field. Assuming a full run.')
            check = False
        else:
            check = mapper[check] # convert Y/N to True/False
        verbose = self.get_input_value('VERBOSE')
        if verbose == '':
            print('* No value provided for the VERBOSE field. Assuming a verbose run.')
            verbose = True
        else:
            verbose = mapper[verbose] # convert Y/N to True/False
            
        # For OUTPUT OPTIONS, assume something if no value is given
        if region_of_interest == '':
            print('* No value provided for the REGION_OF_INTEREST field. Assuming full California run.')
            region_of_interest = 'CA'
        if region_category == '':
            print('* No value provided for the REGION_CATEGORY field. Assuming full California run.')
            region_category = 'STATE'
        if output_resolution == '':
            print('* No value provided for the OUTPUT_RESOLUTION field. Assuming ISRM grid cells.')
            output_resolution = 'ISRM'
        
        return batch_name, run_name, emissions_path, emissions_units, check, verbose, region_of_interest, region_category, output_resolution
    
    def get_region_dict(self):
        ''' Hard-coded dictionary of acceptable values for regions '''
        # Define lists and a dictionary of acceptable inputs
        air_basins = ['GREAT BASIN VALLEYS','LAKE COUNTY','LAKE TAHOE','MOJAVE DESERT',
                      'MOUNTAIN COUNTIES','NORTH CENTRAL COAST','NORTH COAST',
                      'NORTHEAST PLATEAU','SACRAMENTO VALLEY','SALTON SEA','SAN DIEGO COUNTY',
                      'SAN FRANCISCO BAY','SAN JOAQUIN VALLEY','SOUTH CENTRAL COAST','SOUTH COAST']
        
        air_districts = ['AMADOR','ANTELOPE VALLEY','BAY AREA','BUTTE','CALAVERAS','COLUSA',
                         'EL DORADO','FEATHER RIVER','GLENN','GREAT BASIN UNIFIED','IMPERIAL',
                         'KERN','LAKE','LASSEN','MARIPOSA','MENDOCINO','MODOC','MOJAVE DESERT',
                         'MONTEREY BAY UNIFIED','NORTH COAST UNIFIED','NORTHERN SIERRA',
                         'NORTHERN SONOMA','PLACER','SACRAMENTO METRO','SAN DIEGO',
                         'SAN JOAQUIN VALLEY UNIFIED','SAN LUIS OBISPO','SANTA BARBARA','SHASTA',
                         'SISKIYOU','SOUTH COAST','TEHAMA','TUOLUMNE','VENTURA','YOLO-SOLANO']
        
        counties = ['ALAMEDA','ALPINE','AMADOR','BUTTE','CALAVERAS','COLUSA','CONTRA COSTA',
                    'DEL NORTE','EL DORADO','FRESNO','GLENN','HUMBOLDT','IMPERIAL','INYO',
                    'KERN','KINGS','LAKE','LASSEN','LOS ANGELES','MADERA','MARIN','MARIPOSA',
                    'MENDOCINO','MERCED','MODOC','MONO','MONTEREY','NAPA','NEVADA','ORANGE',
                    'PLACER','PLUMAS','RIVERSIDE','SACRAMENTO','SAN BENITO','SAN BERNARDINO',
                    'SAN DIEGO','SAN FRANCISCO','SAN JOAQUIN','SAN LUIS OBISPO','SAN MATEO',
                    'SANTA BARBARA','SANTA CLARA','SANTA CRUZ','SHASTA','SIERRA','SISKIYOU',
                    'SOLANO','SONOMA','STANISLAUS','SUTTER','TEHAMA','TRINITY','TULARE',
                    'TUOLUMNE','VENTURA','YOLO','YUBA']
        
        dacs = counties.copy() # DACs are listed by county, and the names are the same -- ignore for now
        
        region_dict = {'AB':air_basins, 'AD':air_districts, 'C':counties}#, 'DAC':dacs}
        
        return region_dict
    
    def region_check_helper(self):
        ''' Simple helper function for checking the region of interest and region_category inputs '''
        # Import the region dict
        region_dict = self.get_region_dict()
        
        ## First check that the REGION_CATEGORY is okay
        region_category = self.region_category
        
        # Define some possible but acceptable wrong names
        possible_wrong_names_dict = {'AB':['AIR BASIN', 'AIRBASIN', 'AIR_BASIN', 'BASIN'],
                                     'AD':['AIR DISTRICT', 'AIRDISTRICT', 'AIR_DISTRICT', 'DISTRICT'],
                                     'C':['COUNTY', 'COUNTIES']}#,
                                     #'DAC':['DISADVANTAGED COMMUNITY', 'EJ COMMUNITY']} # to add DACs in future update
        
        # Check if region_category is correct
        if region_category not in region_dict.keys() and region_category != 'STATE':
            # Set possible replacement to NA as a holder value and create flag variable
            possible_replacement = 'NA'
            region_category_flag = False
            
            # Loop through the dictionary
            for region in possible_wrong_names_dict.keys():
                if region_category in possible_wrong_names_dict[region]: # If it's in the list, update possible_replacement
                    possible_replacement = region
                    region_category_flag = True
                    print('* Incorrect value provided for the REGION_CATEGORY field, however a close replacement ({}) was found.'.format(possible_replacement))
            
            # Update region_category
            region_category = possible_replacement
            
        else:
            region_category_flag = True
                
        if region_category_flag and region_category != 'STATE':
            ## Second, check the region of interest
            region_of_interest = self.region_of_interest
            
            # Get a list of the options
            region_opts = region_dict[region_category]
            
            # Check if region_of_interest is in region_opts
            if region_of_interest not in region_opts and region_of_interest != 'CA':
                # Set possible replacement to NA as a holder value and create flag variable
                possible_replacement = 'NA'
                region_of_interest_flag = False
                
                # Use the difflib function get_close_matches to find the closest match
                closest_match = difflib.get_close_matches(region_of_interest, region_opts, n=1, cutoff=0.6)
                
                try: # If the list does not come back empty, should be able to grab the first entry 
                    closest_match = closest_match[0]
                    print('* Incorrect value provided for the REGION_OF_INTEREST field, however a close replacement ({}) was found.'.format(closest_match))
                    region_of_interest_flag = True
                    region_of_interest = closest_match
                    
                except: # Do nothing, will return region_of_interest = False
                    pass
                    
            else:
                region_of_interest_flag = True
                
            ## Update the self variables for region_category and region_of_interest
            self.region_category = region_category
            self.region_of_interest = region_of_interest
        
        elif region_category_flag and region_category == 'STATE':
            region_of_interest_flag = True
        
        else:
            region_of_interest_flag = False
        
        return region_category_flag, region_of_interest_flag
    
    def check_inputs(self):
        ''' Once the inputs are imported, check them '''
        ## (1) Check the batch name and run name are strings
        valid_batch_name = type(self.batch_name) == str
        print('* The batch name provided is not valid (must be string).') if not valid_batch_name else ''
        
        valid_run_name = type(self.run_name) == str
        print('* The run name provided is not valid (must be string).') if not valid_run_name else ''
        
        ## (2) If batch and run name are blank, replace with 'isrm_calcs'
        if self.batch_name == '' and self.run_name == '':
            self.batch_name = 'isrm_calcs'
        
        ## (3) Check the emissions path
        valid_emissions_path = self.check_path(file=self.emissions_path)
        print('* The emissions path provided is not valid.') if not valid_emissions_path else ''
        
        ## (4) Check the emissions units
        mass_units = ['ug','g','lb','ton','mt','kg'] # mass units from emissions.py
        time_units = ['s','min','hr','day','yr'] # time units from emissions.py

        # Try to split with the division sign
        tmp_units = self.emissions_units.lower().split('/')

        valid_emissions_units = tmp_units[0] in mass_units and tmp_units[1] in time_units
        print('* The emissions units provided is not valid. Acceptable mass emissions units are '+', '.join(mass_units)\
              +' and acceptable time units are '+', '.join(time_units)) if not valid_emissions_units else ''
        
        ## (5) Check the check flag
        valid_check = type(self.check) == bool
        print('* The check flag provided is not valid. Use Y or N or leave blank.') if not valid_check else ''
        
        ## (6) Check the verbose flag
        valid_verbose = type(self.verbose) == bool
        print('* The verbose flag provided is not valid. Use Y or N or leave blank.') if not valid_verbose else ''
        
        ## (7) Check the region and region_of_interest
        valid_region_category, valid_region_of_interest = self.region_check_helper()
        print('* The region category provided is not valid. Valid options include: AB, AD, and C. For fully state, leave blank.') if not valid_region_category else ''
        print('* The region of interest provided is not valid or cannot be validated. Consult user manual for correct options.') if not valid_region_of_interest else ''
        
        ## (8) Check the output_resolution variable
        valid_output_resolutions = ['ISRM', 'C', 'AB', 'AD']# to add DACs in future update, 'DACS']
        valid_output_resolution = self.output_resolution in valid_output_resolutions
        print('* The output resolution provided is not valid. Valid options include: '+', '.join(valid_output_resolutions)) if not valid_output_resolution else ''
        
        ## Output only one time
        valid_inputs = valid_batch_name and valid_run_name and valid_emissions_path and \
            valid_emissions_units and valid_check and valid_verbose and \
                valid_region_category and valid_region_of_interest and valid_output_resolution

        return valid_inputs
