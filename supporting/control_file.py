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
        self.keywords = ['RUN_NAME','EMISSIONS_FILENAME','EMISSIONS_UNITS',
                         'CHECK_INPUTS','VERBOSE']
        self.blanks_okay = [True, False, False, 
                            True, True]
        
        # Run basic checks on control file
        if self.valid_file:
            self.valid_structure, self.no_incorrect_blanks = self.check_control_file()
        else:
            self.valid_structure, self.no_incorrect_blanks = 'NA'
            
        # If checks are good, import values
        if self.valid_structure and self.no_incorrect_blanks and self.valid_file:
            self.run_name, self.emissions_path, self.emissions_units, self.check, self.verbose = self.get_all_inputs()
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
    
    def get_input_value(self, keyword):
        ''' Gets the input for the given keyword '''
        
        # Iterate through each line of the file to find the keyword
        for line in open(self.file_path):
            re_k = '- '+keyword+':' # Grabs exact formatting
            if re_k in line:
                line_val = line.split(':')[1].strip('\n').strip(' ')
            
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
        run_name = self.get_input_value('RUN_NAME')
        emissions_path = self.get_input_value('EMISSIONS_FILENAME')
        emissions_units = self.get_input_value('EMISSIONS_UNITS')
        
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
        
        return run_name, emissions_path, emissions_units, check, verbose 
    
    def check_inputs(self):
        ''' Once the inputs are imported, check them '''
        ## (1) Check the run name
        valid_run_name = type(self.run_name) == str
        print('* The run name provided is not valid (must be string).') if not valid_run_name else ''
        
        ## (2) Check the emissions path
        valid_emissions_path = self.check_path(file=self.emissions_path)
        print('* The emissions path provided is not valid.') if not valid_emissions_path else ''
        
        ## (3) Check the emissions units
        mass_units = ['ug','g','lb','ton','mt','kg'] # mass units from emissions.py
        time_units = ['s','min','hr','day','yr'] # time units from emissions.py

        # Try to split with the division sign
        tmp_units = self.emissions_units.lower().split('/')

        valid_emissions_units = tmp_units[0] in mass_units and tmp_units[1] in time_units
        print('* The emissions units provided is not valid. Acceptable mass emissions units are '+', '.join(mass_units)\
              +' and acceptable time units are '+', '.join(time_units)) if not valid_emissions_units else ''
        
        ## (4) Check the check flag
        valid_check = type(self.check) == bool
        print('* The check flag provided is not valid. Use Y or N or leave blank.') if not valid_check else ''
        
        ## (5) Check the verbose flag
        valid_verbose = type(self.verbose) == bool
        print('* The verbose flag provided is not valid. Use Y or N or leave blank.') if not valid_verbose else ''
        
        ## Output only one time
        valid_inputs = valid_run_name and valid_emissions_path and \
            valid_emissions_units and valid_check and valid_verbose

        return valid_inputs
