#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tool Utils

@author: libbykoolik
last modified: 2023-06-07
"""

#%% Import useful libraries
from pathlib import Path
import sys
sys.path.insert(0,'./supporting')
sys.path.insert(0,'./scripts')
import argparse
import os
from os import path
import datetime
import geopandas as gpd
import logging


#%% ISRM Tool Utils
def check_setup():
    '''
    Checks that the isrm_health_calculations local clone is set up properly
    
    INPUTS: None
        
    OUTPUTS:
        - valid_setup: a Boolean indicating if the setup is correct or not
        
    '''
    
    ## First, get the current working directory
    cwd = os.getcwd()
    
    ## Set up an error tracker
    errors = 0
    
    ## Next, check that each code file exists where it should
    # Did the scripts and supporting folders get downloaded?
    script_folder_flag = path.exists(path.join(cwd, 'scripts'))
    support_folder_flag = path.exists(path.join(cwd, 'supporting'))
    
    # If scripts folder is there, check each file
    if script_folder_flag:
        
        # Hard-coded files
        scripts = ['environmental_justice_calcs.py', 'health_impact_calcs.py', 'tool_utils.py']
    
        # Check each one individually
        for script in scripts:
            tmp_path = path.join(cwd, 'scripts', script)
            if not path.exists(tmp_path) or not path.isfile(tmp_path):
                print('* Missing script file in the scripts directory: {}'.format(script))
                errors += 1
                
    else: # Tell user everything is missing
        print('* Missing script folder, which should contain files: ', ', '.join(scripts))
        errors += 1
    
    # If scripts folder is there, check each file
    if support_folder_flag:
        
        # Hard-coded files
        supports = ['concentration.py', 'concentration_layer.py', 'control_file.py', 'emissions.py',
                    'health_data.py', 'isrm.py', 'population.py']
        
        # Check each one individually
        for support in supports:
            tmp_path = path.join(cwd, 'supporting', support)
            if not path.exists(tmp_path) or not path.isfile(tmp_path):
                print('* Missing supporting file in the supporting directory: {}'.format(support))
                errors += 1
                
    else: # Tell user everything is missing
        print('* Missing supporting folder, which should contain files: ', ', '.join(supports))
        errors += 1
            
    ## Next, check for the necessary directories
    for sub in ['outputs', 'templates']:
        if not path.exists(path.join(cwd, sub)):
            print('* Missing {} sub-directory in isrm_health_calculations.'.format(sub))
            errors += 1
    
    ## Finally, check if all data are present.
    # First ensure that the data folder exists
    data_folder_flag = path.exists(path.join(cwd, 'data'))
    
    # If there is a data folder, continue with check
    if data_folder_flag:
        necessary_data = ['air_basins.feather', 'air_districts.feather', 'benmap_incidence.feather',
                          'ca_border.feather', 'ca2010.feather', 'counties.feather']
        
        for data in necessary_data:
            tmp_path = path.join(cwd, 'data', data)
            if not path.exists(tmp_path) or not path.isfile(tmp_path):
                print('* Missing data file in the data directory: {}'.format(data))
                errors += 1
                
    else: # Missing data folder
        print('* Missing data sub directory in isrm_health_calculations.')
        errors += 1
        
    ## Define the setup as valid only if there are no errors
    valid_setup = errors == 0
    
    ## Check if ISRM is there, but don't count as error if not there (just alert user)
    if not path.exists(path.join(cwd, 'data', 'CA_ISRM')):
        print('- No CA_ISRM file found in the data directory. Be sure to supply the currect filepath of your ISRM directory when running the tool.')
        
    else:
        isrm_files = ['isrm_geo.feather', 'ISRM_NH3.npy', 'ISRM_NOX.npy', 'ISRM_PM25.npy',
                      'ISRM_SOX.npy', 'ISRM_VOC.npy']
        for f in isrm_files:
            tmp_path = path.join(cwd, 'data', 'CA_ISRM', f)
            if not path.exists(tmp_path) or not path.isfile(tmp_path):
                print('- Found CA_ISRM but the folder is missing the following file: {}'.format(f))
    
    if valid_setup:
        print('Set-up is correctly done. Script files, supporting script files, and key data files are stored in the proper place.')
        
    return valid_setup

def setup_logging(debug_mode):
    '''
    Sets up the log file system for runs.
    
    INPUTS:
        - debug_mode: if true, will output logging statements in debugging mode
        
    OUTPUTS:
        - tmp_logger: a filepath string associated with a temporary log file that will 
          be moved as soon as the output directory is created

    '''
    # Create Temporary Logging File (will be renamed)
    tmp_logger = os.path.join(os.getcwd(),'tmp.txt')
    
    if debug_mode:
        level = logging.INFO
        format = '%(asctime)s %(filename)s:%(lineno)s %(message)s'
        datefmt = '%Y-%m-%d %H:%M:%S'
        handlers = [logging.FileHandler(tmp_logger), logging.StreamHandler()]
        logging.basicConfig(level=level, format=format, handlers=handlers, datefmt=datefmt)
    
        # Suppress all other library warnings and information
        for key in logging.Logger.manager.loggerDict:
            logging.getLogger(key).setLevel(logging.CRITICAL)
            
    else:
        level = logging.INFO
        format = '%(message)s'
        handlers = [logging.FileHandler(tmp_logger), logging.StreamHandler()]
        logging.basicConfig(level=level, format=format, handlers=handlers)
    
        # Suppress all other library warnings and information
        for key in logging.Logger.manager.loggerDict:
            logging.getLogger(key).setLevel(logging.CRITICAL)

    return tmp_logger

def verboseprint(verbose, text):
    '''
    Sets up the verbose printing mechanism. Adding here makes it global.
    
    INPUTS:
        - verbose: if true, will output a print statement to console and log file
        - text: the string that should be printed if running in verbose
        
    OUTPUTS: None
    
    '''
    if verbose:
        logging.info(text)
    else: lambda *a, **k:None
    return
    
def report_version():      
    '''
    Reports the current working version of the tool.
    
    INPUTS: None
        
    OUTPUTS: None
    
    '''

    logging.info('╔════════════════════════════════╗')
    logging.info('║ ISRM Health Calculations Tool  ║')
    logging.info('║ Version 0.8.1                  ║')
    logging.info('╚════════════════════════════════╝')
    logging.info('\n')
    return

def create_output_dir(batch, name):
    ''' 
    Creates the output directory for files generated.
    
    INPUTS:
        - batch: the batch name 
        - name: the run name
        
    OUTPUTS:
        - output_dir: a filepath string for the output directory
        - f_out: a string containing the filename pattern to be used in output files
    
    '''
    # Grab current working directory and the 'outputs' sub folder
    parent = os.getcwd()
    sub = 'outputs'
    
    # Output subdirectory will be named with batch and name
    name_list = ['out',batch,name]
    try: name_list.remove('')
    except ValueError: pass
    outdir = '_'.join(name_list)
    
    # If the directory already exists, add an integer to the end
    path_exists_flag = path.exists(os.path.join(parent,sub,outdir))
    # Use while loop in case there are multiple directories already made
    while path_exists_flag:
        if outdir == '_'.join(name_list):
            n = 0
            outdir_tmp = outdir + '_'        
        else:
            # Need to pull just the last two numbers 
            try:
                n = int(outdir.split('_')[-1]) # Grab the number at the end
                outdir_tmp = outdir[:-2] # Get a temporary output directory
                if outdir_tmp[-1] != '_': # Only one run exists in the directory, need to start at 01
                    n = 0
                    outdir_tmp = outdir_tmp + '_'
            except ValueError: #This means there was no number at the end, so start with 0
                outdir_tmp = outdir + '_'
                n = 0
        # Update to the next n
        next_n = str(n+1).zfill(2)
        outdir = outdir_tmp+next_n
        
        # Check if this path exists
        path_exists_flag = path.exists(os.path.join(parent,sub,outdir))
        
    # Add a new variable f_out that adds more information to output file names
    f_out = outdir[4:]# Cuts off the 'out_' from the start
    
    # Make the directory if it does not already exist
    os.mkdir(os.path.join(parent, sub, outdir))
    output_dir = os.path.join(parent,sub,outdir)
    
    # Print a statement to tell user where to look for files
    logging.info("\n << Output files created will be saved in the following directory: "+output_dir+">>")
    
    return output_dir, f_out

def create_shape_out(output_dir):
    ''' 
    Create additional subdirectory for shapefiles.
    
    INPUTS:
        - output_dir: a filepath string for the output directory
        
    OUTPUTS:
        - shape_out: a filepath string for the shapefile output directory
    
    '''
    # Set up parent and sub
    parent = output_dir
    sub = 'shapes'
    
    # Make the directory and store directory path
    os.mkdir(os.path.join(parent, sub))
    shape_out = os.path.join(parent,sub)
    
    return shape_out

def get_output_region(region_of_interest, region_category, output_geometry_fps, ca_fps):
    ''' 
    Outputs a geodataframe of the output region.
    
    INPUTS:
        - region_of_interest:  the name of the region to be contained in the `output_region`
        - region_category: a string containing the region category for the output region, 
          must be one of 'AB','AD', or 'C' for Air Basins, Air Districts, and Counties
        - output_geometry_fps: a dictionary containing a mapping between `region_category` 
          and the filepaths
        - ca_fps: a filepath string containing the link to the California border shapefile
        
    OUTPUTS:
        - output_region: a geodataframe containing only the region of interest
        
    '''
    if region_of_interest != 'CA': # Check if output region even wanted
        sys.path.append(os.path.realpath('..'))
        
        # Get the filepath
        geo_file_fp = output_geometry_fps[region_category]
    
        # Read in the file as a geodataframe
        geo_file = gpd.read_feather(geo_file_fp)
        
        # Clip to the region of interest
        output_region = geo_file[geo_file['NAME']==region_of_interest]
        
    else: # We want all of the domain, so return all of CA
        output_region = gpd.read_feather(ca_fps)
        
    return output_region