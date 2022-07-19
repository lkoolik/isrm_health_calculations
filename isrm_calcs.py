#!/isrm-env/bin/python python3
# -*- coding: utf-8 -*-
"""
Main Run File

@author: libbykoolik
Last updated: 2022-07-19
"""
#%% Import useful libraries, supporting objects, and scripts
# Useful libraries for main script
from pathlib import Path
import sys
import argparse
import logging
import os
import time
import datetime
import shutil

# Import supporting objects
sys.path.insert(0,'./supporting')
from concentration import concentration
from concentration_layer import concentration_layer
from control_file import control_file
from emissions import emissions
from health_data import health_data
from isrm import isrm
from population import population

# Import supporting scripts
sys.path.insert(0,'./scripts')
from environmental_justice_calcs import *
from health_impact_calcs import *
from tool_utils import *


#%% Use argparse to parse command line arguments
start_time = time.time()

# Initialize the parser object
parser = argparse.ArgumentParser(description="Runs the ISRM-based tool for estimating PM2.5 concentrations and associated health impacts.")

# Add necessary arguments
parser.add_argument("-i", "--inputs", help="control file path", type=str)

# Parse all arguments
args = parser.parse_args()

# Create the log file and update logging configuration
tmp_logger = setup_logging()

# Read control file and create useful variables
cf = control_file(args.inputs)
if not cf.ready:
    sys.exit()
else:
    batch = cf.batch_name
    name = cf.run_name
    emissions_path = cf.emissions_path
    units = cf.emissions_units
    run_health = cf.run_health
    race_stratified = cf.race_stratified
    check = cf.check
    verbose = cf.verbose
    region_of_interest = cf.region_of_interest
    region_category = cf.region_category
    output_resolution = cf.output_resolution

# Create the output directory
output_dir, f_out = create_output_dir(batch, name)

# Move the log file into the output directory
new_logger = os.path.join(output_dir, 'log_'+f_out+'.txt')
os.rename(tmp_logger, new_logger)

# Save a copy of the control file into the output directory
shutil.copy(args.inputs, output_dir)

# Define data variable file paths
isrm_fps = ['./data/ISRM_NH3.npy','./data/ISRM_NOX.npy','./data/ISRM_PM25.npy',
            './data/ISRM_SOX.npy','./data/ISRM_VOC.npy']
isrm_gfp = './data/isrm_geo_test.feather'
population_path = './data/ca2000.feather'
ca_shp_path = './data/ca_border.feather'
output_geometry_fps = {'AB': './data/air_basins.feather',
                       'AD': './data/air_districts.feather',
                       'C': './data/counties.feather'}
hia_input_fps = {'POPULATION': './data/benmap_population_new.feather',
                  'INCIDENCE': './data/benmap_incidence.feather'}

# Define output region based on region_of_interest and region_category
output_region = get_output_region(region_of_interest, region_category, output_geometry_fps, ca_shp_path)

#%% Run Program
if __name__ == "__main__":    
    # If check module is selected, run a file check and then exit without 
    # running calculations
    if check:
        try:
            # Default to verbose since this mode is just for checking files
            emis = emissions(emissions_path, units=units, name=name, load_file=False, verbose=True)
            isrmgrid = isrm(isrm_fps, isrm_gfp, output_region, region_of_interest, load_file=False, verbose=True)
            pop = population(population_path, load_file=False, verbose=True)
            logging.info("\n<< Emissions, ISRM, and population files exist and are able to be imported. >>\n")
        except:
            logging.info("\n<< Correct error messages above before running the program. >>\n")
        quit()
    else: # for now, run everything
        # Create an output directory for shapefiles
        shape_out = create_shape_out(output_dir)
        
        ### CONCENTRATION MODULE
        logging.info('\n ╓────────────────────────────────╖')
        logging.info('║ Beginning Concentration Module ║')
        logging.info('╙────────────────────────────────╜\n')
        
        ## Create emissions and ISRM objects
        emis = emissions(emissions_path, units=units, name=name, load_file=True, verbose=verbose)
        isrmgrid = isrm(isrm_fps, isrm_gfp, output_region, region_of_interest, load_file=True, verbose=verbose)
        
        ## Estimate concentrations
        conc = concentration(emis, isrmgrid, run_calcs=True, verbose=verbose)
        
        ## Create plots and export results
        logging.info("<< Generating Concentration Outputs >>")
        conc.visualize_concentrations('TOTAL_CONC_UG/M3',output_region, output_dir, f_out, ca_shp_path, export=True)
        conc.export_concentrations(shape_out, f_out, detailed=False)
        logging.info("- Concentration files output into: {}.".format(output_dir))
        
        ## Perform concentration-related EJ analyses
        # Create a population object and intersect population with concentrations
        pop = population(population_path, load_file=True, verbose=verbose)
        pop_alloc = pop.allocate_population(isrmgrid.geodata, 'ISRM_ID')
        
        # Create the exposure dataframe and run EJ functions
        logging.info('\n << Beginning Exposure EJ Calculations >>')
        exposure_pctl, exposure_disparity = run_exposure_calcs(conc, pop_alloc, verbose)
        
        # Export results
        plot_percentile_exposure(output_dir, f_out, exposure_pctl, verbose)
        
        ### HEALTH MODULE
        if run_health:
            logging.info('\n ╓────────────────────────────────╖')
            logging.info('║ Beginning Health Impact Module ║')
            logging.info('╙────────────────────────────────╜\n')
            
            # Create health input object
            hia_inputs = health_data(hia_input_fps, verbose=verbose, race_stratified=False)
            
            # Estimate excess mortality
            logging.info('\n << Estimating Excess Mortality for Three Endpoints >>')
            allcause = calculate_excess_mortality(conc, hia_inputs, 'ALL CAUSE', krewski, verbose=verbose)
            ihd = calculate_excess_mortality(conc, hia_inputs, 'ISCHEMIC HEART DISEASE', krewski, verbose=verbose)
            lungcancer = calculate_excess_mortality(conc, hia_inputs, 'LUNG CANCER', krewski, verbose=verbose)            
            
            # Plot and export
            logging.info('\n<< Health Impact Outputs >>')
            visualize_and_export_hia(allcause, ca_shp_path, 'TOTAL', 'ALL CAUSE', output_dir, f_out, shape_out, verbose=verbose)
            visualize_and_export_hia(ihd, ca_shp_path, 'TOTAL', 'ISCHEMIC HEART DISEASE', output_dir, f_out, shape_out, verbose=verbose)
            visualize_and_export_hia(lungcancer, ca_shp_path, 'TOTAL', 'LUNG CANCER', output_dir, f_out, shape_out, verbose=verbose)
            pass
        
        # Final log statements
        logging.info('\n ╓────────────────────────────────╖')
        logging.info('║ Success! Run complete.         ║')
        logging.info('╙────────────────────────────────╜\n')
        logging.info('<< ISRM calculations tool has completed all calculations and exports. >>')
        run_time = round((time.time() - start_time)/60.0,0)
        logging.info('- Total run time: {} minutes'.format(run_time))
        logging.info('- All log statements have been saved into a text file: {}'.format(new_logger))
        
        quit()
        