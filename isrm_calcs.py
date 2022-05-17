#!/isrm-env/bin/python python3
# -*- coding: utf-8 -*-
"""
Main Run File

@author: libbykoolik
Last updated: 2022-05-04
"""
#%% Import useful libraries
from pathlib import Path
import sys
sys.path.insert(0,'./supporting')
from isrm import isrm
from emissions import emissions
from concentration import concentration
from population import population
sys.path.insert(0,'./scripts')
from environmental_justice_calcs import *
import argparse
import os
import datetime

#%% Use argparse to parse command line arguments
# Initialize the parser object
parser = argparse.ArgumentParser(description="Runs the ISRM-based tool for estimating PM2.5 concentrations and associated health impacts.")

# Add necessary arguments
parser.add_argument("-c", "--checkinputs", help="use this run option if you want to check files but not run",
                    action='store_true')
parser.add_argument("-e", "--emissions", help="emissions input file path", type=str)
parser.add_argument("-u", "--units", help="emissions units")
parser.add_argument("-sn", "--scenarioname", help="name of scenario being modeled", type=str)
parser.add_argument("-v", "--verbose", help="include this tag to run the code in verbose mode", action='store_true')
parser.add_argument("-p", "--population", help="population input file path", type=str)

# Parse all arguments
args = parser.parse_args()

# Separate arguments into useful variables
check = args.checkinputs
verbose = args.verbose
if args.emissions: emissions_path = args.emissions
if args.population: population_path = args.population
units = args.units if args.units else 'ug/s'
name = args.scenarioname if args.scenarioname else ''

# Define ISRM Variables
isrm_fp = './data/ca_isrm.ncf'
isrm_gfp = './data/InMAP_gridCells.shp'

def create_output_dir():
    ''' Creates the output directory for files generated '''
    # Grab current working directory and the 'outputs' sub folder
    parent = os.getcwd()
    sub = 'outputs'
    
    # Output subdirectory will be named with current datetime
    now = datetime.datetime.now()
    outdir = 'out_'+name+now.strftime("%Y%m%d_%H%M")
    
    # Make the directory if it does not already exist
    os.mkdir(os.path.join(parent, sub, outdir))
    output_dir = os.path.join(parent,sub,outdir)
    
    # Print a statement to tell user where to look for files
    print("\n<< Output files created will be saved in the following directory: "+output_dir+">>")
    
    return output_dir

#%% Run Program
if __name__ == "__main__":

    # If check module is selected, run a file check and then exit without 
    # running calculations
    if check:
        try:
            # Default to verbose since this mode is just for checking files
            emis = emissions(emissions_path, units=units, name=name, load_file=False, verbose=True)
            isrmgrid = isrm(isrm_fp, isrm_gfp, load_file=False, verbose=True)
            pop = population(population_path, load_file=False, verbose=True)
            print("\n<< Emissions, ISRM, and population files exist and are able to be imported. >>\n")
        except:
            print("\n<< Correct error messages above before running the program. >>\n")
        quit()
    else: # for now, run concentration calculations since no health built
        output_dir = create_output_dir()
        emis = emissions(emissions_path, units=units, name=name, load_file=True, verbose=verbose)
        isrmgrid = isrm(isrm_fp, isrm_gfp, load_file=True, verbose=verbose)
        conc = concentration(emis, isrmgrid, run_calcs=True, verbose=verbose)
        
        print("\n<< Concentrations estimated >>\n")
        conc.visualize_concentrations('TOTAL_CONC_UG/M3',output_dir, export=True)
        conc.export_concentrations(output_dir, detailed=False)
        print("* Concentration files output into: {}.".format(output_dir))
        
        pop = population(population_path, load_file=True, verbose=verbose)
        pop_alloc = pop.allocate_population(isrmgrid.geodata, 'ISRM_ID')
        
        exposure_gdf = create_exposure_df(conc, pop_alloc)
        exposure_disparity = get_overall_disparity(exposure_gdf)
        exposure_pctl = estimate_exposure_percentile(exposure_gdf)
        plot_percentile_exposure(output_dir, exposure_pctl)
        
        quit()
        