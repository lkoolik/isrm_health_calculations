#!/isrm-env/bin/python python3
# -*- coding: utf-8 -*-
"""
Main Run File

@author: libbykoolik
Last updated: 2022-06-09
"""
#%% Import useful libraries
from pathlib import Path
import sys
sys.path.insert(0,'./supporting')
from isrm import isrm
from emissions import emissions
from concentration_layer import concentration_layer
from concentration import concentration
from population import population
from control_file import control_file
sys.path.insert(0,'./scripts')
from environmental_justice_calcs import *
from tool_utils import *
import argparse
import os
import datetime

#%% Use argparse to parse command line arguments
# Initialize the parser object
parser = argparse.ArgumentParser(description="Runs the ISRM-based tool for estimating PM2.5 concentrations and associated health impacts.")

# Add necessary arguments
parser.add_argument("-i", "--inputs", help="control file path", type=str)

# Parse all arguments
args = parser.parse_args()

# Read control file and create useful variables
cf = control_file(args.inputs)
if not cf.ready:
    sys.exit()
else:
    batch = cf.batch_name
    name = cf.run_name
    emissions_path = cf.emissions_path
    units = cf.emissions_units
    check = cf.check
    verbose = cf.verbose

# Define ISRM Variables and population variables
isrm_fps = ['./data/ISRM_NH3.npy','./data/ISRM_NOX.npy','./data/ISRM_PM25.npy',
            './data/ISRM_SOX.npy','./data/ISRM_VOC.npy']
isrm_gfp = './data/isrm_geo.feather'
population_path = './data/ca2000.feather'
ca_shp_path = './data/ca_state.feather'

#%% Run Program
if __name__ == "__main__":

    # If check module is selected, run a file check and then exit without 
    # running calculations
    if check:
        try:
            # Default to verbose since this mode is just for checking files
            emis = emissions(emissions_path, units=units, name=name, load_file=False, verbose=True)
            isrmgrid = isrm(isrm_fps, isrm_gfp, load_file=False, verbose=True)
            pop = population(population_path, load_file=False, verbose=True)
            print("\n<< Emissions, ISRM, and population files exist and are able to be imported. >>\n")
        except:
            print("\n<< Correct error messages above before running the program. >>\n")
        quit()
    else: # for now, run concentration calculations since no health built
        output_dir, f_out = create_output_dir(batch, name)
        emis = emissions(emissions_path, units=units, name=name, load_file=True, verbose=verbose)
        isrmgrid = isrm(isrm_fps, isrm_gfp, load_file=True, verbose=verbose)
        conc = concentration(emis, isrmgrid, run_calcs=True, verbose=verbose)
        
        print("\n<< Concentrations estimated >>\n")
        conc.visualize_concentrations('TOTAL_CONC_UG/M3',output_dir, f_out, ca_shp_path, export=True)
        conc.export_concentrations(output_dir, f_out, detailed=False)
        print("* Concentration files output into: {}.".format(output_dir))
        
        pop = population(population_path, load_file=True, verbose=verbose)
        pop_alloc = pop.allocate_population(isrmgrid.geodata, 'ISRM_ID')
        
        exposure_gdf = create_exposure_df(conc, pop_alloc)
        exposure_disparity = get_overall_disparity(exposure_gdf)
        exposure_pctl = estimate_exposure_percentile(exposure_gdf)
        plot_percentile_exposure(output_dir, f_out, exposure_pctl)
        
        quit()
        