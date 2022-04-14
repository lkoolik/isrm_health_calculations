#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main Run File

@author: libbykoolik
"""
#%% Import useful libraries
from pathlib import Path
import sys
sys.path.append('/Users/libbykoolik/Documents/Research/OEHHA Project/scripts/isrm_health_calculations/supporting')
from isrm import isrm
from emissions import emissions

#%% Import the argparse library to parse command line arguments
import argparse

# Initialize the parser object
parser = argparse.ArgumentParser(description="Runs the ISRM-based tool for estimating PM2.5 concentrations and associated health impacts.")

# Add necessary arguments
parser.add_argument("-c", "--checkinputs", help="use this run option if you want to check files but not run",
                    action='store_true')
parser.add_argument("-e", "--emissions", help="emissions input file path", type=str)
parser.add_argument("-u", "--units", help="emissions units")
parser.add_argument("-sn", "--scenarioname", help="name of scenario being modeled", type=str)

# Parse all arguments
args = parser.parse_args()

# Separate arguments into useful variables
check = args.checkinputs
if args.emissions: emissions_path = args.emissions
if args.units: units = args.units
if args.scenarioname: name = args.scenarioname

#%% Run Program
if __name__ == "__main__":

    # If check module is selected, run a file check and then exit without 
    # running calculations
    if check:
        emis = emissions(emissions_path, load_file=False, verbose=True)
        