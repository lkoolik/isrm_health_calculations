#!/isrm-env/bin/python python3
# -*- coding: utf-8 -*-
"""
Main Run File

@author: libbykoolik
Last updated: 2023-06-07
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
import concurrent.futures
import multiprocessing
import platform

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
  
#%% A few things need to go outside the __main__
#% Use argparse to parse command line arguments
# Initialize the parser object
parser = argparse.ArgumentParser(description="Runs the ISRM-based tool for estimating PM2.5 concentrations and associated health impacts.")

# Add necessary arguments
parser.add_argument("-i", "--inputs", help="control file path", type=str)
parser.add_argument("--debug", help="enable for debugging mode", action="store_true")

# Parse all arguments
args = parser.parse_args()
debug_mode = args.debug   

#%% Run Program
if __name__ == "__main__":        

    # Start the timer
    start_time = time.time()
    
    #% Create the log file and update logging configuration
    tmp_logger = setup_logging(debug_mode) 
    
    # Report the current version
    report_version()

    # Read control file and create useful variables
    cf = control_file(args.inputs)
    if not cf.ready:
        sys.exit()
    else: # load all the control file info
        batch = cf.batch_name
        name = cf.run_name
        emissions_path = cf.emissions_path
        units = cf.emissions_units
        isrm_path = cf.isrm_path
        population_path = cf.population_path
        run_health = cf.run_health
        race_stratified = cf.race_stratified
        check = cf.check
        verbose = cf.verbose
        region_of_interest = cf.region_of_interest
        region_category = cf.region_category
        output_resolution = cf.output_resolution
        output_exposure = cf.output_exposure
        detailed_conc_flag = cf.detailed_conc

    # Create the output directory
    output_dir, f_out = create_output_dir(batch, name)

    # Move the log file into the output directory
    new_logger = os.path.join(output_dir, 'log_'+f_out+'.txt')
    os.rename(tmp_logger, new_logger)
    
    # Save a copy of the control file into the output directory
    shutil.copy(args.inputs, output_dir)

    # Define data variable file paths
    ca_shp_path = './data/ca_border.feather'
    output_geometry_fps = {'AB': './data/air_basins.feather',
                           'AD': './data/air_districts.feather',
                           'C': './data/counties.feather'}
    incidence_fp = './data/benmap_incidence.feather'

    # Define output region based on region_of_interest and region_category
    output_region = get_output_region(region_of_interest, region_category, output_geometry_fps, ca_shp_path)

    # If check module is selected, run a file check and then exit without 
    # running calculations
    if check:
        try:
            # Default to verbose since this mode is just for checking files
            emis = emissions(emissions_path, units=units, name=name, load_file=False, verbose=True)
            isrmgrid = isrm(isrm_path, output_region, region_of_interest, load_file=False, verbose=True)
            pop = population(population_path, load_file=False, verbose=True)
            logging.info("\n<< Emissions, ISRM, and population files exist and are able to be imported. >>\n")
            
        except:
            logging.info("\n<< Correct error messages above before running the program. >>\n")
        quit()
    else: # for now, run everything
    
        # Setting the Multiprocessing startup mode should only be done once
        # in the lifetime of an application.
        multiprocessing.set_start_method('forkserver')
        
        # Create an output directory for shapefiles
        shape_out = create_shape_out(output_dir)
        
        ### CONCENTRATION MODULE
        logging.info('\n')
        logging.info('╓────────────────────────────────╖')
        logging.info('║ Beginning Concentration Module ║')
        logging.info('╙────────────────────────────────╜')
        logging.info('\n')

        ## Create emissions and ISRM objects
        # By using a ThreadPoolExecutor to spin up multiple threads, we can
        # read several files at the same time, and even let the population file
        # keep loading while we compute concentrations based on the emissions
        # and ISRM files. Multiple threads are a good option when the bottleneck
        # is disk I/O, not CPU usage.
        logging.info('<< Loading emissions, ISRM, and population files in parallel. Log messages may appear out of order. >>')        
        verboseprint(verbose, '- Processing for the emissions in verbose mode will be preceeded by [EMISSIONS].')
        verboseprint(verbose, '- Processing for the ISRM grid in verbose mode will be preceeded by [ISRM].')
        verboseprint(verbose, '- Processing for the population data in verbose mode will be preceeded by [POPULATION].')
        logging.info('\n')
        verboseprint(verbose, 'Details about file import:')
        file_reader_pool = concurrent.futures.ThreadPoolExecutor()
        emis_future = file_reader_pool.submit(emissions, emissions_path, units=units, name=name, load_file=True, verbose=verbose)
        isrm_future = file_reader_pool.submit(isrm, isrm_path, output_region, region_of_interest, load_file=True, verbose=verbose)
        pop_future = file_reader_pool.submit(population, population_path, load_file=True, verbose=verbose)
  
        # To run multiple computations at once, we need to create multiple
        # processes instead of threads. Processes take longer to create, but
        # they are necessary when we are paralellizing computation rather than
        # disk I/O.
        executor = concurrent.futures.ProcessPoolExecutor(max_workers=5)
        executor_jobs = []
        
        # We need to wait for the population and ISRM files to load before we
        # start creating the exp_pop_alloc or HIA inputs.
        pop = pop_future.result()
        isrmgrid = isrm_future.result()
        
        # Start estimating the concentrations by creating the exp_pop_alloc object.
        logging.info('\n<< Re-allocating population data to the ISRM grid >>')        
        verboseprint(verbose, '- This step will take some time, so the details about it may pop up in other sections.')
        verboseprint(verbose, '- Notes about this step will be preceded by the tag [POPULATION].')
        logging.info('\n')
        exp_pop_alloc_future = executor.submit(
            pop.allocate_population, pop.pop_exp, isrmgrid.geodata, 'ISRM_ID', False)
        executor_jobs.append(exp_pop_alloc_future)
        
        # Creating HIA inputs takes a long time. Start the process as early as 
        # possible, even if the rest of the health module won't start yet.
        if run_health:
            # Starts computing the HIA inputs now, in a parallel process.
            logging.info('\n<< Beginning to import health calculation inputs in parallel. Log messages may appear out of order. >>')
            verboseprint(verbose, '- Health calculation input details in verbose mode will be preceded by the tag [HEALTH].')
            hia_inputs_future = executor.submit(
                create_hia_inputs, pop, load_file=True, verbose=verbose, geodata=isrmgrid.geodata, incidence_fp=incidence_fp)
            executor_jobs.append(hia_inputs_future)
        
        ## Estimate concentrations
        logging.info('\n<< Estimating concentrations. >>')        
        verboseprint(verbose, '- Notes about this step will be preceded by the tag [CONCENTRATION].')
        logging.info('\n')
        emis = emis_future.result() # This almost always finishes earlier than the other files
        conc = concentration(emis, isrmgrid, detailed_conc_flag, run_calcs=True, verbose=verbose)


        ## Create plots and export results
        logging.info("\n<< Generating Concentration Outputs >>")
        verboseprint(verbose, '- Notes about this step will be preceded by the tag [CONCENTRATION].')
        logging.info('\n')
        conc.visualize_concentrations('TOTAL_CONC_UG/M3', output_region, output_dir, f_out, ca_shp_path, export=True)
        conc.export_concentrations(shape_out, f_out)
        logging.info("- [CONCENTRATION] Concentration files output into: {}.".format(output_dir))
        
        ## Best practice is to close the pool, but this causes an exception on Mac
        # Don't shut it down if it is still importing the health calculation inputs
        if platform.system() != 'Darwin' and run_health == False:
            conc_viz_future.result() # waits until this is done
            executor.shutdown(wait=False) # Closes this pool for any more tasks, but returns immediately
        
        ## Perform concentration-related EJ analyses
        exp_pop_alloc = pop.allocate_population(pop.pop_exp, isrmgrid.geodata, 'ISRM_ID', False)
        verboseprint(verbose, '- [POPULATION] Population data is properly allocated to the ISRM grid and ready for EJ calculations.')
        
        # Create the exposure dataframe and run EJ functions
        logging.info('\n<< Beginning Exposure EJ Calculations >>')
        verboseprint(verbose, '- Notes about this step will be preceded by the tag [EJ].')
        logging.info('\n')
        exposure_gdf, exposure_pctl, exposure_disparity = run_exposure_calcs(conc, exp_pop_alloc, verbose)    
        if output_exposure:
            export_exposure(exposure_gdf, shape_out, f_out)
        
        # Export results
        plot_percentile_exposure(output_dir, f_out, exposure_pctl, verbose)
        
        ### HEALTH MODULE
        if run_health:
            logging.info('\n╓────────────────────────────────╖')
            logging.info('║ Beginning Health Impact Module ║')
            logging.info('╙────────────────────────────────╜')
            logging.info('\n')
            
            # Wait for process creating health input object to finish
            verboseprint(verbose, '- [HEALTH] Waiting for health calculation inputs to finish')
            hia_inputs = hia_inputs_future.result()
            if platform.system() != 'Darwin':
                executor.shutdown(wait=False)
            verboseprint(verbose, '- [HEALTH] Health calculation inputs ready to proceed.')
            
            # Use the same process pool as above to estimate the three health endpoints. 
            with concurrent.futures.ProcessPoolExecutor(max_workers=5) as health_executor:
                logging.info('\n')
                logging.info('<< Beginning Health Calculations >>')
                verboseprint(verbose, '- The tool will estimate excess mortality from three endpoints in parallel. This results in a known bug that suppresses update statements and will be fixed in a future update. Note that this step may take a few minutes.')
                # verboseprint(verbose, '- The tool will estimate excess mortality from three endpoints in parallel. Log statements may appear out of order.')
                # verboseprint(verbose, '- Notes about All-Cause Mortality will be preceded by the tag [ACM].')
                # verboseprint(verbose, '- Notes about Ischemic Heart Disease Mortality will be preceded by the tag [IHD].')
                # verboseprint(verbose, '- Notes about Lung Cancer Mortality will be preceded by the tag [LCM].')
                logging.info('\n')

                # Estimate excess mortality
                logging.info('<< Estimating Excess Mortality for Three Endpoints >>')
                trimmed_conc = conc.detailed_conc_clean[['ISRM_ID','TOTAL_CONC_UG/M3','geometry']]
                pop = hia_inputs.population.groupby('ISRM_ID')[['ASIAN','BLACK','HISLA','INDIG','WHITE','TOTAL']].sum().reset_index()
                
                # Submit each endpoint as its own process to the health_executor
                allcause_future = health_executor.submit(calculate_excess_mortality, trimmed_conc,
                                                         hia_inputs.pop_inc, pop, 'ALL CAUSE', krewski, verbose)
                ihd_future = health_executor.submit(calculate_excess_mortality, trimmed_conc,
                                                         hia_inputs.pop_inc, pop, 'ISCHEMIC HEART DISEASE', krewski, verbose)
                lungcancer_future = health_executor.submit(calculate_excess_mortality, trimmed_conc,
                                                         hia_inputs.pop_inc, pop, 'LUNG CANCER', krewski, verbose)
                                
                # Collect all three results
                allcause = allcause_future.result()
                ihd = ihd_future.result()
                lungcancer = lungcancer_future.result()
                
                # Begin exporting the results in parallel
                logging.info('<< Exporting Health Impact Outputs >>')
                
                ## Set up futures
                allcause_ve_future = health_executor.submit(visualize_and_export_hia, allcause, ca_shp_path, 'TOTAL', 'ALL CAUSE', output_dir, f_out, shape_out, verbose=verbose)
                ihd_ve_future = health_executor.submit(visualize_and_export_hia, ihd, ca_shp_path, 'TOTAL', 'ISCHEMIC HEART DISEASE', output_dir, f_out, shape_out, verbose=verbose)
                lungcancer_ve_future = health_executor.submit(visualize_and_export_hia, lungcancer, ca_shp_path, 'TOTAL', 'LUNG CANCER', output_dir, f_out, shape_out, verbose=verbose)
                
                logging.info('- [HEALTH] Waiting for visualizations and exports to complete...')
                
                ## We don't actually need anything stored, we just need the program to wait until
                ## all three are done before exiting
                tmp = allcause_ve_future.result()
                tmp = ihd_ve_future.result()
                tmp = lungcancer_ve_future.result()
                
                ## Return that everything is done
                logging.info('- [HEALTH] All outputs have been exported!')
        
        # Final log statements to close out
        logging.info('\n')
        logging.info('╓────────────────────────────────╖')
        logging.info('║ Success! Run complete.         ║')
        logging.info('╙────────────────────────────────╜\n')
        logging.info('<< ISRM calculations tool has completed all calculations and exports. >>')
        run_time = round((time.time() - start_time)/60.0,1)
        logging.info('- Total run time: {} minutes'.format(run_time))
        logging.info('- All log statements have been saved into a text file: {}'.format(new_logger))
        logging.shutdown()
        
        quit()