#!/isrm-env/bin/python python3
# -*- coding: utf-8 -*-
"""
Main Run File

@author: libbykoolik
Last updated: 2023-03-15
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



def create_hia_inputs(pop, load_file: bool, verbose: bool, geodata:pd.DataFrame,
                      incidence_fp: str):
    """ Creates the hia_inputs object.
    
    Moving this into a separate function allows us to run this in parallel while
    other functions are running, speeding up the overall execution of the
    application.
    """
    hia_pop_alloc = pop.allocate_population(pop.pop_all, geodata, 'ISRM_ID', True)
    return health_data(hia_pop_alloc, incidence_fp, verbose=verbose, race_stratified=False)
  
  
#%% Run Program
if __name__ == "__main__":    
   
    #% Use argparse to parse command line arguments
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
        logging.info('\n ╓────────────────────────────────╖')
        logging.info('║ Beginning Concentration Module ║')
        logging.info('╙────────────────────────────────╜\n')
        
        ## Create emissions and ISRM objects
        # By using a ThreadPoolExecutor to spin up multiple threads, we can
        # read several files at the same time, and even let the population file
        # keep loading while we compute concentrations based on the emissions
        # and ISRM files. Multiple threads are a good option when the bottleneck
        # is disk I/O, not CPU usage.
        logging.info('- Starting multithreaded file reading...')        
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
        logging.info('- Population and ISRM files loaded.')
        
        # Start creating the exp_pop_alloc object.
        exp_pop_alloc_future = executor.submit(
            pop.allocate_population, pop.pop_exp, isrmgrid.geodata, 'ISRM_ID', False)
        executor_jobs.append(exp_pop_alloc_future)
        
        # Creating HIA inputs takes a long time. Start the process as early as 
        # possible, even if the rest of the health module won't start yet.
        if run_health:
            # Starts computing the HIA inputs now, in a parallel process.
            logging.info('- Starting to create hia_inputs in a separate process')
            hia_inputs_future = executor.submit(
                create_hia_inputs, pop, load_file=True, verbose=verbose, geodata=isrmgrid.geodata, incidence_fp=incidence_fp)
            executor_jobs.append(hia_inputs_future)
        
        ## Estimate concentrations
        emis = emis_future.result() # This almost always finishes earlier than the otehr files
        conc = concentration(emis, isrmgrid, detailed_conc_flag, run_calcs=True, verbose=verbose)
        logging.info("<< Generating Concentration Outputs >>")

        ## Create plots and export results
        conc_viz_future = conc.visualize_concentrations_in_background(executor, 'TOTAL_CONC_UG/M3', output_region, output_dir, f_out, ca_shp_path, export=True)
        executor_jobs.append(conc_viz_future)
        conc.export_concentrations(shape_out, f_out)
        logging.info("- Concentration files output into: {}.".format(output_dir))
        
        ## Best practice is to close the pool, but this causes an exception on Mac
        if platform.system() != 'Darwin':
            executor.shutdown(wait=False) # Closes this pool for any more tasks, but returns immediately
        
        ## Perform concentration-related EJ analyses
        exp_pop_alloc = pop.allocate_population(pop.pop_exp, isrmgrid.geodata, 'ISRM_ID', False)
        logging.info('- exp_pop_alloc ready')
        # exp_pop_alloc = pop.allocate_population(pop.pop_exp, isrmgrid.geodata, 'ISRM_ID', False)
        
        # Create the exposure dataframe and run EJ functions
        logging.info('\n << Beginning Exposure EJ Calculations >>')
        exposure_gdf, exposure_pctl, exposure_disparity = run_exposure_calcs(conc, exp_pop_alloc, verbose)    
        if output_exposure:
            export_exposure(exposure_gdf, shape_out, f_out)
        
        # Export results
        plot_percentile_exposure(output_dir, f_out, exposure_pctl, verbose)
        
        ### HEALTH MODULE
        if run_health:
            logging.info('\n ╓────────────────────────────────╖')
            logging.info('║ Beginning Health Impact Module ║')
            logging.info('╙────────────────────────────────╜\n')
            
            # Wait for process creating health input object to finish
            logging.info('- Waiting for hia_inputs to finish')
            hia_inputs = hia_inputs_future.result()
            logging.info('- hia_inputs ready')
            
            # TODO(jdbus): Use the same process pool as above. Why not? 
            with concurrent.futures.ProcessPoolExecutor(max_workers=5) as health_executor:
                logging.info('Using multiprocessing...')
                # Estimate excess mortality
                logging.info('\n << Estimating Excess Mortality for Three Endpoints >>')
                trimmed_conc = conc.detailed_conc_clean[['ISRM_ID','TOTAL_CONC_UG/M3','geometry']]
                pop = hia_inputs.population.groupby('ISRM_ID')[['ASIAN','BLACK','HISLA','INDIG','WHITE','TOTAL']].sum().reset_index()
                
                call_parameters = [
                        (trimmed_conc, hia_inputs.pop_inc, pop, 'ALL CAUSE', krewski),
                        (trimmed_conc, hia_inputs.pop_inc, pop, 'ISCHEMIC HEART DISEASE', krewski),
                        (trimmed_conc, hia_inputs.pop_inc, pop, 'LUNG CANCER', krewski),
                ]
                futures = [health_executor.submit(calculate_excess_mortality, *params, verbose=verbose) for params in call_parameters]
                allcause, ihd, lungcancer = (futures[0].result(), futures[1].result(), futures[2].result())
                
                # Plot and export
                logging.info('\n<< Health Impact Outputs >>')
                futures = []
                for data, title in ([allcause, 'ALL CAUSE'], [ihd, 'ISCHEMIC HEART DISEASE'], [lungcancer, 'LUNG CANCER']):
                    futures.append(health_executor.submit(visualize_and_export_hia, data, ca_shp_path, 'TOTAL', title, output_dir, f_out, shape_out, verbose=verbose))
                
                # TODO(jdbus): I don't think we need to wait for this, why not let the other stuff start right away?
                # The parent process won't terminate until all the children finish.
                logging.info('Waiting for visualizations and exports to complete...')
                concurrent.futures.wait(futures)
                logging.info('Done!')
        
        concurrent.futures.wait(executor_jobs)
        
        # Final log statements
        logging.info('\n ╓────────────────────────────────╖')
        logging.info('║ Success! Run complete.         ║')
        logging.info('╙────────────────────────────────╜\n')
        logging.info('<< ISRM calculations tool has completed all calculations and exports. >>')
        run_time = round((time.time() - start_time)/60.0,0)
        logging.info('- Total run time: {} minutes'.format(run_time))
        logging.info('- All log statements have been saved into a text file: {}'.format(new_logger))
        
        quit()
        