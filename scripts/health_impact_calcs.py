#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Health Impact Functions

@author: libbykoolik
last modified: 2022-08-10
"""

# Import Libraries
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pyarrow
from scipy.io import netcdf_file as nf
import logging
import os
from os import path
import sys

#%% Health Calculation Helper Functions
def krewski(conc, inc, pop, endpoint):
    ''' 
    Estimates excess mortality from all causes using the Krewski (2009) function 
    
    INPUTS:
        - verbose: a Boolean indicating whether or not detailed logging statements should 
          be printed
        - conc: a float with the exposure concentration for a given geography
        - inc: a float with the background incidence for a given group in a given geography
        - pop: a float with the population estimate for a given group in a given geography
        - endpoint: a string containing either 'ALL CAUSE', 'ISCHEMIC HEART DISEASE', 
          or 'LUNG CANCER'
        
    OUTPUTS: 
        - a float estimating the number of excess mortalities for the `endpoint` across 
          the group in a given geography
    
    '''
    # Beta from Krewski et al (2009) and BenMAP
    beta_dict = {'ALL CAUSE':0.005826891,
                 'ISCHEMIC HEART DISEASE': 0.021511138,
                 'LUNG CANCER': 0.013102826}
    beta = beta_dict[endpoint]
    
    return (1 - (1/np.exp(beta*conc)))*inc*pop

#%% Main Calculation Functions
def calculate_excess_mortality(conc, health_data_obj, endpoint, function, verbose):
    ''' 
    Calculate Excess Mortality 
    
    INPUTS:
        - conc: a float with the exposure concentration for a given geography
        - health_data_obj: a `health_data` object as defined in the `health_data.py` 
          supporting script
        - endpoint: a string containing either 'ALL CAUSE', 'ISCHEMIC HEART DISEASE', 
          or 'LUNG CANCER'
        - function: the health impact function of choice (currently only `krewski` is 
          built out)
        - verbose: a Boolean indicating whether or not detailed logging statements 
          should be printed 
        
    OUTPUTS:
        - pop_inc_conc: a dataframe containing excess mortality for the `endpoint` using
          the `function` provided
    
    '''
    logging.info('- Estimating excess {} mortality from PM2.5. This step may take time.'.format(endpoint.lower()))
    
    # Get the population-incidence  and total concentration
    if verbose:
        logging.info('   - Creating dataframe to combine concentration data with {} mortality BenMAP inputs.'.format(endpoint.lower()))
    conc_hia = conc.detailed_conc_clean[['ISRM_ID','TOTAL_CONC_UG/M3','geometry']].copy()
    pop_inc = health_data_obj.pop_inc.copy().to_crs(conc_hia.crs)
    
    # Merge these on ISRM_ID
    pop_inc_conc = pd.merge(pop_inc, conc_hia[['ISRM_ID','TOTAL_CONC_UG/M3']], on='ISRM_ID')
    
    if verbose:
        logging.info('   - Successfully merged concentrations and {} input data.'.format(endpoint.title()))
        logging.info('   - Estimating {} mortality for each ISRM grid cell.'.format(endpoint.title()))
        
    # Estimate excess mortality
    pop_inc_conc[endpoint] = pop_inc_conc.apply(lambda x: function(x['TOTAL_CONC_UG/M3'],
                                                                   x[endpoint+' INC'],
                                                                   x['POPULATION'],
                                                                   endpoint), axis=1)
    if verbose:
        logging.info('   - Successfully estimated {} mortality for each ISRM grid cell.'.format(endpoint.title()))
        logging.info('   - Determining excess {} mortality by racial/ethnic group.'.format(endpoint.title()))
    
    # Pivot the dataframe to get races as columns
    pop_inc_conc = pop_inc_conc.pivot_table(index='ISRM_ID',columns='RACE', 
                                            values=endpoint, aggfunc='sum', 
                                            fill_value=0)
    if verbose:
        logging.info('   - Performing initial clean up of excess {} mortality data.'.format(endpoint.title()))

    # Add geometry back in
    pop_inc_conc = pd.merge(conc_hia, pop_inc_conc, on='ISRM_ID', how='left')
    pop_inc_conc = pop_inc_conc.fillna(0)
    pop_inc_conc = gpd.GeoDataFrame(pop_inc_conc, geometry='geometry')
    
    # Update column names
    col_rename_dict = {'ASIAN':endpoint+'_ASIAN',
                        'BLACK':endpoint+'_BLACK',
                        'HISLA':endpoint+'_HISLA',
                        'INDIG':endpoint+'_INDIG',
                        'TOTAL':endpoint+'_TOTAL',
                        'WHITE':endpoint+'_WHITE'}
    pop_inc_conc.rename(columns=col_rename_dict, inplace=True)
    
    # Merge the population back in
    if verbose:
        logging.info('   - Adding population data back in for per capita calculations.')
        
    pop = health_data_obj.population.groupby('ISRM_ID')[['ASIAN','BLACK','HISLA','INDIG','WHITE','TOTAL']].sum().reset_index()
    pop_inc_conc = pd.merge(pop_inc_conc, pop, left_on='ISRM_ID', right_on='ISRM_ID', how='left')
    pop_inc_conc = pop_inc_conc.fillna(0)
    
    # Final Clean Up
    if verbose:
        logging.info('   - Performing final clean up.')
        
    pop_inc_conc = pop_inc_conc[['ISRM_ID', 'TOTAL_CONC_UG/M3', 'ASIAN', 'BLACK', 'HISLA',
                                 'INDIG', 'WHITE', 'TOTAL', endpoint+'_ASIAN', endpoint+'_BLACK', 
                                 endpoint+'_HISLA', endpoint+'_INDIG',endpoint+'_TOTAL', 
                                 endpoint+'_WHITE', 'geometry']]
    
    # Print statement
    logging.info('   - {} health impacts calculated.'.format(endpoint.title()))
    
    return pop_inc_conc

#%% Formatting and Exporting Functions
def plot_total_mortality(hia_df, ca_shp_fp, group, endpoint, output_dir, f_out, verbose):
    ''' 
    Plots mortality maps and exports as a png. 
    
    INPUTS:
        - hia_df: a dataframe containing excess mortality for the `endpoint` using the `function`
          provided
        - ca_shp_fp: a filepath string of the California state boundary shapefile
        - group: the racial/ethnic group name
        - endpoint: a string containing either 'ALL CAUSE', 'ISCHEMIC HEART DISEASE', or 
          'LUNG CANCER'
        - output_dir: a filepath string of the location of the output directory
        - f_out: the name of the file output category (will append additional information) 
        - verbose: a Boolean indicating whether or not detailed logging statements should 
          be printed
        
    OUTPUTS:
        - fname: a string filename made by combining the `f_out` with the `group`
          and `endpoint`.
          
    '''
    if verbose:
        logging.info('   - Drawing plot of excess {} mortality from PM2.5 exposure.'.format(endpoint.lower()))
    
    sns.set_theme(context="notebook", style="whitegrid", font_scale=1.25)
    plt.rcParams['patch.linewidth'] = 0
    plt.rcParams['patch.edgecolor'] = 'none'
    plt.rcParams["patch.force_edgecolor"] = False
    
    # Create the output file directory and name string
    fname = f_out + '_' + group + '_' + endpoint + '_excess_mortality.png'
    fname = str.lower(fname)
    fpath = os.path.join(output_dir, fname)
    
    # Read in CA boundary and project hia_df to same coordinates (meters)
    ca_shp = gpd.read_feather(ca_shp_fp)
    hia_df = hia_df.to_crs(ca_shp.crs)
    
    # Clip dataframe to California
    hia_df = gpd.clip(hia_df, ca_shp)
    
    # Determine which column to use
    group = group.upper()
    mortality_col = endpoint + '_' + group
    group_label = group.title()
    
    # Set true zeros to 10^-9 to avoid divide by zero issues
    hia_df.loc[hia_df[group]==0,group] = 10.0**-9.0
    hia_df.loc[hia_df[mortality_col]==0, mortality_col] = 10.0**-9.0

    # Add new columns to hia_df for plotting
    hia_df['POP_AREA_NORM'] = hia_df[group]/hia_df.area*1000.0*1000.0
    hia_df['MORT_AREA_NORM'] = hia_df[mortality_col]/hia_df.area*1000.0*1000.0
    hia_df['MORT_OVER_POP'] = hia_df[mortality_col]/hia_df[group]*100000.0    
    
    # Grab the minimums that do not include the surrogate zeros
    hia_pop_area_min = hia_df.loc[hia_df[group]>10.0**-9.0,'POP_AREA_NORM'].min()
    hia_mort_area_min = hia_df.loc[hia_df[mortality_col]>10.0**-9.0, 'MORT_AREA_NORM'].min()
    
    # Update MORT_OVER_POP to avoid 100% mortality in areas where there is no population
    hia_df.loc[hia_df[group]==hia_df[mortality_col],'MORT_OVER_POP'] = hia_df['MORT_OVER_POP'].min()*0.0001

    # Initialize the figure as three panes
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(1,4, figsize=(22,6))
    
    ## Pane 0: PM2.5 Exposure Concentration
    hia_df.plot(column='POP_AREA_NORM', legend=True,
                legend_kwds={'label':r'Population Density (population/km$^2$)'},
                edgecolor='none', cmap='Greys',
                norm=matplotlib.colors.LogNorm(vmin=hia_pop_area_min,
                                                vmax=hia_df['POP_AREA_NORM'].max()),
                ax=ax0)
    ca_shp.dissolve().plot(edgecolor='black',facecolor='none', linewidth=1,ax=ax0)
    
    ## Pane 1: PM2.5 Exposure Concentration
    hia_df.plot(column='TOTAL_CONC_UG/M3', legend=True,
                legend_kwds={'label':r'PM2.5 Concentration ($\mu$g/m$^3$)'},
                edgecolor='none', cmap='Greys',
                norm=matplotlib.colors.LogNorm(vmin=hia_df['TOTAL_CONC_UG/M3'].min(),
                                                vmax=hia_df['TOTAL_CONC_UG/M3'].max()),
                ax=ax1)
    ca_shp.dissolve().plot(edgecolor='black',facecolor='none', linewidth=1,ax=ax1)
    
    ## Pane 2: Excess Mortality per Area
    hia_df.plot(column='MORT_AREA_NORM', legend=True,
                legend_kwds={'label':r'Excess Mortality (mortality/km$^2$)'},
                edgecolor='none', cmap='Greys',
                norm=matplotlib.colors.LogNorm(vmin=hia_mort_area_min,
                                                vmax=hia_df['MORT_AREA_NORM'].max()),
                ax=ax2)
    ca_shp.dissolve().plot(edgecolor='black',facecolor='none', linewidth=1,ax=ax2)
    
    ## Pane 3: Excess Mortality per Population
    hia_df.plot(column='MORT_OVER_POP',legend=True,
                legend_kwds={'label':r'Mortality per Population (mortality/100 K people)'},
                edgecolor='none', cmap='Greys',
                norm=matplotlib.colors.LogNorm(vmin=hia_df['MORT_OVER_POP'].min(),
                                                vmax=hia_df['MORT_OVER_POP'].max()),
                ax=ax3)
    ca_shp.dissolve().plot(edgecolor='black',facecolor='none', linewidth=1,ax=ax3)

    # Figure Formatting
    ax0.set_title((group_label + ' Population Density').title())
    ax0.xaxis.set_visible(False)
    ax0.yaxis.set_visible(False)

    ax1.set_title((group_label + ' Exposure').title())
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)

    ax2.set_title((group_label + ' ' + endpoint + ' Excess Mortality').title())
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)

    ax3.set_title((group_label + ' ' + endpoint + ' Mortality per 100K').title())
    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)
    fig.tight_layout()
    
    # Export!
    fig.savefig(fpath, dpi=200)
    logging.info('   - Plot of excess {} mortality from PM2.5 exposure output as {}'.format(endpoint.lower(), fname))
    
    return fname

def export_health_impacts(hia_df, group, endpoint, output_dir, f_out, verbose):
    ''' 
    Plots mortality as a shapefile. 
    
    INPUTS:
        - hia_df: a dataframe containing excess mortality for the `endpoint` using the 
          `function` provided
        - group: the racial/ethnic group name
        - endpoint: a string containing either 'ALL CAUSE', 'ISCHEMIC HEART DISEASE', or 
          'LUNG CANCER'
        - output_dir: a filepath string of the location of the output directory
        - f_out: the name of the file output category (will append additional information) 
        - verbose: a Boolean indicating whether or not detailed logging statements should 
          be printed  
        
    OUTPUTS:
        - fname: a string filename made by combining the `f_out` with the `group`
          and `endpoint`.
        
    '''
    if verbose:
        logging.info('   - Exporting excess {} mortality from PM2.5 exposure as a shapefile.'.format(endpoint.lower()))
        
    # Create the output file directory and name string
    fname = f_out + '_' + group + '_' + endpoint + '_excess_mortality.shp'
    fname = str.lower(fname)
    fpath = os.path.join(output_dir, fname)
    
    # Get endpoint shortlabel
    endpoint_labels = {'ALL CAUSE':'ACM_',
                       'ISCHEMIC HEART DISEASE':'IHD_',
                       'LUNG CANCER':'CAN_'}
    l = endpoint_labels[endpoint]
    
    # Update column names
    col_name_dict = {'TOTAL_CONC_UG/M3':'CONC_UG/M3', 
                     'ASIAN':'POP_ASIAN', 
                     'BLACK':'POP_BLACK',
                     'HISLA':'POP_HISLA',
                     'INDIG':'POP_INDIG',
                     'WHITE':'POP_WHITE',
                     'TOTAL':'POP_TOTAL',
                     endpoint+'_ASIAN':l+'ASIAN',
                     endpoint+'_BLACK':l+'BLACK',
                     endpoint+'_HISLA':l+'HISLA',
                     endpoint+'_INDIG':l+'INDIG',
                     endpoint+'_TOTAL':l+'TOTAL',
                     endpoint+'_WHITE':l+'WHITE'}
    
    hia_df.rename(columns=col_name_dict, inplace=True)
    
    # Export
    hia_df.to_file(fpath)
    logging.info('   - Excess {} mortality from PM2.5 exposure output as a shapefile as {}'.format(endpoint.lower(), fname))
    
    return fname

def visualize_and_export_hia(hia_df, ca_shp_fp, group, endpoint, output_dir, f_out, shape_out, verbose):
    ''' 
    Automates this process a bit.
    
    INPUTS:
        - hia_df: a dataframe containing excess mortality for the `endpoint` using the 
          `function` provided
        - ca_shp_fp: a filepath string of the California state boundary shapefile
        - group: the racial/ethnic group name
        - endpoint: a string containing either 'ALL CAUSE', 'ISCHEMIC HEART DISEASE', or 
          'LUNG CANCER'
        - output_dir: a filepath string of the location of the output directory
        - f_out: the name of the file output category (will append additional information) 
        - shape_out: a filepath string for shapefiles
        - verbose: a Boolean indicating whether or not detailed logging statements should 
          be printed      
        
    OUTPUTS:
        - None
    
    '''    
    logging.info('- Visualizing and exporting excess {} mortality.'.format(endpoint.lower()))
    # Plot the map of mortality
    fname = plot_total_mortality(hia_df, ca_shp_fp, group, endpoint, output_dir, f_out, verbose)
    
    # Export the shapefile
    fname = export_health_impacts(hia_df, group, endpoint, shape_out, f_out, verbose)
        
    return #nothing