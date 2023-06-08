#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EJ Functions

@author: libbykoolik
last modified: 2023-03-15
"""

# Import Libraries
import pandas as pd
import geopandas as gpd
import logging
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pyarrow
from scipy.io import netcdf_file as nf
import os
from os import path
import sys
sys.path.append('./scripts')
from tool_utils import *

#%%
def create_exposure_df(conc, isrm_pop_alloc, verbose):
    ''' 
    Create an exposure geodataframe from concentration and population.
    
    INPUTS:
        - conc: concentration object
        - isrm_pop_alloc: population object re-allocated to the ISRM grid cell 
          geometry
        - verbose: a Boolean indicating whether or not detailed logging statements 
          should be printed
          
    OUTPUTS:
        - exposure_gdf: a geodataframe with the exposure concentrations and allocated 
          population by racial group
    
    '''
    # Pull the total concentration from the conc object
    conc_gdf = conc.total_conc.copy()
    conc_gdf.columns = ['ISRM_ID', 'geometry', 'PM25_UG_M3']
    
    # Pull only relevant columns from isrm_pop_alloc
    groups = ['TOTAL', 'ASIAN', 'BLACK', 'HISLA', 'INDIG', 'PACIS', 'WHITE','OTHER']
    isrm_pop_alloc = isrm_pop_alloc[['ISRM_ID']+groups].copy()
    
    # Merge concentration and population data based on ISRM_ID
    exposure_gdf = pd.merge(conc_gdf, isrm_pop_alloc, left_on='ISRM_ID', right_on='ISRM_ID')
    
    # Get PWM columns per group
    verboseprint(verbose, '- [EJ] Estimating population weighted mean exposure for each demographic group.')
    for group in groups:
        exposure_gdf = add_pwm_col(exposure_gdf, group)
        
    return exposure_gdf

def add_pwm_col(exposure_gdf, group):
    ''' 
    Adds an intermediate column that multiplies population by exposure.
    
    INPUTS:
        - conc: concentration object from `concentration.py`
        - isrm_pop_alloc: population object (from `population.py`) re-allocated to the ISRM 
          grid cell geometry
        - verbose: a Boolean indicating whether or not detailed logging statements should be 
          printed
          
    OUTPUTS:
        - exposure_gdf: a geodataframe with the exposure concentrations and allocated population
          by racial group
    
    '''
    # Create a string for the PWM column name
    pwm_col = group+'_PWM'
    
    # Create a column for each ISRM cell that is the group total exposure
    exposure_gdf[pwm_col] = exposure_gdf[group]*exposure_gdf['PM25_UG_M3']
    
    return exposure_gdf

def get_pwm(exposure_gdf, group):
    ''' 
    Estimates the population weighted mean exposure for a given group 
    
    INPUTS:
        - exposure_gdf: a geodataframe with the exposure concentrations and allocated population 
          by racial group
        - group: the racial/ethnic group name
        
    OUTPUTS: 
        - PWM_group: the group-level population weighted mean exposure concentration (float)
    '''
    # Create a string for the PWM column name
    pwm_col = group+'_PWM'
    
    # Estimate the total group-level PWM
    PWM_group = exposure_gdf[pwm_col].sum()/exposure_gdf[group].sum()
    
    return PWM_group

def get_overall_disparity(exposure_gdf):
    ''' 
    Returns a table of overall disparity metrics 
    
    INPUTS:
        - exposure_gdf: a geodataframe with the exposure concentrations and allocated population 
          by racial group
        
    OUTPUTS: 
        - pwm_df: a dataframe containing the PWM, absolute disparity, and relative disparity
          of each group
    
    '''
    # Define racial/ethnic groups of interest
    groups = ['TOTAL', 'ASIAN', 'BLACK', 'HISLA', 'INDIG', 'PACIS', 'WHITE','OTHER']
    
    # Create a dataframe to store information
    pwm_df = pd.DataFrame({'Group':groups}, columns=['Group'])
    
    # Use predefined function to get the group PWMs
    pwm_df['Group PWM'] = pwm_df.apply(lambda x: get_pwm(exposure_gdf, x['Group']), axis=1)
    
    # Calculate Absolute and Relative Disparities 
    pwm_df['Absolute Disparity'] = pwm_df['Group PWM'] - pwm_df.loc[0,'Group PWM']
    pwm_df['Relative Disparity'] = pwm_df['Absolute Disparity']/pwm_df.loc[0,'Group PWM']
    
    return pwm_df

def estimate_exposure_percentile(exposure_gdf, verbose):
    ''' 
    Creates a dataframe of percentiles
    
    INPUTS:
        - exposure_gdf: a geodataframe with the exposure concentrations and allocated population 
          by racial group
        - verbose: a Boolean indicating whether or not detailed logging statements should be printed
        
    OUTPUTS:
        - df_pctl: a dataframe of exposure concentrations by percentile of population exposed 
          by group
    
    '''
    if verbose:
        logging.info('- Estimating the exposure level for each percentile of each demographic group population.')
    
    # Define racial/ethnic groups of interest
    groups = ['TOTAL', 'ASIAN', 'BLACK', 'HISLA', 'INDIG', 'PACIS', 'WHITE','OTHER']
    
    # Create a copy to avoid overwriting, then sort based on PM25 concentration
    df_pctl = exposure_gdf.copy()
    df_pctl.sort_values(by='PM25_UG_M3', inplace=True)
    df_pctl.reset_index(drop=True, inplace=True)
    
    # Iterate through each group to estimate the percentile of exposure
    for group in groups:
        # Create a slice of the percentile dataframe
        df_slice = df_pctl[['PM25_UG_M3',group]].copy()
        
        # Add the cumulative sum of the population
        df_slice.loc[:,'Cumulative_Sum_Pop'] = df_slice.loc[:, group].cumsum()
        
        # Estimate the total population in that group, then divide the cumulative sum
        # to get the percentile
        total_pop_group = df_slice[group].sum()
        df_slice.loc[:, 'Percentile_'+group] = df_slice['Cumulative_Sum_Pop']/total_pop_group
        
        # Add the Percentile column into the main percentile dataframe
        df_pctl.loc[:, group] = df_slice.loc[:, 'Percentile_'+group]
    
    return df_pctl

def run_exposure_calcs(conc, pop_alloc, verbose):
    ''' 
    Run the exposure EJ calculations from one script 
    
    INPUTS:
        - conc: concentration object from `concentration.py`
        - isrm_pop_alloc: population object (from `population.py`) re-allocated to the 
          ISRM grid cell geometry
        - verbose: a Boolean indicating whether or not detailed logging statements should
          be printed
        
    OUTPUTS: 
        - exposure_gdf: a dataframe containing the exposure concentrations and population
          estimates for each group
        - exposure_pctl: a dataframe of exposure concentrations by percentile of population
          exposed by group
        - exposure_disparity: a dataframe containing the PWM, absolute disparity, and relative
          disparity of each group
    
    '''
    exposure_gdf = create_exposure_df(conc, pop_alloc, verbose)
    exposure_disparity = get_overall_disparity(exposure_gdf)
    exposure_pctl = estimate_exposure_percentile(exposure_gdf, verbose)
    
    return exposure_gdf, exposure_pctl, exposure_disparity 

def export_exposure(exposure_gdf, output_dir, f_out):
    ''' 
    Exports the exposure_gdf dataframe as a shapefile 
    
    INPUTS:
        - exposure_gdf: a dataframe containing the exposure concentrations and population
          estimates for each group
        - output_dir: a filepath string of the location of the output directory
        - f_out: the name of the file output category (will append additional information)
        
    OUTPUTS:
        - None
    
    '''
    logging.info('- Exporting exposure geodataframe as a shapefile.')
    # If detailed flag is True, export detailed shapefile
        
    fname = str.lower(f_out + '_exposure_concentrations.shp') # File Name
    fpath = os.path.join(output_dir, fname)
    
    # Update the columns slightly
    exposure_gdf = exposure_gdf[['ISRM_ID', 'PM25_UG_M3', 'TOTAL', 'ASIAN', 'BLACK',
                                 'HISLA', 'INDIG', 'PACIS', 'WHITE', 'OTHER',
                                 'geometry']].copy()

    # Export to file
    exposure_gdf.to_file(fpath)
    logging.info('   - [EJ] Exposure concentrations output as {}'.format(fname))

    return

def plot_percentile_exposure(output_dir, f_out, df_pctl, verbose):
    ''' 
    Creates a percentile plot by group 
    
    INPUTS:
        - output_dir: a filepath string of the location of the output directory
        - f_out: the name of the file output category (will append additional information)
        - df_pctl: a dataframe of exposure concentrations by percentile of population
          exposed by group
        - verbose: a Boolean indicating whether or not detailed logging statements should
          be printed
        
    OUTPUTS:
        - None
    
    '''
    verboseprint(verbose, '- [EJ] Drawing plot of exposure by percentile of each racial/ethnic group.')
    # Define racial/ethnic groups of interest
    groups = ['TOTAL', 'ASIAN', 'BLACK', 'HISLA', 'INDIG', 'PACIS', 'WHITE','OTHER']
    
    # Melt the dataframe for easier use of seaborn
    pctl_melt = pd.melt(df_pctl, id_vars='PM25_UG_M3',
                        value_vars=groups,var_name='Racial/Ethnic Group', 
                        value_name='Percentile')
    
    # Adjust formatting for a prettier plot
    pctl_melt['Percentile'] = pctl_melt['Percentile']*100
    rename_dict = {'TOTAL':'Total', 'ASIAN':'Asian','BLACK':'Black',
                              'HISLA':'Hispanic/Latino', 'INDIG':'Indigenous', 
                              'PACIS':'Pacific Islander', 'WHITE':'White', 
                              'OTHER':'Other'}
    pctl_melt['Racial/Ethnic Group'] = pctl_melt['Racial/Ethnic Group'].map(rename_dict)
    sns.set_theme(context="notebook", style="whitegrid", font_scale=1.75)

    # Initialize the figure
    fig, ax = plt.subplots(figsize=(10,8))
    sns.lineplot(data=pctl_melt, x='Percentile', y='PM25_UG_M3', hue='Racial/Ethnic Group', ci=None, 
                 linewidth=3, palette='deep', ax=ax)
    ax.set(ylabel='PM$_{2.5}$ Exposure ($\mu$g/m$^3$)')
    ax.set_xticks(ticks=[5,25,50,75,95], 
                  labels=['5th','25th','50th','75th','95th'])
    
    # Save the file
    fname =f_out+'_PM25_Exposure_Percentiles.png' # File Name
    fpath = os.path.join(output_dir, fname)
    fig.savefig(fpath, dpi=200)
    logging.info('- [EJ] Exposure concentration by percentile figure output as {}'.format(fname))
    
    return