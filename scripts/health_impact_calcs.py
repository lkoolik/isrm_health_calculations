#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Health Impact Functions

@author: libbykoolik
last modified: 2022-06-29
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
import os
from os import path
import sys

#%% Health Calculation Helper Functions
def krewski(start_age, conc, inc, pop):
    ''' Estimates excess mortality from all causes using the Krewski (2009) function '''
    # Beta from Krewski et al (2009) and BenMAP
    beta = 0.005827
    
    # Only estimate mortality for people older than 30 years old
    if start_age >= 30:
        return (1 - (1/np.exp(beta*conc)))*inc*pop
    else:
        return 0
    
#%% Main Calculation Functions
def intersect_conc_hia(conc_int, pop_inc_int):
    ''' Create an intersection between concentration data and health data '''
    ## Store county population to make sure it doesn't change after intersection
    pop_totals = pop_inc_int.groupby('NAME')[['ASIAN', 'BLACK', 'HISLA', 'INDIG', 'WHITE', 'TOTAL']].sum()

    ## Add the land area as a feature of the dataframes
    pop_inc_int['POP_AREA_M2'] = pop_inc_int.geometry.area/(1000.0*1000.0)
    conc_int['CONC_AREA_M2'] = conc_int.geometry.area/(1000.0*1000.0)

    ## Create an intersect object
    hia_intersect = gpd.overlay(pop_inc_int, conc_int, how='intersection')

    # Allocate using area fractions
    hia_intersect['INT_AREA_M2'] = hia_intersect.geometry.area/(1000.0*1000.0)
    hia_intersect['FRAC'] = hia_intersect['INT_AREA_M2'] / hia_intersect['POP_AREA_M2']

    for col in ['ASIAN', 'BLACK', 'HISLA', 'INDIG', 'WHITE', 'TOTAL']:
        hia_intersect[col] = hia_intersect[col] * hia_intersect['FRAC']
        
    intersect_pop_totals = hia_intersect.groupby('NAME')[['ASIAN', 'BLACK', 'HISLA', 'INDIG', 'WHITE', 'TOTAL']].sum()
    
    return hia_intersect

def calculate_excess_mortality(conc, hia_inputs, function=krewski):
    ''' Description '''
    # Make copies of the concentration layer and population-incidence layer and convert to same CRS
    conc_int = conc.detailed_conc_clean[['ISRM_ID','TOTAL_CONC_UG/M3','geometry']].copy()
    pop_inc_int = hia_inputs.pop_inc.copy().to_crs(conc_int.crs)
    
    # Store the geodata for later
    geo_data = pop_inc_int[['NAME','geometry']].drop_duplicates()
    
    # Perform an intersection on the two datasets
    hia_intersect = intersect_conc_hia(conc_int, pop_inc_int)
    
    # Call the function for each racial group
    for col in ['ASIAN', 'BLACK', 'HISLA', 'INDIG', 'WHITE', 'TOTAL']:
        new_col = 'MORTALITY_'+col
        hia_intersect[new_col] = hia_intersect.apply(lambda x: function(x['START_AGE'],
                                                                        x['TOTAL_CONC_UG/M3'],
                                                                        x['INCIDENCE'],
                                                                        x[col]),
                                                     axis=1)

    ## Clean up the dataframe
    hia_df = hia_intersect[['NAME', 'STATE_NAME','ASIAN','BLACK', 'HISLA', 'INDIG', 
                        'WHITE', 'TOTAL','MORTALITY_ASIAN','MORTALITY_BLACK', 
                        'MORTALITY_HISLA','MORTALITY_INDIG','MORTALITY_WHITE', 
                        'MORTALITY_TOTAL']]

    hia_df = hia_df.groupby(['NAME', 'STATE_NAME'])[['ASIAN','BLACK', 'HISLA', 'INDIG', 
                        'WHITE', 'TOTAL','MORTALITY_ASIAN', 'MORTALITY_BLACK', 
                        'MORTALITY_HISLA','MORTALITY_INDIG', 'MORTALITY_WHITE',
                        'MORTALITY_TOTAL']].sum().reset_index()

    ## Add back in geometry
    hia_df = pd.merge(geo_data, hia_df, left_on='NAME', right_on='NAME')
    
    return hia_df

def plot_total_mortality(hia_df, ca_shp_fp, group, output_dir, f_out, export):
    ''' 
    Plots mortality maps and exports as a png. 
    INPUTS:
        - hia_df: dataframe of health impacts
        - ca_shp_fp: california shapefile filepath
        - group: racial/ethnic group
        - output_dir: directory to output the plot into
        - f_out: filename
    '''
    sns.set_theme(context="notebook", style="whitegrid", font_scale=1.25)
    
    # Create the output file directory and name string
    fname = f_out + '_' + group + '_excess_mortality.png'
    fname = str.lower(fname)
    fpath = os.path.join(output_dir, fname)
    
    # Read in CA boundary and project hia_df to same coordinates (meters)
    ca_shp = gpd.read_feather(ca_shp_fp)
    hia_df = hia_df.to_crs(ca_shp.crs)
    
    # Determine which column to use
    group = group.upper()
    mortality_col = 'MORTALITY_'+group
    group_label = group.title()
    
    # Add new columns to hia_df for plotting
    hia_df['POP_AREA_NORM'] = hia_df[group]/hia_df.area*1000.0*1000.0
    hia_df['MORT_AREA_NORM'] = hia_df[mortality_col]/hia_df.area*1000.0*1000.0
    hia_df['MORT_OVER_POP'] = hia_df[mortality_col]/hia_df[group]*100000.0    

    # Initialize the figure as three panes
    fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(18,6))
    
    ## Pane 1: Total Population per Area
    hia_df.plot(column='POP_AREA_NORM', legend=True,
                legend_kwds={'label':r'Total Population (people/km$^2$)'},
                norm=matplotlib.colors.LogNorm(vmin=hia_df['POP_AREA_NORM'].min(),
                                                vmax=hia_df['POP_AREA_NORM'].max()),
                edgecolor='none', cmap='viridis',
                ax=ax1)
    ca_shp.dissolve().plot(edgecolor='black',facecolor='none', linewidth=1,ax=ax1)
    
    ## Pane 2: Excess Mortality per Area
    hia_df.plot(column='MORT_AREA_NORM', legend=True,
                legend_kwds={'label':r'Excess Mortality (mortality/km$^2$)'},
                norm=matplotlib.colors.LogNorm(vmin=hia_df['MORT_AREA_NORM'].min(),
                                                vmax=hia_df['MORT_AREA_NORM'].max()),
                edgecolor='none', cmap='viridis',
                ax=ax2)
    ca_shp.dissolve().plot(edgecolor='black',facecolor='none', linewidth=1,ax=ax2)
    
    ## Pane 3: Excess Mortality per Population
    hia_df.plot(column='MORT_OVER_POP',legend=True,
                legend_kwds={'label':r'Mortality per Population (mortality/100 K people)'},
                edgecolor='none', cmap='viridis',
                ax=ax3)
    ca_shp.dissolve().plot(edgecolor='black',facecolor='none', linewidth=1,ax=ax3)

    # Figure Formatting
    ax1.set_title(group_label + ' Population')
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)

    ax2.set_title(group_label + ' Excess Mortality')
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)

    ax3.set_title(group_label + ' Mortality per 100K')
    ax3.xaxis.set_visible(False)
    ax3.yaxis.set_visible(False)
    fig.tight_layout()
    
    ## Export 
    if export: 
        fig.savefig(fpath, dpi=200)
    return