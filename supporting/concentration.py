#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Total Concentration Data Object

@author: libbykoolik
last modified: 2022-06-09
"""

# Import Libraries
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file as nf
import os
from os import path
import sys
sys.path.append('/Users/libbykoolik/Documents/Research/OEHHA Project/scripts/isrm_health_calculations/supporting')
from isrm import isrm
from emissions import emissions
from concentration_layer import concentration_layer

#%% Define the Concentration Object
class concentration:
    '''
    Defines a new object for storing and manipulating concentration data.
    
    INPUTS:
        - tbd
        
    '''
    def __init__(self, emis_obj, isrm_obj, run_calcs=True, verbose=False):
        ''' Initializes the Concentration object'''        
        print('\nEstimating concentrations from provided emissions using the ISRM.')
        # Initialize concentration object by reading in the emissions and isrm 
        self.emissions = emis_obj
        self.isrm = isrm_obj
        
        # Get a few other metadata
        self.isrm_id = self.isrm.ISRM_ID
        self.isrm_geom = self.isrm.geometry
        self.crs = self.isrm.crs
        self.name = self.emissions.emissions_name
        self.verbose = verbose
        verboseprint = print if self.verbose else lambda *a, **k:None # for logging
        verboseprint('\nCreating a new concentration object')
                
        # Run concentration calculations
        if run_calcs:
            self.detailed_conc, self.detailed_conc_clean, self.total_conc = self.combine_concentrations()
            verboseprint('\nTotal concentrations are now ready.')
            
    def __str__(self):
        return 'Concentration object created from the emissions from '+self.name + ' and the ISRM grid.'

    def __repr__(self):
        return '< Emissions object created from '+self.name + ' and the ISRM grid.>'

    def run_layer(self, layer):
        ''' Estimates concentratiton for a single layer '''
        # Creates a concentration_layer object for the given layer
        conc_layer = concentration_layer(self.emissions, self.isrm, layer, run_calcs=True, verbose=self.verbose)
        
        # Copies out just the detailed_conc object and adds the LAYER column
        detailed_conc_layer = conc_layer.detailed_conc.copy()
        detailed_conc_layer['LAYER'] = layer
        
        return detailed_conc_layer
 
    
    def combine_concentrations(self):
        ''' 
        Creates a concentration_layer object for each valid layer and then 
        combines them all into three sets of concentration data
        '''
        # Define a concentration layer list for easier appending
        conc_layers = []
        
        # Run each layer if the layer flag is True
        if self.emissions.L0_flag: conc_layers.append(self.run_layer(0))
        if self.emissions.L1_flag: conc_layers.append(self.run_layer(1))
        if self.emissions.L2_flag: conc_layers.append(self.run_layer(2))
        
        # Concatenate these detailed concentration dataframes
        detailed_concentration = pd.concat(conc_layers)
        
        ## Sum each concentration field across ISRM_ID
        # First, need to get rid of unnecessary columns
        detailed_concentration_clean = detailed_concentration[detailed_concentration.columns.drop(list(detailed_concentration.filter(regex='EMISSIONS')))]
        detailed_concentration_clean = detailed_concentration_clean.drop(columns='geometry').copy()
        
        # Add across ISRM IDs
        detailed_concentration_clean = detailed_concentration_clean.groupby(['ISRM_ID']).sum().reset_index()
        
        # Merge back in the geodata
        geodata = self.isrm.geodata.copy()
        detailed_concentration_clean = pd.merge(detailed_concentration_clean, geodata, 
                                                left_on='ISRM_ID', right_on='ISRM_ID')
        
        # Make a final version that is very simple
        total_concentration = detailed_concentration_clean[['ISRM_ID','geometry', 'TOTAL_CONC_UG/M3']].copy()
        
        return detailed_concentration, detailed_concentration_clean, total_concentration
    
    
    def visualize_concentrations(self, var, output_region, output_dir, f_out, ca_shp_fp, export=False):
        ''' Creates map of concentrations using simple chloropleth '''
        # Note to build this out further at some point in the future, works for now
        # Read in CA boundary
        ca_shp = gpd.read_feather(ca_shp_fp)
        ca_prj = ca_shp.to_crs(self.crs)
        
        # Reproject output_region
        output_region = output_region.to_crs(self.crs)
        
        # Create necessary labels and strings
        if var[0:10] == 'CONC_UG/M3':
            pol = 'Emissions of '+var.split('_')[-1]
        else:
            pol = 'All Emissions'
        t_str = 'PM2.5 Concentrations from {}'.format(pol)
        fname = f_out + '_' + pol + '_concentrations.png'
        fname = str.lower(fname)
        fpath = os.path.join(output_dir, fname)
        
        # Grab relevant layer
        c_to_plot = self.detailed_conc_clean[['ISRM_ID','geometry',var]].copy()
        
        # Clip to output region
        c_to_plot = gpd.clip(c_to_plot, output_region)
        
        fig, ax = plt.subplots(1,1)
        c_to_plot.plot(column=var,
                              figsize=(20,10),
                              legend=True,
                              legend_kwds={'label':r'Concentration of PM$_{2.5}$ ($\mu$g/m$^3$)'},
                              ax = ax)
        
        ca_prj.plot(edgecolor='black', facecolor='none', ax=ax)
        
        # Clip to output_region
        minx, miny, maxx, maxy = output_region.total_bounds
        ax.set_xlim(minx, maxx)
        ax.set_ylim(miny, maxy)
        
        ax.set_title(t_str)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        fig.tight_layout()
        
        if export: 
            fig.savefig(fpath, dpi=200)
        return

    def export_concentrations(self, output_dir, f_out, detailed=False,):
        ''' Exports concentration as a shapefile (detailed or total) '''
        
        # If detailed flag is True, export detailed shapefile
        if detailed:
            fname = f_out + '_detailed_concentration.shp' # File Name
            fpath = os.path.join(output_dir, fname)
            
            # Make a copy and change column names to meet shapefile requirements
            gdf_export = self.detailed_conc.copy()
            gdf_export.columns = ['ISRM_ID', 'geometry', 'PM25_UG_S', 'NH3_UG_S',
                                  'VOC_UG_S', 'NOX_UG_S', 'SOX_UG_S', 'fPM_UG_M3', 
                                  'fNH3_UG_M3', 'fVOC_UG_M3', 'fNOX_UG_M3',
                                  'fSOX_UG_M3', 'PM25_UG_M3', 'LAYER']
            
            # Export
            gdf_export.to_file(fpath)
            print('- Detailed concentrations output as {}'.format(fname))
            
        # If detailed flag is False, export only total concentration shapefile
        else:
            fname = str.lower(f_out + '_total_concentration.shp') # File Name
            fpath = os.path.join(output_dir, fname)
            
            # Make a copy and change column names to meet shapefile requirements
            gdf_export = self.total_conc.copy()
            gdf_export.columns = ['ISRM_ID', 'geometry', 'PM25_UG_M3']
            
            # Export
            gdf_export.to_file(fpath)
            print('- Total concentrations output as {}'.format(fname))
        
        return