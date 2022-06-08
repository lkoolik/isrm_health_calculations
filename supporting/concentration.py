#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Concentration Data Object

@author: libbykoolik
last modified: 2022-03-23
"""

# Import Libraries
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
            # Allocate emissions to the ISRM grid
            verboseprint('- Reallocating emissions to the ISRM grid.')
            self.PM25e, self.NH3e, self.VOCe, self.NOXe, self.SOXe = self.process_emissions(self.emissions, self.isrm)
        
            # Estimate concentrations
            verboseprint('- Calculating concentrations of PM25 from each pollutant.')
            self.pPM25 = self.get_concentration(self.PM25e, self.isrm.get_pollutant_layer('PM25'), 0)
            verboseprint('- Concentrations estimated from primary PM2.5.')
            self.pNH4 = self.get_concentration(self.NH3e, self.isrm.get_pollutant_layer('NH3'), 0)
            verboseprint('- Concentrations estimated from ammonia.')
            self.pVOC = self.get_concentration(self.VOCe, self.isrm.get_pollutant_layer('VOC'), 0)
            verboseprint('- Concentrations estimated from organics.')
            self.pNO3 = self.get_concentration(self.NOXe, self.isrm.get_pollutant_layer('NOX'), 0)
            verboseprint('- Concentrations estimated from nitrous oxides.')
            self.pSO4 = self.get_concentration(self.SOXe, self.isrm.get_pollutant_layer('SOX'), 0)
            verboseprint('- Concentrations estimated from sulfur oxides.')
    
            # Add these together at each ISRM grid cell
            self.detailed_conc = self.combine_concentrations(self.pPM25,
                                                             self.pNH4,
                                                             self.pVOC,
                                                             self.pNO3,
                                                             self.pSO4)
            verboseprint('- Detailed concentrations are now ready as .detailed_conc.')
            
            # Create a variable for simplified total concentration in each grid cell
            self.total_conc = self.detailed_conc[['ISRM_ID','geometry','TOTAL_CONC_UG/M3']]
            verboseprint('- Total concentration is now ready as .total_conc.')
            
    def __str__(self):
        return 'Concentration object created from the emissions from '+self.name + ' and the ISRM grid.'

    def __repr__(self):
        return '< Emissions object created from '+self.name + ' and the ISRM grid.>'

    # def buffer_points(self, dist=0.005):
    #     ''' Adds a buffer (in m) to the point type geometries in order to create polygons '''
    #     # First, need to project to coordinates in meters
    #     crs_old = self.emissions.crs
    #     emissions_prj = self.emissions.geometry.copy().to_crs('3310') # California NAD83 Albers (m)
        
    #     # Create buffer of radius dist
    #     emissions_prj['geometry'] = emissions_prj.buffer(dist)
        
    #     # Re-project back to original coordinates
    #     emissions_new_geo = emissions_prj.to_crs(crs_old)
        
    #     # Overwrite point data
        
        
    #     return

    def allocate_emissions(self, emis_layer, isrm_geography):    
        ''' Reallocates the emissions into the ISRM geography using a spatial intersect '''
        
        ## Pre-Process Slightly for Easier Functioning Downstream
        # Deep copy the emissions layer and add an ID field
        emis = emis_layer.copy(deep=True)
        emis['EMIS_ID'] = 'EMIS_'+emis.index.astype(str)
        
        # Re-project the emissions layer into the ISRM coordinate reference system
        emis = emis.to_crs(isrm_geography.crs)
        
        # Store the total emissions from the raw emissions data for later comparison
        old_total = emis['EMISSIONS_UG/S'].sum()
        
        ## Perform Intersect to Reallocate Emissions
        # Get total area of each emissions cell
        emis['area_km2'] = emis.geometry.area/(1000*1000)
        
        # Create intersect object between emis and ISRM grid
        intersect = gpd.overlay(emis, isrm_geography, how='intersection')
        emis_totalarea = intersect.groupby('EMIS_ID').sum()['area_km2'].to_dict()
        
        # Add a total area and area fraction to the intersect object
        intersect['area_total'] = intersect['EMIS_ID'].map(emis_totalarea)
        intersect['area_frac'] = intersect['area_km2'] / intersect['area_total']
        
        # Update the EMISSIONS_UG/S field to scale emissions by the area fraction
        intersect['EMISSIONS_UG/S'] = intersect['area_frac'] * intersect['EMISSIONS_UG/S']  
            
        # Sum over ISRM grid cell
        reallocated_emis = intersect.groupby('ISRM_ID')[['EMISSIONS_UG/S']].sum().reset_index()
        
        ## Preserve all ISRM grid cells for consistent shapes
        reallocated_emis = isrm_geography[['ISRM_ID','geometry']].merge(reallocated_emis,
                                                          how='left',
                                                          left_on='ISRM_ID',
                                                          right_on='ISRM_ID')
        reallocated_emis['EMISSIONS_UG/S'].fillna(0, inplace=True)
        
        ## Confirm that the total has not changed
        assert np.isclose(reallocated_emis['EMISSIONS_UG/S'].sum(), old_total)
        
        return reallocated_emis
    
    def process_emissions(self, emis, isrm_obj):
        pollutants = ['PM25', 'NH3', 'VOC', 'NOX', 'SOX']
        tmp_dct = {}
        for pollutant in pollutants:
            tmp_dct[pollutant] = self.allocate_emissions(emis.get_pollutant_layer(pollutant),
                                                    isrm_obj.geodata)

        return tmp_dct['PM25'], tmp_dct['NH3'], tmp_dct['VOC'], tmp_dct['NOX'], tmp_dct['SOX']
    
    def get_concentration(self, pol_emis, pol_isrm, layer):
        ''' For a given pollutant layer, get the resulting PM25 concentration '''
        # Slice off just the appropriate layer of the ISRM
        pol_isrm_slice = pol_isrm[layer, :, :]
        
        # Concentration is the dot product of emissions and ISRM
        conc = np.dot(pol_emis['EMISSIONS_UG/S'], pol_isrm_slice)
        
        # Convert into a geodataframe
        conc_df = pd.DataFrame(conc, columns=['CONC_UG/M3'], index=pol_emis.index)
        conc_gdf = pol_emis.merge(conc_df, left_index=True, right_index=True)
        
        return conc_gdf
    
    def combine_concentrations(self, pPM25, pNH4, pVOC, pNO3, pSO4):
        ''' Combines concentration from each pollutant into one geodataframe '''
        # Merge to combine into one dataframe
        pol_gdf = pd.merge(pPM25, pNH4, left_on=['ISRM_ID','geometry'], 
                           right_on=['ISRM_ID','geometry'],
                           suffixes=('_PM25','_NH3'))
        
        pol_gdf = pol_gdf.merge(pVOC, left_on=['ISRM_ID','geometry'], 
                           right_on=['ISRM_ID','geometry'],
                           suffixes=('','_VOC'))
                
        pol_gdf = pol_gdf.merge(pNO3, left_on=['ISRM_ID','geometry'], 
                           right_on=['ISRM_ID','geometry'],
                           suffixes=('','_NOX'))
        
                
        pol_gdf = pol_gdf.merge(pSO4, left_on=['ISRM_ID','geometry'], 
                           right_on=['ISRM_ID','geometry'],
                           suffixes=('','_SOX'))
        
        # Quick ugly fix to add the pollutant back onto VOC (otherwise it is dropped)
        pol_gdf.rename(columns={'EMISSIONS_UG/S':'EMISSIONS_UG/S_VOC',
                                'CONC_UG/M3':'CONC_UG/M3_VOC'}, inplace=True)
        
        # Reorder columns for prettiness
        pol_gdf = pol_gdf[['ISRM_ID', 'geometry', 'EMISSIONS_UG/S_PM25',
                           'EMISSIONS_UG/S_NH3', 'EMISSIONS_UG/S_VOC',  
                           'EMISSIONS_UG/S_NOX', 'EMISSIONS_UG/S_SOX', 
                           'CONC_UG/M3_PM25','CONC_UG/M3_NH3', 'CONC_UG/M3_VOC',
                           'CONC_UG/M3_NOX', 'CONC_UG/M3_SOX']]
    
        pol_gdf['TOTAL_CONC_UG/M3'] = pol_gdf['CONC_UG/M3_PM25'] \
                                        + pol_gdf['CONC_UG/M3_NH3'] \
                                        + pol_gdf['CONC_UG/M3_VOC'] \
                                        + pol_gdf['CONC_UG/M3_NOX'] \
                                        + pol_gdf['CONC_UG/M3_SOX']
                                    
        return pol_gdf
    
    def visualize_concentrations(self, var, output_dir, export=False):
        ''' Creates map of concentrations using simple chloropleth '''
        # Note to build this out further at some point in the future, works for now
        # Read in CA boundary
        ca_shp = gpd.read_file('/Users/libbykoolik/Documents/Research/OEHHA Project/scripts/isrm_health_calculations/data/CA_State_TIGER2016.shp')
        ca_prj = ca_shp.to_crs(self.crs)
        
        # Create necessary labels and strings
        if var[0:10] == 'CONC_UG/M3':
            pol = 'Emissions of '+var.split('_')[-1]
        else:
            pol = 'All Emissions'
        t_str = 'PM2.5 Concentrations from {}'.format(pol)
        fname = self.name + '_' + pol + '_concentrations.png'
        fname = str.lower(fname)
        fpath = os.path.join(output_dir, fname)
        
        # Grab relevant layer
        c_to_plot = self.detailed_conc[['ISRM_ID','geometry',var]].copy()
        
        fig, ax = plt.subplots(1,1)
        c_to_plot.plot(column=var,
                              figsize=(20,10),
                              legend=True,
                              legend_kwds={'label':r'Concentration of PM$_{2.5}$ ($\mu$g/m$^3$)'},
                              ax = ax)
        
        ca_prj.plot(edgecolor='black', facecolor='none', ax=ax)
        
        ax.set_title(t_str)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        fig.tight_layout()
        
        if export: 
            fig.savefig(fpath, dpi=200)
        return

    def export_concentrations(self, output_dir, detailed=False,):
        ''' Exports concentration as a shapefile (detailed or total) '''
        
        # If detailed flag is True, export detailed shapefile
        if detailed:
            fname = self.name + '_detailed_concentration.shp' # File Name
            fpath = os.path.join(output_dir, fname)
            
            # Make a copy and change column names to meet shapefile requirements
            gdf_export = self.detailed_conc.copy()
            gdf_export.columns = ['ISRM_ID', 'geometry', 'PM25_UG_S', 'NH3_UG_S',
                                  'VOC_UG_S', 'NOX_UG_S', 'SOX_UG_S', 'fPM_UG_M3', 
                                  'fNH3_UG_M3', 'fVOC_UG_M3', 'fNOX_UG_M3',
                                  'fSOX_UG_M3', 'PM25_UG_M3']
            
            # Export
            gdf_export.to_file(fpath)
            print('- Detailed concentrations output as {}'.format(fname))
            
        # If detailed flag is False, export only total concentration shapefile
        else:
            fname = str.lower(self.name + '_total_concentration.shp') # File Name
            fpath = os.path.join(output_dir, fname)
            
            # Make a copy and change column names to meet shapefile requirements
            gdf_export = self.total_conc.copy()
            gdf_export.columns = ['ISRM_ID', 'geometry', 'PM25_UG_M3']
            
            # Export
            gdf_export.to_file(fpath)
            print('- Total concentratitons output as {}'.format(fname))
        
        return