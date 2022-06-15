#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tool Utils

@author: libbykoolik
last modified: 2022-06-15
"""

#%% Import useful libraries
from pathlib import Path
import sys
import argparse
import os
import datetime
from os import path

#%% ISRM Tool Utils

def create_output_dir(batch, name):
    ''' Creates the output directory for files generated '''
    # Grab current working directory and the 'outputs' sub folder
    parent = os.getcwd()
    sub = 'outputs'
    
    # Output subdirectory will be named with batch and name
    name_list = ['out',batch,name]
    try: name_list.remove('')
    except ValueError: pass
    outdir = '_'.join(name_list)
    
    # If the directory already exists, add an integer to the end
    path_exists_flag = path.exists(os.path.join(parent,sub,outdir))
    # Use while loop in case there are multiple directories already made
    while path_exists_flag:
        if outdir == '_'.join(name_list):
            n = 0
            outdir_tmp = outdir + '_'        
        else:
            # Need to pull just the last two numbers 
            try:
                n = int(outdir.split('_')[-1]) # Grab the number at the end
                outdir_tmp = outdir[:-2] # Get a temporary output directory
                if outdir_tmp[-1] != '_': # Only one run exists in the directory, need to start at 01
                    n = 0
                    outdir_tmp = outdir_tmp + '_'
            except ValueError: #This means there was no number at the end, so start with 0
                outdir_tmp = outdir + '_'
                n = 0
        # Update to the next n
        next_n = str(n+1).zfill(2)
        outdir = outdir_tmp+next_n
        
        # Check if this path exists
        path_exists_flag = path.exists(os.path.join(parent,sub,outdir))
        
    # Add a new variable f_out that adds more information to output file names
    f_out = outdir[4:]# Cuts off the 'out_' from the start
    
    # Make the directory if it does not already exist
    os.mkdir(os.path.join(parent, sub, outdir))
    output_dir = os.path.join(parent,sub,outdir)
    
    # Print a statement to tell user where to look for files
    print("\n<< Output files created will be saved in the following directory: "+output_dir+">>")
    
    return output_dir, f_out