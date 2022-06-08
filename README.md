# ISRM Health Calculations
A repository of scripts used for converting emissions to concentrations and health impacts using the ISRM for California.

*Libby Koolik, UC Berkeley*

Last modified June 7, 2022

## Table of Contents
* Purpose and Goals ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#purpose-and-goals))
* Methodology ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#methodology))

----

## Purpose and Goals
The ultimate goal of this repository is to create a pipeline for estimating disparities in health impacts associated with incremental changes in emissions. Annual average PM<sub>2.5</sub> concentrations are estimated using the [InMAP Source Receptor Matrix](https://www.pnas.org/doi/full/10.1073/pnas.1816102116).

**This is an ongoing effort and is not yet complete.**

## Methodology ##
TBD

----

## Code Details ##
Below is a brief table of contents for the Code Details section of the Readme.
* Requirements ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#requirements))
* Supporting Code ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#supporting-code))
   * `control_file.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#control_filepy))
   * `emissions.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#emissionspy))

### Requirements

### Supporting Code
To streamline calculations and increase functionality of the code, python `classes` were created. These class definitions are saved in the `supporting` folder of the repository. The following sections outline how each of these classes work.

#### `control_file.py`
The `control_file` object is used to check and read the control file for a run:

*Inputs*
* `file_path`: the file path of the control file

*Attributes*
* `valid_file`: a Boolean indicating whether or not the control file path is valid
* `keywords`: a hardcoded list of the keywords that should be present in the control file
* `blanks_okay`: a hardcoded list of whether each keyword can be blank (based on order of `keywords`)
* `valid_structure`, `no_incorrect_blanks`: Boolean keywords based on internal checks of the control file format
* `run_name`: a string representing the run name preferred by the user
* `emissions_path`: a string representing the path to the emissions input file
* `emissions_units`: a string representing the units of the emissions data
* `check`: a Boolean indicating whether the program should run, or if it should just check the inputs (useful for debugging)
* `verbose`: a Boolean indicating whether the user wants to run in verbose mode

*Simple Functions*
* `get_file_path`: returns the file path

#### `emissions.py`
The `emissions` object is primarily built off of `geopandas`. It has the following attributes:

*Inputs*
* `file_path`: the file path of the raw emissions data
* `units`: units associated with the emissions (e.g., μg/s)
* `name`: a plain English name tied to the emissions data, either provided or automatically generated from the filepath
* `details_to_keep`: any additional details to be preserved throughout the processing (e.g., sector, fuel type) (not fully built out yet)
* `filter_dict`: filters the emissions inputs based on inputted dictionary (not fully built out yet)
* `load_file`: a Boolean indicating whether or not the file should be loaded (for debugging)
* `verbose`: a Boolean indicating whether or not detailed logging statements should be printed

*Emissions Attributes*
* `valid_file`: a Boolean indicating whether or not the file provided is valid
* `valid_units`: a Boolean indicating whether or not emissions units are compatible with the program
* `valid_emissions`: a Boolean indicating whether or not emissions passed required tests
* `file_type`: the type of file being used to provide raw emissions data (for now, only `.shp` is allowed)
* `geometry`: geospatial information associated with the emissions input
* `crs`: the inherent coordinate reference system associated with the emissions input
* `emissions_data`: complete, detailed emissions data from the source
* `emissions_data_clean`: simplified emissions in each grid cell
* `PM25`: primary PM<sub>2.5</sub> emissions in each grid cell
* `NH3`: ammonia emissions in each grid cell
* `VOC`: volatile organic compound emissions in each grid cell
* `NOX`: nitrous oxide emissions in each grid cell
* `SOX`: sulfur oxide emissions in each grid cell

*Simple Functions*
* `get_file_path`: returns the file path
* `get_name`: returns the name associated with the emissions (`emissions_name`)
* `get_unit_conversions`: returns two dictionaries of built-in unit conversions
* `check_path`: uses the `path` library to check if the provided `file_path` exists and if the file is a file
* `check_units`: checks that the provided units are valid against the `get_unit_conversions` dictionaries
* `convert_units`: converts units from provided units to μg/s using the unit dictionaries built-in
* `visualize_emissions`: creates a simple map of emissions for a provided pollutant
* `buffer_emis`: adds a buffer of distance `dist` to the emissions geography
