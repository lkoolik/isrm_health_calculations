# ISRM Health Calculations
A repository of scripts used for converting emissions to concentrations and health impacts using the ISRM for California.

*Libby Koolik, UC Berkeley*

Last modified March 23, 2022

## Purpose and Goals
The ultimate goal of this repository is to create a pipeline for estimating disparities in health impacts associated with incremental changes in emissions. Annual average PM<sub>2.5</sub> concentrations are estimated using the [InMAP Source Receptor Matrix](https://www.pnas.org/doi/full/10.1073/pnas.1816102116).

**This is an ongoing effort and is not yet complete.**

## Methodology ##
TBD

## Code Details ##
### Requirements

### Supporting Code
To streamline calculations and increase functionality of the code, python `classes` were created. These class definitions are saved in the `supporting` folder of the repository. The following sections outline how each of these classes work.

#### `emissions.py`
The `emissions` object is primarily built off of `geopandas`. It has the following attributes:

*Metadata Attributes*
* `file_path`: the file path of the raw emissions data
* `file_type`: the type of file being used to provide raw emissions data (for now, only `.shp` is allowed)
* `emissions_name`: a plain English name tied to the emissions data, either provided or automatically generated from the filepath
* `details_to_keep`: any additional details to be preserved throughout the processing (e.g., sector, fuel type). This is not fully built out yet.
* `verbose`: a Boolean indicating whether or not detailed logging statements should be printed
* `units`: units associated with the emissions (e.g., μg/s)
* `valid_file`: a Boolean indicating whether or not the file provided is valid
* `valid_units`: a Boolean indicating whether or not emissions units are compatible with the program
* `valid_emissions`: a Boolean indicating whether or not emissions passed required tests

*Emissions Attributes*
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

asdfasdf
