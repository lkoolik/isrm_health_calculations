# ISRM Health Calculations
A repository of scripts used for converting emissions to concentrations and health impacts using the ISRM for California.

*Libby Koolik, UC Berkeley*

Last modified July 7, 2022

## Table of Contents
* Purpose and Goals ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#purpose-and-goals))
* Methodology ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#methodology))
* Code Details ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#code-details))

----

## Purpose and Goals
The Intervention Model for Air Pollution (InMAP) is a powerful first step towards lowering key technical barriers by making simplifying assumptions that allow for streamlined predictions of PM<sub>2.5</sub> concentrations resulting from emissions-related policies or interventions.\[[*](https://doi.org/10.1371/journal.pone.0176131)\] InMAP performance has been validated against observational data and WRF-Chem, and has been used to perform source attribution and exposure disparity analyses.\[[*](https://doi.org/10.1126/sciadv.abf4491), [*](https://doi.org/10.1073/pnas.1816102116), [*](https://doi.org/10.1073/pnas.1818859116)\] The InMAP Source-Receptor Matrix (ISRM) was developed by running the full InMAP model tens of thousands of times to understand how a unit perturbation of emissions from each grid cell affects concentrations across the grid. However, both InMAP and the ISRM require considerable computational and math proficiency to run and an understanding of various atmospheric science principles to interpret. Furthermore, estimating health impacts requires additional knowledge and calculations beyond InMAP. Thus, a need arises for a standalone and user-friendly process for comparing air quality health disparities associated with various climate change policy scenarios.

The ultimate goal of this repository is to create a pipeline for estimating disparities in health impacts associated with incremental changes in emissions. Annual average PM<sub>2.5</sub> concentrations are estimated using the [InMAP Source Receptor Matrix](https://www.pnas.org/doi/full/10.1073/pnas.1816102116) for California.

## Methodology ##
The ISRM Health Calculation model works by a series of two modules. First, the model estimates annual average change in PM<sub>2.5</sub> concentrations as part of the **Concentration Module**. Second, the excess mortality resulting from the concentration change is calculated in the **Health Module**.

### Concentration Module Methodology ###
The InMAP Source Receptor Matrix (ISRM) links emissions sources to changes in receptor concentrations. There is a matrix layer for each of the five precursor species: primary PM<sub>2.5</sub>, ammonia (NH<sub>3</sub>), oxides of nitrogen (NOx), oxides of sulfur (SOx), and volatile organic compounds (VOC). For each of these species, the ISRM matrix dimensions are: 3 elevations by 21,705 sources by 21,705 receptors. The three elevations of release height within the ISRM are:
* Less than 57 meters
* Between 57 and 140 meters
* Greater than 760 meters.

The units of each cell within the ISRM are micrograms per meter cubed per microgram per second, or concentration per emissions. 

The concentration module has the following steps. Details about the code handling each step are described in the Code Details([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#code-details)) section below.

1. **Preprocessing**: the tool will load the emissions shapefile and perform a series of formatting checks and adjustments. Any updates will be reported through the command line. Additionally, the ISRM layers will be imported as an object. The tool will also identify how many of the ISRM layers are required for concentration calculations.

For each layer triggered in the preprocessing step: 

2. **Emissions Re-Allocation**: the tool will re-grid emissions to the ISRM grid.
   1. The emissions shape and the ISRM shape are intersected.
   2. Emissions for the intersection object are allocated from the original emissions shape by the percent of the original emissions area that is contained within the intersection.
   3. Emissions are summed by ISRM grid cell.
   4. Note: for point source emissions, a small buffer is added to each point to allocate to ISRM grid cells.
3. **Matrix Multiplication**: Once the emissions are re-gridded to the ISRM grid, they are multiplied by the ISRM grid level for the corresponding layer. 

Once all layers are done:

4. **Sum all Concentrations**: concentrations of PM<sub>2.5</sub> are summed by ISRM grid cell.


----

## Code Details ##
Below is a brief table of contents for the Code Details section of the Readme.
* Requirements ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#requirements))
* `isrm_calcs.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#isrm_calcspy))
* Supporting Code ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#supporting-code))
   * `concentration_layer.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#concentration-layerpy))
   * `concentration.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#concentrationpy))
   * `control_file.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#control_filepy))
   * `emissions.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#emissionspy))
   * `health_data.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#health_datapy))
   * `isrm.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#isrmpy))
   * `population.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#populationpy))
* Scripts ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#scripts))
   * `environmental_justice_calcs.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#environmental-justice-calcspy))
   * `health_impact_calcs.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#health-impact-calcspy))
   * `tool_utils.py` ([*](https://github.com/lkoolik/isrm_health_calculations/blob/main/README.md#tool-utilspy))

### Requirements
The code is written in Python 3. The library requirements are included in this repository as `requirements.txt`. For completeness, they are reproduced here:
* attrs==21.4.0
* certifi==2021.10.8
* click==8.1.2
* click-plugins==1.1.1
* cligj==0.7.2
* cycler==0.11.0
* DateTime==4.5
* Fiona==1.8.21
* fonttools==4.32.0
* geopandas==0.10.2
* kiwisolver==1.4.2
* matplotlib==3.5.1
* munch==2.5.0
* numpy==1.22.3
* packaging==21.3
* pandas==1.4.2
* pathlib==1.0.1
* Pillow==9.1.0
* pyarrow==7.0.0
* pyparsing==3.0.8
* pyproj==3.3.0
* python-dateutil==2.8.2
* pytz==2022.1
* Rtree==1.0.0
* scipy==1.8.0
* seaborn==0.11.2
* Shapely==1.8.1.post1
* six==1.16.0
* zope.interface==5.4.0

Python libraries can be installed by running `pip install -r requirements.txt` on a Linux/Mac command line.  

### `isrm_calcs.py`
Text


### Supporting Code
To streamline calculations and increase functionality of the code, python `classes` were created. These class definitions are saved in the `supporting` folder of the repository. The following sections outline how each of these classes work.

#### `concentration_layer.py` 
Text goes here

*Inputs*
text

*Attributes*
text

*Simple Functions*
text

#### `concentration.py` 
Text goes here

*Inputs*
text

*Attributes*
text

*Simple Functions*
text

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

#### `health_data.py` 
Text goes here

*Inputs*
text

*Attributes*
text

*Simple Functions*
text

#### `isrm.py` 
Text goes here

*Inputs*
text

*Attributes*
text

*Simple Functions*
text

#### `population.py` 
Text goes here

*Inputs*
text

*Attributes*
text

*Simple Functions*
text

### Scripts
To streamline calculations and increase functionality of the code, python scripts were created for major calculations/operations. Scripts are saved in the `scripts` folder of the repository. The following sections outline the contents of each script file, and how the functions inside them work.

#### `environmental_justice_calcs.py` 
Text goes here

*Function #1*
text

#### `health_impact_calcs.py` 
Text goes here

*Function #1*
text

#### `tool_utils.py` 
Text goes here

*Function #1*
text
