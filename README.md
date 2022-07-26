# ISRM Health Calculations
A repository of scripts used for converting emissions to concentrations and health impacts using the ISRM for California.

*Libby Koolik, UC Berkeley*

Last modified July 18, 2022

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
The `isrm_calcs.py` script is the main script file that drives the tool. This script operates the command line functionality, defines the health impact calculation objects, calls each of the supporting functions, and outputs the desired files. The `isrm_calcs.py` script is not split into functions or objects, instead, it is run through two sections: (1) Initialization and (2) Run Program.

#### Initialization
In the initialization section of `isrm_calcs.py`, the parser object is created in order to interface with the command line. The parser object is created using the `argparse` library. Currently, the only arguments accepted by the parser object are `-i` for input file and `-h` for help. 

Once the parser is defined, the control file object is created using `control_file.py` class object. A number of metadata variables are defined from the control file. 

Next, a number of internally saved data file paths are saved.

Finally, the output_region is defined based on the `get_output_region` function defined in `tool_utils.py`. The output region is then stored for use in later functions.

#### Run Program
The run program section of the code is split into two modes. If the CHECK_INPUTS flag is given, the tool will run in check mode, where it will check that each of the inputs is valid and then quit. If the CHECK_INPUTS flag is not given, the tool will run the full program. 

It will start by creating an output directory using the `create_output_dir` function from `tool_utils.py`. It will also create a shapefile subdirectory within the output folder directory using `create_shape_out`. 

Then, the tool will begin the concentration module. This starts by defining an emissions object and an isrm object using the `emissions.py` and `isrm.py` supporting class objects. The concentrations will be estimated using the `concentration.py` object, which relies on the `concentration_layer.py` object. The concentrations will then be output as a map of total exposure concentration and a shapefile with detailed exposure information.

Next, the tool will run environmental justice exposure calculations using the `create_exposure_df`, `get_overall_disparity`, and `estimate_exposure_percentile` functions from the `environmental_justice_calcs.py` file. The exposure percentiles will then be plotted and exported using the `plot_percentile_exposure` function.

Finally, if indicated by the user, the tool will begin the health module. It will create the health input object using the `health_data.py` library and then estimate the three endpoints of excess mortality using `calculate_excess_mortality` from the `health_impact_calcs` file. Each endpoint will then be mapped and exported using `visualize_and_export_hia`.


### Supporting Code
To streamline calculations and increase functionality of the code, python `classes` were created. These class definitions are saved in the `supporting` folder of the repository. The following sections outline how each of these classes work.

#### `concentration_layer.py` 
The `concentration_layer` object runs ISRM-based calculations using a single vertical layer of the ISRM grid. The object inputs an emissions object (from `emissions.py`), the ISRM object (from `isrm.py`), and the layer number corresponding to the vertical layer of the ISRM grid. The object then estimates concentrations at ground-level resulting from emissions at that vertical layer release range.

*Inputs*
* `emis_obj`: the emissions object, as defined by `emissions.py`
* `isrm_obj`: the ISRM object, as defined by `isrm.py`
* `layer`: the layer number (0, 1, or 2)

*Attributes*
* `isrm_id`: a Series of all ISRM grid cell IDs
* `receptor_id`: a Series of all receptor IDs
* `isrm_geom`: the geometry (geographic attributes) of the ISRM grid
* `crs`: the coordinate reference system associated with the ISRM grid
* `name`: a string representing the run name preferred by the user
* `check`: a Boolean indicating whether the program should run, or if it should just check the inputs (useful for debugging)
* `verbose`: a Boolean indicating whether the user wants to run in verbose mode

*Calculated Attributes*
* `PM25e`, `NH3e`, `VOCe`, `NOXe`, `SOXe`: geodataframes of the emissions (for each pollutant) from that layer re-allocated onto the ISRM grid
* `pPM25`, `pNH4`, `pVOC`, `pNO3`, `pSO4`: geodataframes of the concentrations from each primary pollutant from the emissions of that pollutant in that `layer`
* `detailed_conc`: geodataframe containing columns for each primary pollutant's contribution to the total ground-level PM<sub>2.5</sub> concentrations

*Simple Functions*
* `allocate_emissions`: inputs the emissions layer and the ISRM geography, and re-allocates the emissions to the ISRM geography using an area-based allocation procedure
* `cut_emissions`: inputs the pollutant geodataframe from the emissions object and slices it based on the minimum and maximum release heights (minimum inclusive, maximum exclusive) associated with the ISRM vertical layer
* `process_emissions`: for each of the five primary pollutants, runs `cut_emissions` and then `allocate_emissions` to return the geodataframes of emissions of each primary pollutant released in the `layer` allocated to the ISRM grid
* `get_concentration`: for a pollutant's emission layer (`POLe`), the ISRM matrix for that pollutant, and the `layer` ID, estimates the concentration at ground-level for the primary pollutant (`pPOL`)
* `combine_concentrations`: merges together all five of the primary pollutant concentration geodataframes (`pPOL`) and adds them together to get total ground-level concentrations resulting from emissions released in that `layer`

#### `concentration.py` 
The `concentration` object runs ISRM-based calculations for each of the vertical layer's of the ISRM grid by processing individual `concentration_layer` objects. The object inputs an emissions object (from `emissions.py`) and the ISRM object (from `isrm.py`). The object then estimates total concentrations at ground-level resulting from emissions.

*Inputs*
* `emis_obj`: the emissions object, as defined by `emissions.py`
* `isrm_obj`: the ISRM object, as defined by `isrm.py`

*Attributes*
* `isrm_id`: a Series of all ISRM grid cell IDs
* `isrm_geom`: the geometry (geographic attributes) of the ISRM grid
* `crs`: the coordinate reference system associated with the ISRM grid
* `name`: a string representing the run name preferred by the user
* `check`: a Boolean indicating whether the program should run, or if it should just check the inputs (useful for debugging)
* `verbose`: a Boolean indicating whether the user wants to run in verbose mode

*Calculated Attributes*
* `detailed_conc`: geodataframe of the detailed concentrations at ground-level combined from all three vertical layers
* `detailed_conc_clean`: simplified geodataframe of the detailed concentrations at ground-level combined from all three vertical layers
* `total_conc`: geodataframe with total ground-level PM<sub>2.5</sub> concentrations across the ISRM grid

*Simple Functions*
* `run_layer`: estimates concentrations for a single layer by creating a `concentration_layer` object for that layer
* `combine_concentrations`: checks for each of the layer flags in the `emissions` object, and then calls the `run_layer` function for each layer that is flagged. Then, combines the concentrations from each layer flagged into the three concentration geodataframes described above
* `visualize_concentrations`: draws a map of concentrations for a variable (`var`) and exports it as a PNG into an output directory (`output_dir`) of choice
* `export_concentrations`: exports concentrations as a shapefile into an output directory (`output_dir`) of choice

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

*Attributes*
* `valid_file`: a Boolean indicating whether or not the file provided is valid
* `valid_units`: a Boolean indicating whether or not emissions units are compatible with the program
* `valid_emissions`: a Boolean indicating whether or not emissions passed required tests
* `file_type`: the type of file being used to provide raw emissions data (for now, only `.shp` is allowed)
* `geometry`: geospatial information associated with the emissions input
* `crs`: the inherent coordinate reference system associated with the emissions input
* `emissions_data`: complete, detailed emissions data from the source
* `emissions_data_clean`: simplified emissions in each grid cell

*Calculated Attributes*
* `PM25`: primary PM<sub>2.5</sub> emissions in each grid cell
* `NH3`: ammonia emissions in each grid cell
* `VOC`: VOC compound emissions in each grid cell
* `NOX`: NOx emissions in each grid cell
* `SOX`: SOx emissions in each grid cell
* `L0_flag`, `L1_flag`, `L2_flag`, `linear_interp_flag`: Booleans indicating whether each layer should be calculated based on emissions release heights

*Simple Functions*
* `get_file_path`: returns the file path
* `get_name`: returns the name associated with the emissions (`emissions_name`)
* `get_unit_conversions`: returns two dictionaries of built-in unit conversions
* `check_path`: uses the `path` library to check if the provided `file_path` exists and if the file is a file
* `check_units`: checks that the provided units are valid against the `get_unit_conversions` dictionaries
* `load_emissions`: detects the filetype of the emissions file and calls the appropriate load function
* `load_shp`: loads the emissions data from a shapefile
* `check_height`: checks that the height column is present in the emissions file; if not, assumes emissions are released at ground-level
* `check_emissions`: runs a number of checks on the emissions data to ensure data are valid before running anything
* `map_pollutant_names`: replaces pollutant names if they are not found in the emissions data based on near-misses (e.g., PM2.5 for PM25)
* `filter_emissions`: filters the emissions based on the `filter_dict` input
* `check_geo_types`: checks what geometries are present in the emissions shapefile (e.g., points, polygons, multipolygons); if points exist, uses `buffer_emis` to convert to polygons
* `buffer_emis` converts points to polygons by adding a buffer of `dist` 
* `clean_up`: simplifies the emissions data by removing unnecessary dimensions, converting units as appropriate, and updating the column names
* `convert_units`: converts units from provided units to μg/s using the unit dictionaries built-in
* `split_polutants`: converts the emissions layer into separate objects for each pollutant
* `which_layers`: determines the `L0_flag`, `L1_flag`, `L2_flag`, and `linear_interp_flag` variables based on the HEIGHT column of the emissions data
* `visualize_emissions`: creates a simple map of emissions for a provided pollutant
* `get_pollutant_layer`: pulls a single pollutant layer based on `pol_name`

#### `health_data.py` 
The `health_data` object stores and manipulates built-in health data (population and incidence rates) from BenMAP. It inputs a dictionary of filepaths and two Boolean run options (`verbose` and `race_stratified`) to return dataframes of population, incidence, and combined population-incidence information (`pop_inc`).

*Inputs*
* `filepath_dict`: a dictionary with the filepaths of each input feather file, which must have "POPULATION" and "INCIDENCE" as keys
* `verbose`: a Boolean indicating whether or not detailed logging statements should be printed
* `race_stratified`: a Boolean indicating whether race-stratified incidence rates should be used

*Attributes*
* `population_source`: string containing the source of the population data
* `incidence_source`: string containing the source of the incidence data

*Calculated Attributes*
* `population`: a geodataframe containing the raw population data from BenMAP
* `incidence`: a geodataframe containing the raw incidence data from BenMAP
* `pop_inc`: a geodataframe containing the combined population and incidence data based on the requested geographies

*Simple Functions*
* `load_data`: reads in the population and incidence data from feather files
* `update_pop`: updates the population dataset by melting (unpivot) and renaming columns
* `update_inc`: updates the incidence dataset by pivoting columns around endpoints and renaming columns
* `get_incidence_lookup`: creates a small incidence lookup table based on the name and age ranges
* `get_incidence_pop`: helper function that returns the incidence for a given name, race, age range, and endpoint
* `make_incidence_lookup`: creates a lookup dictionary using the `get_incidence_pop` function for each endpoint
* `incidence_by_age`: creates a smaller incidence table for merging by calling `get_incidence_lookup` for each endpoint
* `combine_pop_inc`: creates the `pop_inc` dataframe by doing a spatial merge on the population and incidence data and then using lookup tables to determine the appropriate values

#### `isrm.py` 
The `isrm` object loads, stores, and manipulates the ISRM grid data. 

*Inputs*
* `isrm_fps`: a list of filepath strings for the NH3, NOx, PM25, SOX, and VOC paths, respectively
* `isrm_gfp`: a filepath string for the geometry feather of the ISRM grid
* `output_region`: a geodataframe of the region for results to be output, as calculated by `get_output_region` in `tool_utils.py`
* `region_of_interest`: the name of the region contained in the `output_region`
* `load_file`: a Boolean indicating whether or not the file should be loaded (for debugging)
* `verbose`: a Boolean indicating whether or not detailed logging statements should be printed

*Attributes*
* `nh3_path`, `nox_path`, `pm25_path`, `sox_path`, `voc_path`: the filepath strings for each of the primary pollutant ISRM variables
* `valid_file`: a Boolean indicating whether or not the file provided is valid
* `valid_geo_file`: a Boolean indicating whether the ISRM geometry file provided is valid
* `geodata`: a geodataframe containing the ISRM feather file information
* `crs`: the inherent coordinate reference system associated with the ISRM geometry
* `geometry`: geospatial information associated with the ISRM geometry

*Calculated Attributes*
* `receptor_IDs`: the IDs associated with ISRM receptors within the `output_region`
* `receptor_geometry`: the geospatial information associated with the ISRM receptors within the `output_region`
* `PM25`, `NH3`, `NOx`, `SOX`, `VOC`: the ISRM matrices for each of the primary pollutants

*Simple Functions*
* `check_path`: checks if the files exist at the paths specified (both data and geo files)
* `load_and_cut`: loads the numpy layers for a pollutant and trims the columns of each vertical layer's matrix to only include the `receptor_IDs` within the `output_region`
* `load_geodata`: loads the feather file into a geopandas dataframe
* `clip_isrm`: clips the ISRM receptors to only the relevant ones based on the `output_region` (i.e., returns the `receptor_IDs` and `receptor_geometry` objects)
* `get_pollutant_layer`: returns the ISRM matrix for a single pollutant
* `map_isrm`: simple function for mapping the ISRM grid cells

#### `population.py` 
The `population` object stores detailed Census tract-level population data for the environmental justice exposure calculations.

*Inputs*
* `file_path`: the file path of the raw population data
* `load_file`: a Boolean indicating whether or not the file should be loaded (for debugging)
* `verbose`: a Boolean indicating whether or not detailed logging statements should be printed

*Attributes*
* `valid_file`: a Boolean indicating whether or not the file provided is valid
* `geometry`: geospatial information associated with the emissions input
* `pop_data`: complete, detailed population data from the source
* `crs`: the inherent coordinate reference system associated with the emissions input
* `pop_gdf`: a geodataframe containing the population information with associated spatial information

*Simple Functions*
* `check_path`: checks to see if the file exists at the path specified and returns whether the file is valid
* `load_population`: loads the population data based on the file extension
* `load_feather`: loads the population feather data using geopandas and post-processes
* `allocate_population`: reallocates population into new geometry using a spatial intersect

### Scripts
To streamline calculations and increase functionality of the code, python scripts were created for major calculations/operations. Scripts are saved in the `scripts` folder of the repository. The following sections outline the contents of each script file, and how the functions inside them work.

#### `environmental_justice_calcs.py` 
The `environmental_justice_calcs` script file contains a number of functions that help calculate exposure metrics for environmental justice analyses.

1. `create_exposure_df`
   1. Inputs:
      * `conc`: concentration object from `concentration.py`
      * `isrm_pop_alloc`: population object (from `population.py`) re-allocated to the ISRM grid cell geometry
      * `verbose`: a Boolean indicating whether or not detailed logging statements should be printed
   2. Outputs
      * `exposure_gdf`: a geodataframe with the exposure concentrations and allocated population by racial group
   3. Methodology:
      1. Pulls the total concentration from the concentration object
      2. Grabs the population by racial/ethnic group from the population object
      3. Merges the concentration and population data based on the ISRM ID
      4. Adds the population weighted mean exposure as a column of the geodataframe using `add_pwm_col`
2. `add_pwm_col`
   1. Inputs:
      * asdfasdf
   2. Outputs:
      * asdfasdf
   3. Methodology:
      * asdf
      * asdf
3. `get_pwm`
   1. Inputs:
      * asdfasdf
   2. Outputs:
      * asdfasdf
   3. Methodology:
      * asdf
      * asdf
4. `get_overall_disparity`
   1. Inputs:
      * asdfasdf
   2. Outputs:
      * asdfasdf
   3. Methodology:
      * asdf
      * asdf
5. `estimate_exposure_percentile`
   1. Inputs:
      * asdfasdf
   2. Outputs:
      * asdfasdf
   3. Methodology:
      * asdf
      * asdf
6. `run_exposure_calcs`
   1. Inputs:
      * asdfasdf
   2. Outputs:
      * asdfasdf
   3. Methodology:
      * asdf
      * asdf
7. `plot_percentile_exposure`
   1. Inputs:
      * asdfasdf
   2. Outputs:
      * asdfasdf
   3. Methodology:
      * asdf
      * asdf

#### `health_impact_calcs.py` 
Text goes here

*Function #1*
text

#### `tool_utils.py` 
Text goes here

*Function #1*
text
