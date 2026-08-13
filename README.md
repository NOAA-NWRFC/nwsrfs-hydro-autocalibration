# NWRFC Autocalibration Framework

## Description
This repository contains a version of the Northwest River Forecast Center (NWRFC) autocalibration tool for parameterizing the National Weather Service River Forecast System (NWSRFS) models using an evolving dynamically dimensioned search (EDDS). NWSRFS, originally developed in the late 1970s, remains a core component of the NWS Community Hydrologic Prediction System (CHPS).  This framework supports simultaneous calibration of a suite of NWSRFS models across multiple zones, including: SAC-SMA, SNOW-17, Unit Hydrograph, LAG-K, CHANLOSS, and CONS_USE.  See the [NWSRFS documentation](https://www.weather.gov/owp/oh_hrl_nwsrfs_users_manual_htm_xrfsdocpdf) for more detail on each individual model.

**Language:** R  
**Package Dependency:** [nwsrfs-hydro-models R package](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models/nwsrfs_r)  


## Prerequisites

1. Install [R](http://r-project.org) (Version 4.5 or later). 

2. Install [pixi](https://pixi.prefix.dev/latest/)

**NOTES:**

The code has been tested only with a 6-hour timestep. Use with other timesteps may require additional configuration and validation.

## Example Calibrations
There are five basin directories included in this repo that serve as examples which utilize all the features of the auto calibration tool.

**Example Basins**
| NWSLI ID  |  Name | USGS # | Zones | Description |
|--------|-------|-------|-------|-------------|
| FSSO3  | Nehalem at Foss, OR|14301000 |1     | Rain-dominated (CAMELS) |
| SAKW1  | Sauk nr Sauk, WA|12189500 |2     |  Rain/Snow-dominated, LAG-K example (CAMELS)|
| SFLN2  | Salmon Falls nr San Jacinto, NV|13105000 |2     | Arrid basin, CONS_USE and CHANLOSS example |
| WCHW1  | Sauk ab White Chuck, WA|12186000|2     | Rain/Snow-dominated, routing reach to SAKW1 (CAMELS) |
| WGCM8  | MF Flathead nr W Glacier, MT|12358500 |2     | Snow-dominated (CAMELS) |

*supporting files are stored in the `runs/` directory in the repo

**Archived sample Data**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18829935.svg)](https://doi.org/10.5281/zenodo.18829935)

Input data for the example basins is archived on Zenodo for reproducibility. 

## Set up the auto calibration scripts and the environmennt

```bash
    # get repo with auto-calibration scripts
    # sample data is in the repo already, but is also archived on zenodo
    # https://doi.org/10.5281/zenodo.18829935
    git clone https://github.com/NOAA-NWRFC/nwsrfs-hydro-autocalibration.git
    cd nwsrfs-hydro-autocalibration
    # set up environment
    pixi install
    pixi run setup
    # now you are ready to run the workflow below
```

### Example Workflow 
 We recommend that you complete at least 4 cross validation runs in addition to the full period of record run to evaluate the calibration for any potential issues.
 
    # period of record run
    ./run-controller.R --dir runs/2zone --objfun lognse_kge --basin WGCM8

    # cross validation
    ./run-controller.R --dir runs/2zone --objfun lognse_kge --basin WGCM8 --cvfold 1
    ./run-controller.R --dir runs/2zone --objfun lognse_kge --basin WGCM8 --cvfold 2
    ./run-controller.R --dir runs/2zone --objfun lognse_kge --basin WGCM8 --cvfold 3
    ./run-controller.R --dir runs/2zone --objfun lognse_kge --basin WGCM8 --cvfold 4

    # postprocessing
    ./postprocess.R --dir runs/2zone --basins WGCM8

    # cross-validation
    ./cv-plots.R --dir runs/2zone --basins WGCM8

Any of the example basins could be swapped into this same workflow.

## Required Directory Structure

Refer to the example basins in the `runs/` directory for the expected directory structure and file formats.

```
[LID]/
  ├── flow_daily_[LID].csv             # Daily average flow observations (optional)
  ├── flow_instantaneous_[LID].csv     # Instantaneous flow observations (optional)
  ├── forcing_por_[LID]-[zone #].csv   # Forcing data for each zone (MAP, MAT, PTPS)
  ├── pars_default.csv                 # Default parameter file (-99 indicates the parameter will be optimized)
  ├── pars_limits.csv                  # Upper/lower limits for parameters that are optimized
  ├── [optional files...]
```

**Optional Files:**
- `forcing_validation_cv_[fold #]_[LID]-[zone #].csv`: Forcing data for cross-validation folds. Note that the data for each cross validation fold must be created manually by subsetting your data, but any number of folds is possible. as long as they split the data into even groups. 
- `upflow_[RR LID].csv`: Upstream flow data for routing reach (LAG-K model). A reach may have more than one routed upstream flow input. 

**Notes:**
- `LID`: 5-character basin ID (e.g., `FSSO3`). Note that this is an arbitrary basin identification code, you may swap in any unique 5 character alphanumeric identifier. 
- `zone #`: Numeric zone ID (at least one required)
- `fold #`: Numeric ID for cross-validation fold, starting at 1. 
- `RR LID`: Upstream reach LID (e.g., `WCHW1` for LAG-K optimization)
- Need at least one daily or instanteous flow file for autocalibration 

## Autocalibration Steps

### 1. `run-controller.R`

The `run-controller.R` script is run to create an optimized parameter file (`pars_optimal.csv`).  

    usage: run-controller.R [--] [--help] [--por] [--overwrite] [--lite]
          [--dir DIR] [--basin BASIN] [--objfun OBJFUN] [--optimizer
          OPTIMIZER] [--cvfold CVFOLD] [--num_cores NUM_CORES]

    Auto-calibration run controller

    flags:
      -h, --help        show this help message and exit
      -p, --por         Do a period of record run [default]
      -ov, --overwrite  Don't create new results dir, overwrite the first exising one
      -l, --lite        Testing run with 1/2 the total optimizer iteration

    optional arguments:
      -d, --dir         Input directory path
      -b, --basin       Basin name
      -o, --objfun      Objective function name [default: nselognse_NULL]
      --optimizer       Optimzer to use {edds [default],pso,dds} [default:
                        edds]
      -c, --cvfold      CV fold to run (integer 1-4) [default: none]
      -n, --num_cores   Number of cores to allocate for run, FULL uses all
                        available cores -2 [default: FULL]


**Example:**
```bash
./run-controller.R --dir runs/1zone --objfun kge_NULL --basin FSSO3
```

**Notes:**

- The script can only calibrate one basin at a time, although multiple can be done manually by utilizing only a portion of avaiable cores (`--num_cores #`) and running multiple scripts simultaneously (i.e. poor man's parallelization).
- Multiple runs can use the same input directory; results are placed in squentially numbered output directories, `results_por_01`, `results_por_02`, etc.
- For cross validation runs: use `--cvfold [#]`. The cross validation fold number must have a corresponding forcing file in the basin directory.  See [Required Directory Structure](#required-directory-structure) section for more information on the cross validation forcing file.
- Light run (half the number of optimizer iterations): use `--lite`. Note that this will take less real time to run but may result in a lower quality parameter set. 
- To overwrite last results directory, i.e. don't increment the output directory: use `--overwrite` (ignored if no results exist).
- To control number of utilized CPU cores: `--cores [#]` or `--cores FULL` (uses all available minus 2).
- The number of iterations is set to 5000 (or 2500 for a lite run) and is intentionally not user editiable based on extensive testing of how many iterations will produce stable parameter sets. Note that it is still important to manage equally viable parameter sets (i.e. equifinal solutions) by setting appropriate ranges for optimized parameter sets. This is particularly important for multi zone basins. 
- Increasing the number of cores used will not speed up the calibration but may make the calibration converge faster or come to a better overall solution. 
- The optimizer supports calibration with daily average, instantaneous, or both flow types.

#### Objective Function Argument (`--objfun`)

**Objective Function Naming Convention:**  
Format: `<daily_metric>_<instant_metric>`  
Default: `nselognse_NULL`

**Available Objective Functions:**
| <daily_metric>_<instant_metric> | Additional Options|
|----------------------------|---------------------------------|
| nselognse_NULL             | lognse_r2                       |
| nsenpbias_NULL             | kgelf_kge                       |
| kge_NULL                   | lognse_npbias95th              |
| NULL_nselognse             | kge_npbias2000q                |
| NULL_kge                   | npbias050607m_kge              |
| lognse_nse                 | nse.25wlognse.75w_npbias99th   |
| lognse_kge                 | lognse.4W_nse1112010203m.6W    |

To create a custom objective function, edit the [obj_fun.R](https://github.com/geoffrey-walters/nwrfc-hydro-evolvingDDS/blob/main/obj_fun.R) file and add your own function.

**Notes:**
1. Do not include `_obj` portion of function name in argument.
2. Functions must accept `results_daily`, `results_inst` as inputs.
3. Selection of objective function should consider availability of daily and instanteous flow observations.
4. Errors in custom functions should produce descriptive messages when running `run-controller.R` starting with  
   `"Objective Function had the following error, exiting:"`
5.   The [obj_fun.R](https://github.com/geoffrey-walters/nwrfc-hydro-evolvingDDS/blob/main/obj_fun.R) file includes additional guidance.
 
### 2. `postprocess.R`

Once `run-controller.R` has been ran and the `pars_optimal.csv` file has been created in a run directory, `postprocess.R` can be run to create simulation timeseries csv files and other supporting tables.  Data files are output into `<dir>/<basin>/<run>` and plots are output into `<dir>/<basin>/<run>/plots`.

    usage: postprocess.R  [--] [--help] [--dir DIR] [--basins BASINS] 
          [--run RUN]

    Auto-calibration postprocessor

    flags:
      -h, --help        show this help message and exit

    optional arguments:
      -d, --dir         Input directory path containing basin directories
      -b, --basins      Basins to run
      -r, --run         Output directories to postprocess

**Example:**
```bash
./postprocess.R --dir runs/2zone --basins SFLN2 --run results_cv_3_01
```

**Notes:**
- The script processes all completed runs for all basins in the `--dir` path by default unless `--basins` or `--run` is specified.

### 3. `cv-plots.R`

`cv-plots.R`, generates plots comparing CV metrics with results from a stationary bootstrap of the POR run. Plots are output into `<dir>/cv_plots`.

    usage: cv-plot.R  [--] [--help] [--cleanup] [--dir DIR]
          [--basins BASINS] 

    Creates CV Plots

    flags:
      -h, --help        show this help message and exit
      -c, --cleanup     Option to delete results directories which are not 
                        used for CV analysis
                        
    optional arguments:
      -d, --dir         Input directory path containing basin directories
      -b, --basins      Basins to run


**Example:**
```bash
./cv-plots.R --dir runs/2zone --basins WGCM8 SAKW1
```

**Notes:**
- When multiple POR or CV runs exist, the runs with the highest KGE score are used to develop the CV plots.
- If the `--basins` argument is omitted then it runs all available basins. 
- Bootstrapping draws `x`-year samples from POR (where `x` = average CV fold length).
  - 8,000 bootstrap iterations performed.
  
## Additional Info

#### Help

Use `--help` to view argument options:
```bash
./run-controller.R --help
./postprocess.R --help
./cv-plots.R --help
```

#### CHANLOSS Model

- Multiple CHANLOSS models can exist per basin.
- Overlapping time ranges are averaged.
- Start/end months can span across years (e.g., Nov–Feb = 11–2).
- `cl_type` must be 1 or 2:
  - 1: VARP adjustment (multiplier on simulation)
  - 2: VARC adjustment (subtracted from simulation)

### Forcings

Each zone requires:
- MAP: Mean Areal Precipitation
- MAT: Mean Areal Temperature
- PTPS: % Precipitation as Snow

**Note:** Forcing data is **end of timestep**; simulated flow is **start of timestep**. Observational time series may need shifting forward/back one timestep to comply with this requirement. 

### AdjustQ

For LAG-K calibration, upstream flows are derived using [AdjustQ](https://publicwiki.deltares.nl/display/FEWSDOC/AdjustQ).  

See [nwsrfs-hydro-models](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models) for equivalent [R](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models/nwsrfs_r) or [Python](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models/nwsrfs_py) utility as part of their respective packages.

### Forcing Climatological Corrections

In the NWRFC autocalibration scheme, mid-month climatological adjustment factors are optimized independently for each forcing variable—precipitation, temperature, precipitation typing, and potential evaporation.  To disable climatological corrections for MAT, MAP, PTPS, or PET:

1. Remove any related lines in `pars_limits.csv` (e.g., `mat_shift_[LID]`)
2.  Set the following in pars_default.csv for each variable:
```
[forcing]_scale = 1
[forcing]_p_redist = 0
[forcing]_std = 10
[forcing]_shift = 0
```
     
## Credits and References

Please cite the following work when using this tool:

Walters, G., C. Bracken, B. Gillies, et al. 2026. “A Comprehensive Calibration Framework for the Northwest River Forecast Center.” *JAWRA Journal of the American Water Resources Association* 62, no. 2: e70112. [https://doi.org/10.1111/1752-1688.70112](https://doi.org/10.1111/1752-1688.70112)


If adapting this code, please credit this repository as the original source. 

## Acknowledgment
The traditional dynamically dimensioned search (DDS) algorithm builds on original code by David Kneis ([david.kneis@tu-dresden.de](mailto:david.kneis@tu-dresden.de)). See: [dds.r GitHub](https://github.com/dkneis/mcu/blob/master/R/dds.r)

## Legal Disclaimer

This is a scientific product and does not represent official communication from NOAA or the U.S. Department of Commerce. All code is provided "as is."

See full disclaimer: [NOAA GitHub Policy](https://github.com/NOAAGov/Information)
 \
 \
 \
<img src="https://www.weather.gov/bundles/templating/images/header/header.png" alt="NWS-NOAA Banner">

[National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [National Weather Service](https://www.weather.gov/) | [Northwest River Forecast Center](https://www.nwrfc.noaa.gov/rfc/)


