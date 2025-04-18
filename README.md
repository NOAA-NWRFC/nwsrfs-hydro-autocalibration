# NWRFC Autocalibration Framework

## Description
This repository contains a version of the Northwest River Forecast Center (NWRFC) autocalibration tool for parameterizing the National Weather Service River Forecast System (NWSRFS) models using an evolving dynamically dimensioned search (EDDS). NWSRFS, originally developed in the late 1970s, remains a core component of the NWS Community Hydrologic Prediction System (CHPS).  This framework supports simultaneous calibration of a suite of NWSRFS models across multiple zones, including: SAC-SMA, SNOW17, Unit Hydrograph, LAGK, CHANLOSS, and CONS_USE.  See the [NWSRFS documentation](https://www.weather.gov/owp/oh_hrl_nwsrfs_users_manual_htm_xrfsdocpdf) for more detail on each individual model.

**Language:** R  
**Package Dependency:** [nwrfc-hydro R package](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models)  


## Prerequisites

1. Install [R](http://r-project.org). 

2. Install these R packages: 

        install.packages(c('xfun','import','devtools'))

3. Install the `rfchydromodels` R package which requires a Fortran complier. This package has been tested with [gfortran](https://gcc.gnu.org/wiki/GFortran). See [here](https://cran.r-project.org/bin/macosx/tools/) for an easy option on MacOS.
    
From R:

        devtools::install_github('NOAA-NWRFC/nwsrfs-hydro-models',subdir='rfchydromodels')

or from the command line:

        git clone https://github.com/NOAA-NWRFC/nwsrfs-hydro-models.git
        cd nwsrfs-hydro-models
        R CMD INSTALL rfchydromodels

4. The autocalibration scripts will try to install a number of R packages when run. If this fails you may need to install the packages manually. 

**NOTES:**

1. Due to its use of fork-based parallelism, the tool is not compatible with Windows systems.

2. The code has been tested only with a 6-hour timestep. Use with other timesteps may require additional configuration and validation.

## Example Calibrations
There are five basin directories included in this repo that serve as examples which utilize all the features of the auto calibration tool.

**Example Basins**
| NWSLI ID  |  Name | USGS # | Zones | Description |
|--------|-------|-------|-------|-------------|
| FSSO3  | Nehalem at Foss, OR|14301000 |1     | Rain-dominated (CAMELS) |
| SAKW1  | Sauk nr Sauk, WA|12189500 |2     |  Rain/Snow-dominated, LAGK example (CAMELS)|
| SFLN2  | Salmon Falls nr San Jacinto, NV|13105000 |2     | Arrid basin, CONS_USE and CHANLOSS example |
| WCHW1  | Sauk ab White Chuck, WA|12186000|2     | Rain/Snow-dominated, routing reach to SAKW1 (CAMELS) |
| WGCM8  | MF Flathead nr W Glacier, MT|12358500 |2     | Snow-dominated (CAMELS) |

*supporting files are stored in the `runs/` directory

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
    ./cv_plots.R --dir runs/2zone --basins WGCM8

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
- `forcing_validation_cv_[fold #]_[LID]-[zone #].csv`: Forcing data for cross-validation folds.
- `upflow_[RR LID].csv`: Upstream flow data for routing reach (LAGK model).

**Notes:**
- `LID`: 5-character basin ID (e.g., `FSSO3`)
- `zone #`: Numeric zone ID (at least one required)
- `fold #`: Numeric ID for cross-validation fold
- `RR LID`: Upstream reach LID (e.g., `WCHW1` for LAGK optimization)
- Need at least one daily or instanteous flow file for autocalibration 

## Autocalibration Steps

### 1. `run-controller.R`

The `run-controller.R` script is run to create a optimized parameter file (`pars_optimal.csv`).  

    usage: run-controller.R [--] [--help] [--por] [--overwrite] [--lite]
          [--dir DIR] [--basin BASIN] [--objfun OBJFUN] [--optimizer
          OPTIMIZER] [--cvfold CVFOLD] [--num_cores NUM_CORES]

    Auto-calibration run controller

    flags:
      -h, --help        show this help message and exit
      -p, --por         Do a period of record run [default]
      -ov, --overwrite  Don't create new results dir, overwrite
      -l, --lite        Testing run with 1/2 the total optimizer iteration

    optional arguments:
      -d, --dir         Input directory path
      -b, --basin       Basin name
      -o, --objfun      Objective function name [default: nselognse_NULL]
      --optimizer       Optimzer to use {edds [default],pso,dds} [default:
                        edds]
      -c, --cvfold      CV fold to run (integer 1-4) [default: none]
      -n, --num_cores   Number of cores to allocate for run, FULL uses all
                        availavble cores -2 [default: FULL]


**Example:**
```bash
./run-controller.R --dir runs/1zone --objfun kge_NULL --basin FSSO3
```

**Notes:**

- The script can only calibrate one basin a time, although multiple can be done manually by utilizing only a portion of avaiable cores.
- Multiple runs can use the same directory; results are placed in squentially numbered output directories, `results_por_01`, `results_por_02`, etc.
- For cross validation runs: use `--cvfold [#]`. The cross validation fold number must have a correspond forcing file in the basin directory.  See [Required Directory Structure](#required-directory-structure) section for more information on the cross validation forcing file.
- Light run (fewer optimizer iterations): use `--lite`.
- Overwrite last results directory, i.e. cont increment the output directory: use `--overwrite`.
- Control number of cores: `--cores [#]` or `--cores Full` (uses all available minus 2).
- The number of iterations is set to 5000 (or 2500 for a lite run) and not user editiable intentionally based on extensive testing of what will produce stable parameter sets. 
- Increasing the number of cores used will not speed up the calibration but may make the calibration converge faster or come to a better overall solution. 
- Supports calibration with daily, instantaneous, or both flow types.

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

To create a custom objective function, edit the [obj_fun.R](https://github.com/geoffrey-walters/nwrfc-hydro-evolvingDDS/blob/main/obj_fun.R) file.

**Notes:**
1. Do not include `_obj` portion of function name in argument.
2. Functions must accept `results_daily`, `results_inst` as inputs.
3. Selection of objective function should consider availability of daily and instanteous flow observations.
4. Errors in custom functions produce descriptive messages when running `run-controller.R` starting with  
   `"Objective Function had the following error, exiting:"`
5.   The [obj_fun.R](https://github.com/geoffrey-walters/nwrfc-hydro-evolvingDDS/blob/main/obj_fun.R) file includes additional guidance.
 
### 2. `postprocess.R`

Once `run-controller.R` has been ran and the `pars_optimal.csv` file has been created in a run directory, `postprocess.R` can be run to create simulation timeseries csv files and other supporting tables.  Data files are outut into `<dir>/<basin>` and plots are output into `<dir>/<basin>/plots`.

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
- The script processes all completed runs in `--dir` path by default unless `--basins` or `--run` is specified.

### 3. `cv-plots.R`

`cv-plots.R`, generates plots comparing CV metrics with results from a stationary bootstrap of the POR run.

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

For LAGK calibration, upstream flows are derived using [AdjustQ](https://publicwiki.deltares.nl/display/FEWSDOC/AdjustQ).  

See [nwrfc-hydro R package](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models) for [equivalent Python code](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models/blob/main/py-rfchydromodels/utilities/adjustq.py).

## Credits and References

Please cite the following work when using this tool:

Walters, G., Bracken, C., et al., "A comprehensive calibration framework for the Northwest River Forecast Center." Unpublished manuscript, Submitted 2025, JAWRA Journal of the American Water Resources Association

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


