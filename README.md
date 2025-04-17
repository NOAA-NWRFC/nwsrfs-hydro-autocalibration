# NWRFC Autocalibration Framework

**Description**  
This repository contains a version of the Northwest River Forecast Center (NWRFC) autocalibration tool for parameterizing the National Weather Service River Forecast System (NWSRFS) models using an evolving dynamically dimensioned search (EDDS). NWSRFS, originally developed in the late 1970s, remains a core component of the NWS Community Hydrologic Prediction System (CHPS).  This framework supports simultaneous calibration of a suite of NWSRFS models across multiple zones, including: SAC-SMA, SNOW17, Unit Hydrograph, LAGK, CHANLOSS, and CONS_USE.  See the [NWSRFS documentation](https://www.weather.gov/owp/oh_hrl_nwsrfs_users_manual_htm_xrfsdocpdf) for more detail on each individual model.

**Language:** R  
**Package Dependency:** [nwrfc-hydro R package](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models)  
**Limitations:**
1. Due to its use of fork-based parallelism, the tool is not compatible with Windows systems.
2. The code has been tested only with a 6-hour timestep. Use with other timesteps may require additional configuration and validation.

**Acknowledgment:**  
The traditional dynamically dimensioned search (DDS) algorithm builds on original code by David Kneis ([david.kneis@tu-dresden.de](mailto:david.kneis@tu-dresden.de)). See: [dds.r GitHub](https://github.com/dkneis/mcu/blob/master/R/dds.r)

## Autocalibration Steps

### Required Directory Structure

See example directories and required file formats in `runs/1zone` and `runs/2zone`.

```
[LID]/
├── flow_daily_[LID].csv             # Daily average flow observations (optional)
├── flow_instantaneous_[LID].csv     # Instantaneous flow observations (optional)
├── forcing_por_[LID]-[zone #].csv   # Forcing data for each zone (MAP, MAT, PTPS)
├── pars_default.csv                 # Default parameter file (-99 indicates optimization)
├── pars_limits.csv                  # Upper/lower limits for parameters
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

### `run-controller.R`

Creates the optimized parameter file (`pars_optimal.csv`).

**Example:**
```bash
./run-controller.R --dir runs/1zone --objfun kge_NULL --basin FSSO3
```

**Notes:**
- Only one basin can be calibrated at a time.
- Multiple runs can use the same directory; results are placed in `results_por_01`, `results_por_02`, etc.
- CV runs: use `--cvfold [#]`.
- Light run (fewer iterations): use `--lite`.
- Overwrite last results directory: use `--overwrite`.
- Control number of cores: `--cores [#]` or `--cores Full` (uses all available minus 2).
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

For custom objective functions, refer to comments in [obj_fun.R](https://github.com/geoffrey-walters/nwrfc-hydro-evolvingDDS/blob/main/obj_fun.R).  

**Notes:**
1. Do not include `_obj` portion of function name in argument.
2. Functions must accept `results_daily`, `results_inst` as inputs.
3. Selection of objective function should consider availability of daily and instanteous flow observations.
4. Errors in custom functions produce descriptive messages starting with  
   `"Objective Function had the following error, exiting:"`
 
### `postprocess.R`

Processes output to create simulation time series and supporting tables.

```bash
./postprocess.R --dir runs/2zone --basins SFLN2 --run results_cv_3_01
```

**Notes:**
- Processes all completed runs in `--dir` path by default unless `--basins` or `--run` is specified.

### `cv-plots.R`

Analyzes and visualizes results from POR and CV runs.

```bash
./cv-plots.R --dir runs/2zone --basins WGCM8 SAKW1
```
**Notes:**
- When multiple POR or CV runs exist, selects the run with the highest KGE score .
- Plots compare CV metrics vs. stationary bootstrap from POR.
- Bootstrapping draws `x`-year samples from POR (where `x` = average CV fold length).
  - 8,000 bootstrap iterations performed.
  
## Additional Info

#### Help

Use `--help` to view argument options:
```bash
./run-controller.R --help
```

### Example Basins

| Basin  |  Name |Zones | Description |
|--------|-------|-------|-------------|
| FSSO3  | Nehalem at Foss, OR|1     | Rain-dominated (CAMELS) |
| SAKW1  | Sauk nr Sauk, WA|2     |  Rain/Snow-dominated, LAGK example (CAMELS)|
| SFLN2  | Salmon Falls nr San Jacinto, NV|2     | Arrid basin, CONS_USE and CHANLOSS example |
| WCHW1  | Sauk ab White Chuck, WA|2     | Rain/Snow-dominated, routing reach to SAKW1 (CAMELS) |
| WGCM8  | MF Flathead nr W Glacier, MT|2     | Snow-dominated (CAMELS) |

*supporting files are stored in the `runs/` directory

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

## Credits and references

Please cite the following work when using this tool:

Walters, G., Bracken, C., et al., "A comprehensive calibration framework for the Northwest River Forecast Center." Unpublished manuscript, Submitted 2025, JAWA Journal of the American Water Resources Association

If adapting this code, please credit this repository as the original source. 

## Legal disclaimer

This is a scientific product and does not represent official communication from NOAA or the U.S. Department of Commerce. All code is provided "as is."

See full disclaimer: [NOAA GitHub Policy](https://github.com/NOAAGov/Information)
 \
 \
 \
<img src="https://www.weather.gov/bundles/templating/images/header/header.png" alt="NWS-NOAA Banner">

[National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [National Weather Service](https://www.weather.gov/) | [Northwest River Forecast Center](https://www.nwrfc.noaa.gov/rfc/)


