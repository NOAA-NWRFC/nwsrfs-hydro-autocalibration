# NWRFC Autocalibration Framework

**Description** This repository contains a version of the Northwest River Forecast Center (NWRFC) Autocalibration tool for parameterization of the National Weather Service River Forecast System Models (NWSRFS) using an evolving dynamically dimensioned search (edds).  NWSRFS  was originally developed in the late 1970s and remains a core component of the NWS Community Hydrologic Prediction System (CHPS). Suite of NWRFS models that can be calibrated simultaneously, with multiple zones, include:  SAC-SMA, SNOW17, Unit Hydrograph, LAGK, CHANLOSS, and CONS_USE.  See NWRFS documentation for more detail about each individual model [link](https://www.weather.gov/owp/oh_hrl_nwsrfs_users_manual_htm_xrfsdocpdf) .

**Language:**   R\
**Package Dependency:**  [nwrfc-hydro R package](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models)\
**Limitations:**  **(1)** Due to usage of fork implemenation of parallelism, does not work on a windows system. **(2)** Code has been tested only with a 6-hour time step. Use with other time steps may require additional configuration or validation.\
**Acknowledgement:**  Traditional dynamically dimensioned search (dds) alogorithm developed by building on original code by David Kneis [david.kneis@@tu-dresden.de], [dds.r github](https://github.com/dkneis/mcu/blob/master/R/dds.r)

## Autocalb steps

### Required directory structure

See example directories and required file formats to utilize auto calibration code in `runs/1zone` and `runs/2zone`.

A general directory structure is as follows:  

- [LID]
  - flow_daily_[LID].csv *[daily average flow observations]* **and/or**
  - flow_instantaneous_[LID].csv *[instanteous flow observations]*
  - forcing_por_[LID]-[zone #].csv *[model forcing data for each zone precipitation, MAP, MAT, and PTPS, for the period of record]*
  - pars_default.csv *[default parameter files (values of -99 indicate parameter to be optimized)]*
  - pars_limits.csv *[upper and lower limits for parameters to be optimized]*\
**<ins>Optional Files<ins>**
  - forcing_validation_cv_[fold #]_[LID]-[zone #].csv *[model forcing data for each zone,MAP, MAT, and PTPS, for a cross validation period]*
  - upflow_[RR LID].csv *[upstream flows for a routing reach/LAGK model]*

Notes: 
1. LID: Five alphanumeric characters (ex:  FSSO3)
2. zone #: Unique numeric character specifying zone. Autocalibration tool can accept multiple zones but requires at least one
3. fold #: Unique numeric character specifying fold number for cross validation (CV). Autocalibration tool can accept folds and is insenstitive to overlap between any folds
4. RR LID:  Five alphanumeric characters (ex:  WCHW1) representing a upstream reach for optimization of lagk parameters\

### run-controller.R 
run-controller.R script is run to create a optimized parameter file (pars_optimal.csv).  Example:
```
./run-controller.R --dir runs/1zone --objfun kge_NULL --basin FSSO3
```
* Only one basin can be calibrated at a time.
* Multiple runs of run-controller.R can be executed referencing the same directory with preprocess files. run-controller.R automatically gives each run a unique directory name (results_por_01, results_por_02 , results_por_03, etc.) unless specfied otherwise.
* A period of record (POR) run is the default.
* If you want to perform a CV run use the `--cvfold` option, enter the fold number you would like to run.
* Use the `--lite` option if you want to make a run with half the amount of iterations of a standard run. This is a good approach for initial configuration testing.
* You can use the `--overwrite` option to delete and rewrite the last results directory created for a POR or a particular CV run. This is helpful when doing multiple runs for testing and want to prevent a buildup of results directories.
* You can specify how many cores to use to perform the parallel optimization using the cores option.  The default option is eight cores. Option `--cores Full` will use all avalilable system cores minus 2
* Calibation can be conducted using daily observed flow, instanteous observed flow, or both   
* The NWRFC objfun naming convention can be interpreted as:  `<daily_metric>_<inst_metric>`
* The default is objfun  `nselognse_NULL`, which is read as the daily flow is using sum (nse, lognse) and `NULL` indicates instantaneous flow metrics are not included in the objective function score. 

#### Additional notes on `--obj_fun` Argument

Any calibrator interested in making their own objective function should read the comments at the beginning of the [obj_fun.R](https://github.com/geoffrey-walters/nwrfc-hydro-evolvingDDS/blob/main/obj_fun.R) script, which covers NWRFC naming conventions of the objective functions and strategies.   Below is a table of available objective functions as part of this repository:

Available Objective Functions ||
----------|-------------|
nselognse_NULL | lognse_r2
nsenpbias_NULL | kgelf_kge
kge_NULL | lognse_npbias95th
NULL_nselognse | kge_npbias2000q
NULL_kge | npbias050607m_kge
lognse_nse | nse.25wlognse.75w_npbias99th
lognse_kge | lognse.4W_nse1112010203m.6W 

1. Do not add the _obj part of the function name when using the `--objfun` argument:
   - CORRECT: `--objfun kge_NULL`
   - INCORRECT: `--objfun kge_NULL_obj`
2. Your custom function must accept both results_daily, results_inst as input.
3. If you attempt to use an objective function which relies on instantaneous data and no sub daily observations are available, a code error will be thrown.
4. If there is a code issue with the custom objective function, a descriptive code error will be thrown which starts with `“Objective Function had the following error, exiting:”`
 
### postprocess.R

Once run-controller.R has been ran and the pars_optimal.csv file has been created in a run directory, postprocess.R can be used to run create simulation timeseries csv files and other supporting tables.  Example:

```
./postprocess.R --dir runs/2zone --basins WGCM8 SAKW1
```

Notes:
- Multiple basins can be processed at once.
- Unless specified otherwise with the `--basins` argument, all completed auto calibration runs in the specified `--dir` argument will be processed
- Unless specified otherwise with the `--run` argument, all completed auto calibration runs in each basin directory will be processed

### cv-plots.R

Once a run directory contains both postprocessed period of record (POR) and cross validation (CV) runs, the cross validation analysis can be executed.  Example:

```
./cv-plots.R --dir runs/2zone --basins SFLN2
```

Notes:
- Multiple basins can be processed at once.
- Unless specified otherwise with the `--basins` argument, all completed auto calibration runs in the specified `--dir` argument will be processed
- If multiple versions POR or the same CV autocalibrations exists in a run directirory, the code will select the version with the highest KGE score.
- Plots compare metrics from the CV validation period and a stationary bootstrapping exercise using the POR run
  - During bootstapping, for each sample, x number of water year simulation and observations are random selected from the POR record run and metric are calculated.  x is equal to the average length of the CV periods.
  - Bootstrappig consists of 8,000 iterations
  
### ETC

#### Help

When executing all of the script from the command line, you can use the `--help` argument to review the required arguments for a run.  Example:

```
./run-controller.R --help
```
The arguments which have a default options are listed.

#### Example basin for calibration
Five example basins have been set up and have all the required file to run through the autocalibration tools.  Basin files can be found in the `runs/` directory, and are as follows:

FSSO3 (1zone):  Nehalem at Foss, rain dominated coastal basin in Oregon.  CAMELS basin
SAKW1 (2zone):  Sauk near Sauk, example of LAGK optization, with upstream routing reach from WCHW1.  CAMELS basin
SFLN2 (2zone):  Salmon Falls near San Jacinto, example of CONS_USE and CHANLOSS optimization 
WCHW1 (2zone):  Sauk above White Chuck, upstream of SAKW1.  CAMELS basin  
WGCM8 (2zone):  MF Flathead near West Glacier, snow dominated basin in Montana, CAMELS basin

#### CHANLOSS

- There can be mutiple CHANLOSS model associated with basin
- Start and end time can overlap between CHANLOSS model. The average value of the optimized values between models will be used for the overlap periods
- Start and End time can cross from one year to the next, for example a start of 11 and end of 2 would span from November to Febuary
- The `cl_type` parameter can be either 1 or 2.  A basin using the CHANLOSS model can only have one `cl_type`
  - 1:  VARP CHANLOSS adjustment, parameter value is equivalent to a multiplication factor applied to the sim
  - 2:  VARC CHANLOSS adjustment, parameter value is subtracted from the sim

### Forcings

There are three required forcings for each zone:  Mean areal precipitation (MAP), meean areal temperature (MAT), and percent precipitation as snow (PTPS).

Be aware that forcing data is considered end of timestep, where as simulated flow is considered beginning of time step.  When using tool, depending on how forcing data are constructed, timeseries may need to be shifted back one time step.  

### AdjustQ

For calibration of the LAGK model, the NWRFC derives the upstream routing reach timeseries using the [AdjustQ procedure](https://publicwiki.deltares.nl/display/FEWSDOC/AdjustQ).  The [nwrfc-hydro R package](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models) has [equivalent code](https://github.com/NOAA-NWRFC/nwsrfs-hydro-models/blob/main/py-rfchydromodels/utilities/adjustq.py) to the deltaires transformation written in python.  

## Credits and references

Please use the folowing journal article for referencing this work:

Walters, G., Bracken, C., et al., "A comprehensive calibration framework for the Northwest River Forecast Center." Unpublished manuscript, Submitted 2025, JAWA Journal of the American Water Resources Association

If you wish to use or adapt the code in this repository, please make sure that your new repository credits this one as the original source of the code. 

## Legal disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

[NOAA GitHub Policy](https://github.com/NOAAGov/Information)
 \
 \
 \
<img src="https://www.weather.gov/bundles/templating/images/header/header.png" alt="NWS-NOAA Banner">

[National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [National Weather Service](https://www.weather.gov/) | [Northwest River Forecast Center](https://www.nwrfc.noaa.gov/rfc/)


