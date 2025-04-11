# nwrfc-hydro-evolvingDDS

**Description** NWRFC Autocalibration tool of the National Weather Service River Forecast System Models (NWSRFS) using an evolving dynamically dimensioned search (edds). 

Primary Language:   R

Dependencies:  [nwrfc-hydro R package](https://github.com/geoffrey-walters/nwrfc-hydro-models)

Limitations:  Due to usage of fork implemenation of parallelism, does not work on a windows system

Acknowledgement:  Traditional dynamically dimensioned search (dds) alogorithm developed by building on original code by David Kneis [david.kneis@@tu-dresden.de], [dds.r github](https://github.com/dkneis/mcu/blob/master/R/dds.r)

## Autocalb Steps

### Required directory structure

See example directories and required file formats to utilize auto calibration code in `runs/1zone` and `runs/2zone`.

A general directory structure is as follows:  

- [LID]
  - flow_daily_[LID].csv **and/or**
  - flow_instantaneous_[LID].csv
  - forcing_por_[LID]-[zone #].csv
  - pars_default.csv
  - pars_limits.csv\
**<ins>Optional Files<ins>**
  - forcing_validation_cv_[fold #]_[LID]-[zone #].csv
  - upflow_[RR LID].csv

Notes: 
1. LID: Five alphanumeric characters (ex:  FSSO3)
2. zone #: Unique numeric character specifying zone. Autocalibration tool can accept multiple zones but requires at least one
3. fold #: Unique numeric character specifying fold number for cross validation (CV). Autocalibration tool can accept folds and is insenstitive to overlap between any folds
4. RR LID:  Five alphanumeric characters (ex:  WCHW1) representing a upstream reach for optimization of lagk parameters
5. Autocalibration can calibrate on either daily flow, instanteous flow, or both

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
* You can specify how many cores to use to perform the parallel optimization using the cores option.  The default option is eight cores. Option `--cores` Full will use all avalilable system cores minus 2   
* The NWRFC objfun naming convention can be interpreted as:  `<daily_metric>_<inst_metric>`
* The default is `nselognse_NULL`, which is read as the daily flow is using sum (nse, lognse) and `NULL` indicates instantaneous metrics are not included in the objective function score.  The table below lists base objective functions available.

#### Notes on `--obj_fun` Argument

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

1. Do not add the _obj part of the function name when using the `--objfun` argument:  CORRECT: `--objfun kge_NULL`, INCORRECT: `--objfun kge_NULL_obj`
2. Your custom function must accept the following inputs: results_daily, results_inst.
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

### ETC

When executing all of the script from the command line, you can use the `--help` argument to review the required arguments for a run.  Example:

```
./run-controller.R --help
```
The arguments which have a default option are listed.

