# nwrfc-hydro-evolvingDDS
NWRFC Autocalibration tool using an evolving dynamically dimensioned search

## Autocalb Steps

### run-controller.R 
run-controller.R script is run to optimize the parameters:
Example: ./run-controller.R --dir runs/1zone --objfun kge_NULL --basin FSSO3

a. Only one basin can be calibrated at a time.
b. Multiple runs of run-controller.R can be executed referencing the same directory with preprocess files. run-controller.R automatically gives each run a unique directory name (results_por_01, results_por_02 , results_por_03, etc.)
c. The objfun naming convention can be interpreted as: 
  _<daily_metric>_<inst_metric> _
  The default is nselognse_NULL, which is read as the daily flow is using sum (nse, lognse) and NULL indicates instantaneous metrics are not included in the objective function score.  The table below lists base objective functions available.

objfun
nselognse_NULL
lognse_r2
nsenpbias_NULL
kgelf_kge
kge_NULL
lognse_npbias95th
NULL_nselognse 
kge_npbias2000q
NULL_kge
npbias050607m_kge
lognse_nse
nse.25wlognse.75w_npbias99th
lognse_kge
lognse.4W_nse1112010203m.6W 

Any calibrator interested in making their own objective function should read the comments at the beginning of the obj_fun.R script, which covers naming conventions of the objective functions and strategies. It is very important to follow the naming convention, otherwise unexpected results may occur. Additional notes about objfun argument:

i.  Do not add the _obj part of the function name when using the objfun argument:  CORRECT: --objfun kge_NULL, INCORRECT: --objfun kge_NULL_obj
ii.  Your custom function must accept the following inputs: results_daily, results_inst.
iii.  If you attempt to use an objective function which relies on instantaneous data and no sub daily observations are available, a code error will be thrown.  
iv.  When writing your own objective function, use the section in obj_fun.R script titled “custom objective function space.”
v.  When writing your own objective function, the easiest strategy is to copy an existing base objective function and alter as needed to fit your needs.
vi.  If there is a code issue with the custom objective function, a descriptive code error will be thrown which starts with “Objective Function had the following error, exiting:”

d.A POR run is the default.
e.If you want to perform a CV run use the cvfold option, enter the fold number you would like to run.
f. Use the lite option if you want to make a run with half the amount of iterations of a standard run. This is a good approach for initial configuration testing, but should not be used for any runs being shared for evaluation.
g. You can use the overwrite option to delete and rewrite the last results directory created for a POR or a particular CV# run. This is helpful when doing multiple runs for testing and want to prevent a buildup of results directories.
h. You can specify how many cores to use to perform the parallel optimization using the cores option. In most cases the default option is adequate. Option --cores Full will use all avalilable system cores minus 2    

### postprocess.R

Ex: ./postprocess.R --dir runs/2zone --basins WGCM8 SAKW1

Multiple basins can be processed at once.

### ETC

Note: When executing all of the script from the command line, you can use the help argument to review the required arguments for a run.  
    _Ex:  ./run-controller.R --help_
The arguments which have a default option are listed when using the help option.
