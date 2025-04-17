#!/usr/bin/env Rscript

# Written by Cameron Bracken and Geoffrey Walters (2025)
# Please see the LICENSE file for license information

# install.packages(c('xfun','import','devtools'))
xfun::pkg_load2('magrittr','dplyr', 'data.table', 'dtplyr','hydroGOF',
                'digest','lubridate','readr','tibble','ggthemes',
                'crayon','argparser','rtop')
xfun::pkg_attach2('ggplot2')
import::from(dplyr, filter, select, summarise, group_by, ungroup, mutate,
             arrange, lag,lead,full_join,right_join,bind_rows)
select <- dplyr::select
import::from(data.table, as.data.table, data.table, fread, merge.data.table, copy,
             melt.data.table, rbindlist, dcast.data.table, fwrite, setnames, nafill)
import::from(dtplyr, lazy_dt)
import::from(rtop,sceua)
import::from(hydroGOF,gof,NSE,rmse,mNSE,rSD,pbias,rPearson,KGE,KGElf,KGEnp)
import::from(digest, digest)
import::from(lubridate, ymd_hms)
import::from(readr, read_delim, read_csv, cols)
import::from(tibble,tibble,as_tibble,rownames_to_column)
import::from(ggthemes, colorblind_pal)
import::from(crayon, '%+%', green, bold, bgGreen, blue, inverse, red, bgCyan)
import::from(stringr, str_subset, str_detect)
import::from(argparser,arg_parser,add_argument,parse_args)
import::from(tidyr,fill)
import::from(parallel,detectCores,makeCluster,clusterSetRNGStream,
             clusterCall,nextRNGStream,stopCluster)
import::from(vctrs,vec_fill_missing)
##############################################################################
#!!To ensure correct parameter/forcings/upstream flow get mapped correctly to
#intended zone/process. The default or optimal parameter file, forcing list, and
#upstream flow list MUST be alphabetized by zone or basin names prior to call wrapper.r
#or rfchydromodels functions. 

import::from(rfchydromodels, sac_snow_uh, sac_snow_uh_lagk, lagk, sac_snow, sac_snow_states,
             uh, consuse, chanloss, fa_nwrfc)
source('wrappers.R')
source('obj_fun.R')

parser = arg_parser("Auto-calibration run controller", hide.opts = TRUE)

# by default ArgumentParser will add an help option
parser = add_argument(parser, "--dir", help = "Input directory path")
parser = add_argument(parser, "--basin", help = "Basin name")
parser = add_argument(parser, "--objfun", default='nselognse_NULL', help = "Objective function name")
parser = add_argument(parser, "--optimizer", default='edds', help = "Optimzer to use {edds [default],pso,dds}")
parser = add_argument(parser, "--cvfold", default=NA_integer_, help = "CV fold to run (integer 1-4)")
parser = add_argument(parser, "--num_cores", default="FULL", help = "Number of cores to allocate for run, FULL uses all availavble cores -2")
parser = add_argument(parser, "--por", flag=TRUE, help = "Do a period of record run [default]")
parser = add_argument(parser, "--overwrite", flag=TRUE, help = "Don't create new results dir, overwrite",short = '-ov')
parser = add_argument(parser, "--lite",flag=TRUE, help = "Testing run with 1/2 the total optimizer iteration")


args = parse_args(parser)

# standard objective functions are:
# <daily_metric>_<inst_metric>
# Add examples

#args = list(dir='runs/opt_run/1zone',lite=TRUE,num_cores=8,
#             basin='TLMO3',objfun='lognse_kge',cvfold=4,por=FALSE,overwrite=TRUE)
#args = list(dir='runs/opt_run/2zone',basin='SFLN2',lite=TRUE,num_cores=8,
#             objfun='nselognse_NULL',cvfold=NA,por=TRUE,overwrite=TRUE)
#args = list(dir='runs/opt_run/2zone',basin='AUBW1', lite=TRUE, num_cores='FULL',
#              objfun='nselognse_NULL',cvfold=NA,por=TRUE,overwrite=TRUE)
#args = list(dir='runs/opt_run/routing_reach',basin='NREW1', lite=TRUE, num_cores='FULL',
#            objfun='nselognse_NULL',cvfold=1,por=FALSE,overwrite=TRUE)

# all basins are sub dirs under here
output_dir = args$dir
basin = args$basin
fold = as.integer(args$cvfold)
por = args$por
obj_fun = args$objfun
optimizer = 'edds'
overwrite = args$overwrite
lite = args$lite
n_cores = args$num_cores
dt_hours = 6

#SET UP NUMBER OF CORES TO USE FOR OPTIMIZATION
#if a number is passed then convert from string to number
if(suppressWarnings(!is.na(as.numeric(n_cores)))) n_cores=as.numeric(n_cores)
#stop run if the number of cores argument is not understood or too many are allocated
#if FULL argument is passed then utilize all cores but 2
if(is.character(n_cores)&toupper(n_cores)!='FULL'){
  stop('num_cores argument not understood, exiting')
}else if((is.character(n_cores)&toupper(n_cores)=='FULL')|(n_cores>detectCores() - 2)){
  n_cores=detectCores() - 2
}

basin_dir = file.path(output_dir,basin)
#CHECK IF DIR ARGUMENT EXISTS
if(!file.exists(basin_dir)){
  stop(basin_dir,' does not exist, exiting')
}

#CHECK IF OBJECTIVE FUNCTION EXISTS
if(!exists(paste0(obj_fun,'_obj'))){
  stop('objective function argument chosen does not exist, exiting')
}
#CHECK IF INST OBSERVED EXIST AND IF OBJECTIVE FUNCTION IS USING IT
logic_inst = ifelse(file.exists(file.path(basin_dir,paste0('flow_instantaneous_',basin,'.csv'))),
                    TRUE, FALSE)

if(!grepl("_NULL$", obj_fun) & !logic_inst){
  stop('cannot use instantaneous flow in objective function if data does not exist, exiting')
}

# if a fold was specified do cv unless por was specified then ignore cv
cv = if(is.na(fold) | por) FALSE else TRUE
if(is.na(fold)) por=TRUE

#CHECK IF CV FOLD ARGUMENT IS VALID
if((cv)&(!file.exists(file.path(basin_dir,paste0('forcing_validation_cv_',
                                                 fold,'_',substr(basin,1,5),'-1.csv'))))){
  stop('fold cv_',fold,' does not exist, exiting')
}

#
settings_message = ""
settings_message = paste0(settings_message,"Basin: ", basin,'\n')
settings_message = paste0(settings_message,"Input Directory: ", output_dir,'\n')
settings_message = paste0(settings_message,"Objective Function: ", obj_fun,'\n')
settings_message = paste0(settings_message,"Optimizer: ", optimizer,'\n')
settings_message = paste0(settings_message,"Iterations: ", ifelse(lite,'2500','5000'),'\n')
settings_message = paste0(settings_message,"Cores Utilized for Parallelization: ", n_cores,'\n')
settings_message = paste0(settings_message,"Cross-Validation: ", cv,'\n')
if(cv){
    settings_message = paste0(settings_message,'CV fold: ', fold,'\n')
} else{
    settings_message = paste0(settings_message,'Period of Record Run: TRUE\n')
}
message()
message(settings_message)


##########################################################################
####################Directory Setup
##########################################################################

# will reside inside of the basin dir
# will get an auto-incremented numeric postfix, ie results_01, results_02
results_dir_prefix = 'results'

if(cv) results_dir_prefix = paste0(results_dir_prefix,'_cv_',fold)
if(por) results_dir_prefix = paste0(results_dir_prefix,'_por')

# will reside inside of the results dir
plot_dir = 'plots'

# add a numeric postfix to the results dir
results_dirs = list.files(basin_dir,paste0(results_dir_prefix,'*'))
results_dir = ''; postfix_i = 0
while(results_dir == ''){
  postfix_i = postfix_i + 1
  postfix = paste0('_',sprintf('%02d',postfix_i))
  potential_dir = paste0(results_dir_prefix,postfix)
  if(!(potential_dir %in% results_dirs) | overwrite) results_dir = potential_dir
}

if(overwrite){unlink(file.path(basin_dir, results_dir), recursive = TRUE)}

output_path = file.path(basin_dir, results_dir)
plot_path = file.path(output_path, plot_dir)
dir.create(output_path, showWarnings = FALSE)
dir.create(plot_path, showWarnings = FALSE)

cat(blue$bold('Outputting results to:',output_path),'\n')

#######################################################################
####################Input Setup
######################################################################

default_pars = fread(file.path(basin_dir,'pars_default.csv'))
default_pars = default_pars[order(zone,name)]
zones = default_pars$zone |> unique() |> as.character() |> sort() |> str_subset('-')
n_zones = length(zones)
#create a dummy zone
if(n_zones==0) zones = paste0(basin,'-1')

# consuse zones
cu_zones = zones |> str_subset('-CU') |> sort()
n_cu_zones = length(cu_zones)
cu = ifelse(length(cu_zones)>0,TRUE,FALSE)

# Get limits file
limits = fread(file.path(basin_dir,paste0('pars_limits.csv')))
lower = limits$lower; names(lower) = limits$name
upper = limits$upper; names(upper) = limits$name

# Get forcings
forcing_raw = list()
for(zone in zones){
  forcing_raw[[zone]] = fread(file.path(basin_dir,paste0('forcing_por_',zone,'.csv')))
}

### Get Streamflow data
#daily flow
obs_daily = fread(file.path(basin_dir,paste0('flow_daily_',basin,'.csv')))[,-'Source']

#Add all timesteps to the daily record
obs_daily = obs_daily |>
  right_join(forcing_raw[[1]] |>
               as.data.table() |>
               group_by(year,month,day) |>
               summarise(hour=mean(hour),
                         .groups='drop') |>
               select(-hour),
             by=c('year','month','day')) |>
  as.data.table()
obs_daily = obs_daily[order(year,month,day)]


#inst flow
if(logic_inst){

  obs_inst = fread(file.path(basin_dir,paste0('flow_instantaneous_',basin,'.csv')))[,-'Source']

  #Add all timesteps to the instantaneous record
  obs_inst = obs_inst |>
    right_join(forcing_raw[[1]] |>
                 as.data.table() |>
                 select(-c("map_mm", "mat_degc","ptps")),
               by=c('year','month','day','hour')) |>
    as.data.table()

}else{obs_inst = NULL}

# flow data for upstream points to route down (if it exists)
upflow_files = list.files(basin_dir,'upflow_*',full.names = TRUE) |> sort()
n_upstream = length(upflow_files)
upflow =  NULL
if(n_upstream > 0){
  upstream_lids = gsub('upflow_','',gsub('.csv','',basename(upflow_files)))
  upflow = list()
  for(u in 1:n_upstream){
    upflow[[upstream_lids[u]]] = fread(upflow_files[u])|>
      right_join(forcing_raw[[1]] |>
                   as.data.table() |>
                   dplyr::select(-c("map_mm", "mat_degc","ptps")),
                 by=c('year','month','day','hour')) |>
      mutate(flow_cfs = vctrs::vec_fill_missing(flow_cfs, max_fill = 4)) |>
      as_tibble()
      }
}

# if we are doing cross validation, set the observations during the cv period to NA
# that way the objective function will not get computed for that period
if(cv){
  fold_tsadj_inst = read_csv(file.path(basin_dir,paste0('forcing_validation_cv_',fold,'_',zone,'.csv')),
                           col_types = cols()) |>
    select(-c("map_mm", "mat_degc","ptps")) |>
    mutate(set_na = TRUE)

  if(!is.null(obs_inst)){
    obs_inst = obs_inst |>
      full_join(fold_tsadj_inst,
                by=c('year','month','day','hour')) |>
      as.data.table()
    obs_inst[set_na==TRUE,flow_cfs:=NA]
    obs_inst$set_na = NULL
    
  }

  fold_tsadj_daily = fold_tsadj_inst |>
    group_by(year,month,day,set_na) |>
    summarise(hour=mean(hour),
              .groups='drop') |>
    select(-hour)

  obs_daily = obs_daily |>
    full_join(fold_tsadj_daily,
              by=c('year','month','day')) |>
    as.data.table()
  obs_daily[set_na==TRUE,flow_cfs:=NA]
  obs_daily$set_na = NULL
}

#If no zones, then the run is a routing reach only.  Revert dummy zone and zones back to NULL
if(n_zones == 0)
  zone = zones = NULL

#Final check if there are observation to calibrate to
if(!grepl("_NULL$", obj_fun) & sum(!is.na(obs_inst$flow_cfs))<=1){
  stop('There is no observed instantaneous flow data within calibration period, exiting')
}else if(!grepl("^NULL_", obj_fun) & sum(!is.na(obs_daily$flow_cfs))<=1){
  stop('There is no observed daily flow data within calibration period, exiting')
}

#Test objective function for any errors
tryCatch({
  obj_test = get(paste0(obj_fun,'_obj'))(obs_daily |> mutate(sim_flow_cfs = flow_cfs+1),
                                      if(logic_inst){
                                             obs_inst |> mutate(sim_flow_cfs = flow_cfs+1)
                                        }else{obs_inst})
},error = function(e){stop(paste("Objective Function had  the following error, exiting:",e))
  }
)


########################################################################
####################Optimization
#########################################################################

#Track run time
ptm = proc.time()

out = run_controller_edds(lower, upper, basin, dt_hours, default_pars, obs_daily, obs_inst,
                         forcing_raw, upflow, obj_fun, n_zones,cu_zones,n_cores,lite)

p_optimal = out["p_best"][[1]]
#dput(p_optimal,file=file.path(output_path,'p_optimal_dds.txt'))
f_trace = out["f_trace"][[1]]
p_trace = out["p_trace"][[1]]

seeds = sapply(out[["seeds"]], function(seed) paste(seed, collapse = ","))

if(length(unique(seeds))!=length(seeds)){
  stop('Parallel process not properly being seeded, exiting, please notify Geoffrey')
}

#Can uncomment in case of error with non unique seeds
#cat("Unique seeds:", length(unique(seeds)), "\n")
#cat("Total workers:", length(seeds), "\n")

#############################################################################
####################Create output plot
#############################################################################

f_trace$iter <- 1:nrow(f_trace)
p_t <- ggplot(f_trace[-(1:100), ]) +
  geom_point(aes(iter, obj_fun)) +
  labs(main = "Auto Calibration Evolution", x = "Iteration", y = "Objective Function Score") +
  theme_bw()
ggsave(sprintf("%s/objfun_trace.pdf", plot_path), p_t, width = 10, height = 8)

p_trace_plot=p_trace |>
  mutate(f_trace) |>
  as.data.table() |>
  melt.data.table(id.vars = c('iter','obj_fun'))

p_trace_plot = p_trace_plot |> lazy_dt() |> group_by(variable) |>
  mutate(value=(value-min(value))/(max(value)-min(value))) |>
  ungroup() |> as.data.table()

ptraj = ggplot(p_trace_plot)+
  geom_line(aes(iter,value,color=variable))+
  #geom_text(aes(Iter,value,label=variable),data=traj[Iter==max(traj$Iter)],
  #          hjust=1,size=3)+
  theme_minimal()+
  theme(legend.position = 'best') +
  scale_color_manual(values=colorRampPalette(colorblind_pal()(8))(length(lower)))
# ptraj
ggsave(sprintf('%s/trajectories_best.pdf',plot_path),ptraj,width=10,height=8)

##########################################################################
####################Write Outputfiles
##########################################################################

###Get optimized parameters
optimal_pars = update_params(p_optimal,names(lower),default_pars)
if(n_zones>0){
  forcing = fa_nwrfc(dt_hours, forcing_raw, optimal_pars)
}else{
  forcing = forcing_raw
}

# if there are any consuse zones, replace those params with the zone 1 params
if(cu){
  zone1_name = names(forcing)[1]
  optimal_pars = update_cu_params(optimal_pars,zone1_name,cu_zones)
}
optimal_pars = optimal_pars[order(zone,name)]

#get optimized daily simulation
optimal_sim_daily = model_wrapper(p_optimal,names(lower),dt_hours,optimal_pars,obs_daily,obs_inst,
                                  forcing_raw,upflow,obj_fun,n_zones,cu_zones,return_flow =TRUE)

#write csv output
write.csv(optimal_pars,file=file.path(output_path,'pars_optimal.csv'),row.names = FALSE,quote=FALSE)
write.csv(f_trace,file=file.path(output_path,'objfun_trace.csv'),row.names = FALSE,quote=FALSE)
write.csv(p_trace,file=file.path(output_path,'pars_trace.csv'),row.names = FALSE,quote=FALSE)

# stats from optimal run
optimal_gof = with(optimal_sim_daily,gof(sim_flow_cfs,flow_cfs))
end = (proc.time() - ptm)['elapsed']
hr = floor(end/60/60)
min = floor((end-hr*60*60)/60)
sec = round(end-hr*60*60-min*60)
settings_message = paste0(settings_message,"Run Time: ", round(end/60,1),"min",'\n')

# run summary file
file.path(output_path,'run_summary.txt') |> sink()
cat('Objective Function Score:',
      as.character(round(tail(f_trace[,'obj_fun'],n=1),3)))# |> print
cat('\n','\n','--STATISTICS-- ','\n')
optimal_gof |> print()
cat('\n','--OPTIMAL PARAMETERS-- ','\n')
options(crayon.enabled = FALSE)
tibble::tibble(names(p_optimal),p_optimal) |> arrange(names(p_optimal)) |> print(n=Inf,)
options(crayon.enabled = TRUE)
sink()

# run settings file
cat(settings_message,file=file.path(output_path,'run_settings.txt'))
file.path(output_path,'run_settings.txt') |> sink(append=TRUE)
cat('\n','--PARAMETER LIMITS--','\n')
as.matrix(limits |> arrange(zone,name) |> as.data.table()) |> print()
sink()

##############################################################################
####################Final Screen Output
#############################################################################

tibble::tibble(names(p_optimal),p_optimal) |> print(n=Inf)
optimal_gof |> print()

cat(bgCyan('Run took ',hr,' hr ',min,' min ',sec,' sec'),'\n')