#!/usr/bin/env Rscript

xfun::pkg_load2(c('magrittr','dplyr', 'data.table', 'dtplyr','hydroGOF',
                   'lubridate','readr','tibble','ggthemes','stringr',
                   'gtable','gridExtra','knitr','kableExtra','rfchydromodels',
                   'tidyr','egg', 'git2r','reshape2','ncdf4','sf','ggmap',
                   'mapproj', 'plotly'))
select <- dplyr::select
xfun::pkg_attach2('ggplot2')
import::from(dplyr, filter, select, summarise, group_by, ungroup, mutate,left_join,
             slice_max, rename, pull,rename,mutate_if,n)
import::from(plotly, ggplotly)
import::from(stringr, str_detect, str_subset)
import::from(argparser,arg_parser,add_argument,parse_args)
import::from(crayon, inverse, green)
import::from(data.table,rbindlist)
import::from(readr,read_csv)
import::from(hydroGOF,NSE,pbias,rPearson,KGE)
import::from(tibble, tibble, as_tibble)
import::from(tidyr, pivot_longer, pivot_wider)
import::from(dplyr, bind_rows)
import::from(parallel,detectCores,makeCluster,clusterSetRNGStream,
             clusterCall,stopCluster)

parser = arg_parser("Create CV Plots", hide.opts = TRUE)
# by default ArgumentParser will add an help option
parser = add_argument(parser, "--dir", help = "Input directory path containing basin directories")
parser = add_argument(parser, "--basins", default=NA_character_, help = "Basins to run", nargs=Inf)
parser = add_argument(parser, "--cleanup",flag=TRUE, help = "Option to delete results directories which are not used for CV analysis")
args = parse_args(parser)

#args = list(dir='runs/opt_run/1zone',basins='TLMO3',cleanup=FALSE)

runs_dir = args$dir
basins = args$basins
delete_dir = args$cleanup

#directory to save the figures
plot_dir = file.path(runs_dir, 'cv_plots')

#inside of the basin dir
results_dir_prefix = 'results'

# check if the user specified the basins, if not run them all
basins = if(any(is.na(basins))){
  list.dirs(runs_dir, full.names = FALSE, recursive = FALSE)
}else{
  basins
}

# loop through each basin
for(basin in basins){

  cat(inverse$blue(basin),'\n')

  #####################################
  #Get list of cv and por results
  ####################################

  basin_dir = file.path(runs_dir,basin)
  results_dirs = list.files(basin_dir,paste0(results_dir_prefix,'*'))

  #detect if there are any cv run results, if not skip loop
  if(any(str_detect(results_dirs,'cv'))){
    cv_results_dir = results_dirs |> str_subset('cv')
  }else{
    cv_results_dir = NULL
  }

  #for each cv run results, make sure the needed files are available, in not
  #remove from the analysis
  if(!is.null(cv_results_dir)){
    for(i in length(cv_results_dir):1){
      if(!all(c('optimal_daily.csv',
                'validation_daily.csv',
                'optimal_6hr_inst.csv',
                'validation_6hr_inst.csv',
                'run_settings.txt') %in% 
              list.files(file.path(basin_dir,cv_results_dir[i]),'*'))){
        cat('\t',paste0('Not using ',cv_results_dir[i],' as no postprocessing files'),'\n')
        cv_results_dir = cv_results_dir[-i]
      }
    }
  }

  #check if there are any valid cv run results, if not skip loop
  if(length(cv_results_dir)>0){
    unique_cv = unique(substr(gsub('results_cv_','',cv_results_dir), start = 1 , stop = 1 ))
  }else{
    cat(inverse$red('\t',paste0('!!Skipping ',basin,' as no CV runs available!!')),'\n')
    unique_cv = cv_results_dir = NULL
    next
  }

  #Add a blank print out row
  cat('\n')

  #detect if there are any cv run results, if not skip loop
  if(any(str_detect(results_dirs,'por'))){
    por_results_dir = results_dirs |> str_subset('por')
  }else{
    por_results_dir = NULL
  }

  #for each por run results, make sure the needed file is available, in not
  #remove from the analysis
  if(!is.null(por_results_dir)){
    for(j in length(por_results_dir):1){
      if(!all(c('optimal_daily.csv',
                'optimal_6hr_inst.csv',
                'run_settings.txt') %in% 
              list.files(file.path(basin_dir,por_results_dir[i]),'*'))){
        cat('\t',paste0('Not using ',por_results_dir[j],' as no postprocessing files'),'\n')
        por_results_dir = por_results_dir[-j]
      }
    }
  }

  #detect if there are any por run results, stop run if not
  if(length(por_results_dir)==0){
    cat(inverse$red('\t',paste0('!!Skipping ',basin,' as no POR runs available!!')),'\n')
    por_results_dir = NULL
    next
  }
  
  obj_fun_check = c()
  #detect if all the cv and por runs are using the same objective function
  for(check_dir in c(cv_results_dir,por_results_dir)){
    run_setting_check = read.table(file.path(basin_dir, check_dir, 'run_settings.txt'),
                                   sep = ':',nrows=9,col.names=c('Run','Setting'))[-(1:2),]
    obj_fun_check = c(obj_fun_check, 
                      run_setting_check[run_setting_check$Run=='Objective Function','Setting'])
  }
  
  if(length(unique(obj_fun_check))>1){
    cat(inverse$red('\t',paste0('!!Skipping ',basin,' because not all run used the same objective function!!')),'\n')
    next
  }
  
  #Check if objective function is using subdaily flow
  logic_inst = ifelse(!grepl("_NULL$", unique(obj_fun_check)),TRUE,FALSE)

  #########################################
  #Identify best performing cv and por runs
  ########################################

  ###CV Runs ####
  cv_results_runs_list = list()
  #loop through each cv fold
  for(fold in unique_cv){

    #get all runs associated with the current fold
    fold_results_dir = cv_results_dir |> str_subset(paste0('cv_',fold))

    #loop through each run associated with each fold
    for(fold_run in fold_results_dir){
      #open the calibration period simulation file
      calib = read_csv(file.path(basin_dir,fold_run,'optimal_daily.csv'), show = FALSE)[-(1:366),] |>
        select(-imputed) |>
        mutate(period = 'calb',
               run = fold_run,
               fold = fold)
      #open the validation period simulation file
      valid = read_csv(file.path(basin_dir,fold_run,'validation_daily.csv'), show = FALSE) |>
        select(-imputed) |>
        mutate(period = 'valid',
               run = fold_run,
               fold = fold)
      #if cv fold 1, remove the warm up period
      if(fold == "1"){
        valid = valid[-(1:366),]
      }
      #bind the calibration and validation simulations to a list
      cv_results_runs_list[[fold_run]] = rbind(calib,valid)
    }
    #combine all output in the master list to a single datatable
    cv_results_runs = rbindlist(cv_results_runs_list)
  }

  #calculate the combo metric by month, run, fold, period.
  #sum by month combining the calb/valid peiods
  #then for each fold find the run with the max score
  cv_results_best = cv_results_runs |>
    as_tibble() |>
    group_by(month,run,fold,period) |>
    summarise(kge = KGE(sim_flow_cfs,flow_cfs),
          nkge = 1/(2-kge), .groups='drop') |>
    group_by(run,fold) |>
    summarise(nkge_sum = sum(nkge), .groups='drop') |>
    group_by(fold) |>
    slice_max(nkge_sum, n = 1) |>
    pull(run)

  #Get a table of only the best cv daily sim results
  cv_results_runs_best = cv_results_runs |>
    as_tibble() |>
    filter(run %in% cv_results_best)

  cat('\t Using the following CV best runs:  ')
  cat(paste(cv_results_best), sep=', ')
  cat('\n')

  #Optional delete directories which aren't tagged as best
  if(delete_dir){
    delete_dir_list = setdiff(cv_results_dir,cv_results_best)
    for (d_dir in delete_dir_list){
      unlink(file.path(basin_dir,d_dir), recursive = TRUE)
    }
  }

  ###POR Runs ####
  por_results_runs_list = list()
  #loop through each por
  for(por_run in por_results_dir){
    #open the calibration period simulation file
    calib = read_csv(file.path(basin_dir,por_run,'optimal_daily.csv'), show = FALSE)[-(1:366),] |>
      select(-imputed) |>
      mutate(period = 'calb',
             run = por_run)
    #bind the calibration and validation simulations to a list
    por_results_runs_list[[por_run]] = calib
  }

  #combine all output in the master list to a single datatable
  por_results_runs = rbindlist(por_results_runs_list)

  #calculate the combo metric by month, run, fold, period.
  #sum by month combining the calb/valid peiods
  #then for each fold find the run with the max score
  por_results_best = por_results_runs |>
    as_tibble() |>
    group_by(month,run) |>
    summarise(kge = KGE(sim_flow_cfs,flow_cfs),
              nkge = 1/(2-kge), .groups='drop') |>
    group_by(run) |>
    summarise(nkge_sum = sum(nkge), .groups='drop') |>
    slice_max(nkge_sum, n = 1) |>
    pull(run)

  #Get a table of only the best cv daily sim results
  por_results_runs_best = por_results_runs |>
    as_tibble() |>
    filter(run %in% por_results_best) |>
    mutate(wyear=ifelse(month>=10,year+1,year))

  cat('\t Using the following POR best run:  ',por_results_best, '\n\n')

  #Optional delete directories which aren't tagged as best
  if(delete_dir){
    delete_dir_list = setdiff(por_results_dir,por_results_best)
    for (d_dir in delete_dir_list){
      unlink(file.path(basin_dir,d_dir), recursive = TRUE)
    }
  }
  
  #########################################
  #Format 95th percentile stats
  ########################################
  
  #Identify if we are using the daily or instantaneous observed data csv
  peak_file = file.path(basin_dir, paste0(ifelse(logic_inst,
                                                 'flow_instantaneous_',
                                                 'flow_daily_'),basin,'.csv'))
  
  #Get teh 95th percentile threshold
  peak_th =  read_csv(peak_file, show = FALSE) |>
    as_tibble() |>
    mutate(wyear=ifelse(month>=10,year+1,year)) |>
    filter(wyear > 1981) |>
    pull('flow_cfs') |>
    quantile(.95,na.rm = TRUE,names = FALSE)
  
  
  if(logic_inst){
    #If using instantaneous data for the 95th percentile, load the cv and por
    #best run's *_6hr_inst.csv files
    cv_results_runs_best_95th_list = list()
    for(k in 1:length(unique_cv)){
      fold = unique_cv[k]
      fold_run =  cv_results_best[k]
      #open the subdaily cv file
      subdaily_valid = read_csv(file.path(basin_dir,fold_run,'validation_6hr_inst.csv'), show = FALSE) |>
        mutate(period = 'valid',
               run = fold_run,
               fold = fold,
               wyear=ifelse(month>=10,year+1,year)) |>
        filter(wyear > 1981 & flow_cfs>=peak_th) |>
        select(-wyear)
      #bind the validation instantaneous data from the best CV runs
      cv_results_runs_best_95th_list[[fold_run]] = subdaily_valid
    }
    #combine all the instantaneous cv data to a single datatable
    cv_results_runs_best_95th = rbindlist(cv_results_runs_best_95th_list)
    
    #get instantaneous por data
    por_results_runs_best_95th= read_csv(file.path(basin_dir,por_results_best,'optimal_6hr_inst.csv'), show = FALSE) |>
      mutate(period = 'calb',
             run = por_results_best,
             wyear=ifelse(month>=10,year+1,year)) |>
      filter(wyear > 1981 & flow_cfs>=peak_th) 
  }else{
    #If using daily data for the 95th percentile, use the already loaded and
    #formated cv and por best runs
    cv_results_runs_best_95th_list = NULL
    cv_results_runs_best_95th  = cv_results_runs_best |>
      filter(flow_cfs>=peak_th)
    por_results_runs_best_95th = por_results_runs_best |>
      filter(flow_cfs>=peak_th)
  }
  
  ##############################################
  #Calculate CV NNSE, Pbias, R2 and KGE scores
  #############################################

  #Metrics using all daily data
  cv_all_daily = cv_results_runs_best |>
                  filter(period == 'valid') |>
                  group_by(fold) |>
                  summarise(nse = NSE(sim_flow_cfs,flow_cfs),
                            pbias = pbias(sim_flow_cfs,flow_cfs),
                            kge = KGE(sim_flow_cfs,flow_cfs),
                            nnse = 1/(2-nse),
                            npbias = 1-abs(pbias)/100,
                            nkge = 1/(2-kge),
                            R2 = rPearson(sim_flow_cfs,flow_cfs)^2,.groups='drop') |>
                  select(-c(nse,pbias,kge)) |>
                  pivot_longer(-c(fold),names_to ='metric') |>
                  mutate(type = 'CV \n run',
                         temporal = 'All Daily')

  #Metrics for monthly binned statistics
  cv_monthly_binned = cv_results_runs_best |>
    filter(period == 'valid') |>
    group_by(fold,month) |>
    summarise(nse = NSE(sim_flow_cfs,flow_cfs),
              pbias = pbias(sim_flow_cfs,flow_cfs),
              kge = KGE(sim_flow_cfs,flow_cfs),
              nnse = 1/(2-nse),
              npbias = 1-abs(pbias)/100,
              nkge = 1/(2-kge),
              R2 = rPearson(sim_flow_cfs,flow_cfs)^2,.groups='drop') |>
    select(-c(nse,pbias,kge)) |>
    group_by(fold) |>
    summarise(nnse = sum(nnse)/12,
              npbias = sum(npbias)/12,
              nkge = sum(nkge)/12,
              R2 = sum(R2)/12,.groups='drop') |>
    pivot_longer(-c(fold),names_to ='metric') |>
    mutate(type = 'CV \n run',
           temporal = 'Monthly Binned')

  #Metrics for April-September water supply statistics
  cv_apr_sep_ws = cv_results_runs_best |>
    filter(month %in% 4:9 & period == 'valid') |>
    group_by(fold,year) |>
    summarise(volume=sum(flow_cfs*86400)*2.29569e-5/1000,
              sim_volume=sum(sim_flow_cfs*86400)*2.29569e-5/1000,
              .groups='drop') |>
    group_by(fold) |>
    summarise(nse = NSE(sim_volume,volume),
              pbias = pbias(sim_volume,volume),
              kge = KGE(sim_volume,volume),
              nnse = 1/(2-nse),
              npbias = 1-abs(pbias)/100,
              nkge = 1/(2-kge),
              R2 = rPearson(sim_volume,volume)^2,.groups='drop') |>
    select(-c(nse,pbias,kge)) |>
    pivot_longer(-c(fold),names_to ='metric') |>
    mutate(type = 'CV \n run',
           temporal = 'Seasonal Volume\n Apr-Sep')
    
  #Metrics for january-july water supply statistics
  cv_jan_jul_ws = cv_results_runs_best |>
    filter(month %in% 1:7 & period == 'valid') |>
    group_by(fold,period,year) |>
    summarise(volume=sum(flow_cfs*86400)*2.29569e-5/1000,
              sim_volume=sum(sim_flow_cfs*86400)*2.29569e-5/1000,
              .groups='drop') |>
    group_by(fold) |>
    summarise(nse = NSE(sim_volume,volume),
              pbias = pbias(sim_volume,volume),
              kge = KGE(sim_volume,volume),
              nnse = 1/(2-nse),
              npbias = 1-abs(pbias)/100,
              nkge = 1/(2-kge),
              R2 = rPearson(sim_volume,volume)^2,.groups='drop') |>
    select(-c(nse,pbias,kge)) |>
    pivot_longer(-c(fold),names_to ='metric') |>
    mutate(type = 'CV \n run',
           temporal = 'Seasonal Volume\n Jan-Jul')
  
  #Need a try/catch for cases where inst data is used and there is none available
  #for the cv period
  tryCatch({
    cv_q95 = cv_results_runs_best_95th |>
      filter(period == 'valid') |>
      group_by(fold) |>
      summarise(nse = NSE(sim_flow_cfs,flow_cfs),
                pbias = pbias(sim_flow_cfs,flow_cfs),
                kge = KGE(sim_flow_cfs,flow_cfs),
                nnse = 1/(2-nse),
                npbias = 1-abs(pbias)/100,
                nkge = 1/(2-kge),
                R2 = rPearson(sim_flow_cfs,flow_cfs)^2,.groups='drop') |>
      select(-c(nse,pbias,kge)) |>
      pivot_longer(-c(fold),names_to ='metric') |>
      mutate(type = 'CV \n run',
             temporal = paste('Q>.95 Highflows\n',ifelse(logic_inst,'Instantaneous','Daily')))
  },error = function(e){cv_q95 <<-NULL
  })
    
  #combine all cv validation metric performance scores for different temporal periods
  cv_perfomance = rbind(cv_all_daily, cv_monthly_binned, cv_apr_sep_ws,
                     cv_jan_jul_ws,cv_q95)

  ##############################################
  #POR run Stationary Bootstrappings
  #############################################
  
  #Filter out non complete wy
  wy_population = por_results_runs_best |> 
    group_by(wyear) |> 
    summarise(count=n(),.groups='drop') |>
    filter(360<count) |>
    pull(wyear)
  
  cv_avg_length = cv_results_runs_best |> 
    mutate(wyear=ifelse(month>=10,year+1,year)) |> 
    filter(period == 'valid' & wyear %in% wy_population ) |>
    filter(month == 1 & day == 1) |>
    group_by(fold) |>
    summarise(wy_count=n(),.groups='drop') |>
    summarise(avg_wy_len = round(mean(wy_count),0)) |>
    pull(avg_wy_len)
 
  #Stationary Bootstrap Looper Function ######
  bs_looper = function(sim, sim_95th, wy_pool, sample_size,inst_logic,bs_iter=1000){
  
    bs_performance_worker = list()
    
    for(i in 1:bs_iter){
      #Get random sample of 10 wy
      sample_wy = sample(wy_pool,sample_size)
      sample_sim = sim |> filter(wyear %in% sample_wy)
      sample_sim_95th = sim_95th |> filter(wyear %in% sample_wy)
      
      #Metrics using all daily data
      all_daily = sample_sim |>
        group_by(period) |>
        summarise(nse = NSE(sim_flow_cfs,flow_cfs),
                  pbias = pbias(sim_flow_cfs,flow_cfs),
                  kge = KGE(sim_flow_cfs,flow_cfs),
                  nnse = 1/(2-nse),
                  npbias = 1-abs(pbias)/100,
                  nkge = 1/(2-kge),
                  R2 = rPearson(sim_flow_cfs,flow_cfs)^2,.groups='drop') |>
        select(-c(nse,pbias,kge)) |>
        pivot_longer(-c(period),names_to ='metric') |>
        mutate(type = 'Stationary \nBootstrapping',
               temporal = 'All Daily')
      
      #Metrics for monthly binned statistics
      monthly_binned = sample_sim |>
        group_by(period,month) |>
        summarise(nse = NSE(sim_flow_cfs,flow_cfs),
                  pbias = pbias(sim_flow_cfs,flow_cfs),
                  kge = KGE(sim_flow_cfs,flow_cfs),
                  nnse = 1/(2-nse),
                  npbias = 1-abs(pbias)/100,
                  nkge = 1/(2-kge),
                  R2 = rPearson(sim_flow_cfs,flow_cfs)^2,.groups='drop') |>
        select(-c(nse,pbias,kge)) |>
        group_by(period) |>
        summarise(nnse = sum(nnse)/12,
                  npbias = sum(npbias)/12,
                  nkge = sum(nkge)/12,
                  R2 = sum(R2)/12,.groups='drop') |>
        pivot_longer(-c(period),names_to ='metric') |>
        mutate(type = 'Stationary \nBootstrapping',
               temporal = 'Monthly Binned')
      
      #Metrics for April-September water supply statistics
      apr_sep_ws = sample_sim |>
        filter(month %in% 4:9) |>
        group_by(period,year) |>
        summarise(volume=sum(flow_cfs*86400)*2.29569e-5/1000,
                  sim_volume=sum(sim_flow_cfs*86400)*2.29569e-5/1000,
                  .groups='drop') |>
        group_by(period) |>
        summarise(nse = NSE(sim_volume,volume),
                  pbias = pbias(sim_volume,volume),
                  kge = KGE(sim_volume,volume),
                  nnse = 1/(2-nse),
                  npbias = 1-abs(pbias)/100,
                  nkge = 1/(2-kge),
                  R2 = rPearson(sim_volume,volume)^2,.groups='drop') |>
        select(-c(nse,pbias,kge)) |>
        pivot_longer(-c(period),names_to ='metric') |>
        mutate(type = 'Stationary \nBootstrapping',
               temporal = 'Seasonal Volume\n Apr-Sep')
      
      #Metrics for january-july water supply statistics
      jan_jul_ws = sample_sim |>
        filter(month %in% 1:7) |>
        group_by(period,year) |>
        summarise(volume=sum(flow_cfs*86400)*2.29569e-5/1000,
                  sim_volume=sum(sim_flow_cfs*86400)*2.29569e-5/1000,
                  .groups='drop') |>
        group_by(period) |>
        summarise(nse = NSE(sim_volume,volume),
                  pbias = pbias(sim_volume,volume),
                  kge = KGE(sim_volume,volume),
                  nnse = 1/(2-nse),
                  npbias = 1-abs(pbias)/100,
                  nkge = 1/(2-kge),
                  R2 = rPearson(sim_volume,volume)^2,.groups='drop') |>
        select(-c(nse,pbias,kge)) |>
        pivot_longer(-c(period),names_to ='metric') |>
        mutate(type = 'Stationary \nBootstrapping',
               temporal = 'Seasonal Volume\n Jan-Jul')
      
      #Need a try/catch for cases where inst data is used and there is none available
      #for sampled water years
      tryCatch({
        q95 = sample_sim_95th |>
          group_by(period) |>
          summarise(nse = NSE(sim_flow_cfs,flow_cfs),
                    pbias = pbias(sim_flow_cfs,flow_cfs),
                    kge = KGE(sim_flow_cfs,flow_cfs),
                    nnse = 1/(2-nse),
                    npbias = 1-abs(pbias)/100,
                    nkge = 1/(2-kge),
                    R2 = rPearson(sim_flow_cfs,flow_cfs)^2,.groups='drop') |>
          select(-c(nse,pbias,kge)) |>
          pivot_longer(-c(period),names_to ='metric') |>
          mutate(type = 'Stationary \nBootstrapping',
                 temporal = paste('Q>.95 Highflows\n',ifelse(inst_logic,'Instantaneous','Daily')))
      },error = function(e){q95 <<-NULL
      })
      
      #combine all cv validation metric performance scores for different temporal periods
      bs_iter_perfomance = rbind(all_daily, monthly_binned, apr_sep_ws,
                            jan_jul_ws,q95)
      
      bs_performance_worker = rbind(bs_performance_worker, bs_iter_perfomance)
    }
  return(bs_performance_worker)
  }
  
  #Parallel Registration
  n_cores = ifelse(detectCores()-2>=8,8,detectCores()-2)
  my_cluster = makeCluster(n_cores,type='FORK') #'PSOCK'
  # Seeding using L'Ecuyer-CMRG, isseed = NULL doc states initializes random process
  RNGkind("L'Ecuyer-CMRG")
  clusterSetRNGStream(cl = my_cluster, iseed = NULL)
  # Set a fresh random seed for the master process
  set.seed(sample.int(.Machine$integer.max, 1))
  (stream <- .Random.seed)
  
  cat('\t Performing Stationary Bootstrapping..... ')
  
  bs_parallel = clusterCall(cl = my_cluster,
                              bs_looper,
                              sim = por_results_runs_best,
                              sim_95th = por_results_runs_best_95th,
                              wy_pool = wy_population,
                              sample_size = cv_avg_length,
                              inst_logic = logic_inst,
                              bs_iter = 1000)
    
    
  stopCluster(cl = my_cluster)
  
  bs_performance= bind_rows(bs_parallel) |> 
    select(-period) |>
    mutate(fold = 'BS')
  
  cat('Done \n')
  
  #combine all cv runs and stationary bootstrapping performance results
  master_performance  = rbind(bs_performance,cv_perfomance)

 
  #change columns from character to factors
  #master_degrade$type = as.factor(master_degrade$type)
  master_performance =  master_performance |> mutate_if(is.character,as.factor)
  master_performance$metric = factor(master_performance$metric,levels=c('nnse','npbias','R2','nkge'))

  ##############
  #Plot figure
  ###############

  #check if the plot_dir exists, if not create
  if(!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings = FALSE)

  g = ggplot() +
    geom_violin(data=master_performance |> filter(fold == 'BS'),aes(x=temporal,y=value,color=type),fill='white') +
    geom_jitter(data=master_performance |> filter(fold != 'BS'),
                aes(x=temporal,y=value,fill=fold),height=0,width=.05,shape=21,size=2) +
    scale_color_manual("",values=rep('black',5)) +
    facet_wrap(vars(metric),ncol=1,scales = 'free_y') +
    ylab(expression(paste("Metric Score")))+xlab('') + ggtitle(basin) +
    labs(caption = paste(c(cv_results_best,por_results_best),collapse=', ')) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          axis.ticks.x.bottom=element_line(colour = "black", linewidth=.5),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 10))

  filename=paste0(basin,'_CV_Plot.png')
  ggsave(filename,g,path = plot_dir,width = 11, height = 9,bg = 'white')
  cat(green(paste('\t',file.path(plot_dir,filename), 'saved', sep=' '),'\n'))
 }
