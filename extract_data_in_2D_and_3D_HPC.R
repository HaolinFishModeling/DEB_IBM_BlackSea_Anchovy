# extract data from preprocessed data 20250120

# 1. extract data in 2D---------------
# functions to extract data frame for the whole grid
extract_current_Temp_2D<-function(updated_date){
  updated_date <- as.Date(updated_date)
  path1<-"/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/Temp/"
  data <- readRDS(paste0(path1,'Temp_',updated_date,'.rds'))
  
  return(data)
}

extract_current_Food_2D<-function(updated_date){
  updated_date <- as.Date(updated_date)
  path1<-"/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/Food/"
  data <- readRDS(paste0(path1,'Food_',updated_date,'.rds'))
  
  return(data)
}

extract_next_Temp_2D<-function(updated_date){
  updated_date <- as.Date(updated_date)
  if(updated_date == as.Date('2023-12-31')){ # last day for hindcast simulation
    next_date <- updated_date
  }else{
    next_date <- updated_date + days(1)
  }
  
  path1<-"/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/Temp/"
  data <- readRDS(paste0(path1,'Temp_',next_date,'.rds'))
  
  return(data)
}

extract_next_Food_2D<-function(updated_date){
  updated_date <- as.Date(updated_date)
  if(updated_date == as.Date('2023-12-31')){ # last day for hindcast simulation
    next_date <- updated_date
  }else{
    next_date <- updated_date + days(1)
  }
  path1<-"/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/Food/"
  data <- readRDS(paste0(path1,'Food_',next_date,'.rds'))
  
  return(data)
}

# extract data for PTM 
extract_current_uo_2D<-function(updated_date){ # always hourly
  updated_date <- as.Date(updated_date)
  path1<-paste0('/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/uo/')
  data <- readRDS(paste0(path1,'uo_',updated_date,'.rds'))
  
  return(data)
}

extract_next_uo_2D<-function(updated_date,var_name){ # always hourly
  updated_date <- as.Date(updated_date)
  if(updated_date == as.Date('2023-12-31')){ # last day for hindcast simulation
    next_date <- updated_date
  }else{
    next_date <- updated_date + days(1)
  }
  path1<-paste0('/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/uo/')
  data <- readRDS(paste0(path1,'uo_',next_date,'.rds'))
  
  return(data)
}

extract_current_vo_2D<-function(updated_date){ # always hourly
  updated_date <- as.Date(updated_date)
  path1<-paste0('/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/vo/')
  data <- readRDS(paste0(path1,'vo_',updated_date,'.rds'))
  
  return(data)
}

extract_next_vo_2D<-function(updated_date,var_name){ # always hourly
  updated_date <- as.Date(updated_date)
  if(updated_date == as.Date('2023-12-31')){ # last day for hindcast simulation
    next_date <- updated_date
  }else{
    next_date <- updated_date + days(1)
  }
  path1<-paste0('/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/vo/')
  data <- readRDS(paste0(path1,'vo_',next_date,'.rds'))
  
  return(data)
}

extract_current_sd_2D<-function(updated_date){  # always hourly
  updated_date <- as.Date(updated_date)
  path1<-'/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/Stokes_Drift/'
  potential.sd <- readRDS(paste0(path1,'Stokes_drift_mean_daily_data_',updated_date,'.rds'))
  
  return(potential.sd)
}

extract_next_sd_2D<-function(updated_date){  # always hourly
  updated_date <- as.Date(updated_date)
  if(updated_date == as.Date('2023-12-31')){ # last day for hindcast simulation
    next_date <- updated_date
  }else{
    next_date <- updated_date + days(1)
  }
  path1<-'/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/Stokes_Drift/'
  potential.sd.next <- readRDS(paste0(path1,'Stokes_drift_mean_daily_data_',next_date,'.rds'))
  
  return(potential.sd.next)
}

# 2. extract full grid in 3D ----------------
# we don't need depth level here, because we only need to consider the full grid on the initial of the day
extract_current_3D<-function(updated_date,var_name){
  # var_name = c(MES, MIC,Temp,uo,vo,wo,Sal,rho)
  updated_date <- as.Date(updated_date)
  path1<-paste0('/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/3D_all_data_frame_fst/',var_name,'/',var_name,'_',updated_date,'.fst')
  potential.var <- read_fst(path1)
  
  return(potential.var)
}

extract_next_3D<-function(updated_date,depth_level,var_name){
  # var_name = c(MES, MIC,Temp,uo,vo,wo,Sal,rho)
  updated_date <- as.Date(updated_date)
  next_date<-updated_date+days(1)
  path1<-paste0('/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/3D_all_data_frame_fst/',var_name,'/',var_name,'_',next_date,'.fst')
  potential.var <- read_fst(path1)
  
  return(potential.var)
}

define_depth_level<-function(depth){
  
  index <- findInterval(depth, all_depth_levels_w, rightmost.closed = TRUE)
  return(index)
}
