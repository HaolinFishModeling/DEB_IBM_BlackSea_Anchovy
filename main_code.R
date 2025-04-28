# Population model construction 20250306
#############################################
####
###     load in basic package
####
########################################

# deal with raster grid
library(raster)
library(sf)
library(lubridate)
library(tidyverse)
library(dplyr)
library(terra)
library(sp)
library(fst)

# parallel calculation
library(parallel)

#############################################
####
###     Determine basic setting 
####
########################################
input.dimension <- 1
# 【very important】
# 1 = 2D
# 2 = 3D

input.data.location<- 2
# 【very important】
# 1 = all Temp, food data are from cell center without interpolation but PTM still uses interpolation
# 2 = all Temp, food, and PTM data are using interpolation at the specific location of the super individuals

input.func <- 1
# 【very important】
# 1 = DEB with Kinesis model
# 2 = only DEB model
# 3 = only Kinesis model

input.stoke.drift<-1
# 【very important】
# 1 = include stokes drift - wave effects
# 2 = without stokes drift - wave effects

input.DD.mortality<-2
# 【very important】
# 1 = add D-D mortality only for late-larvae
# 2 = add D-D mortality only for juvenile
# 3 = add D-D mortality for both late-larvae and juvenile

# add warning information
if(input.dimension==1 & input.stoke.drift == 1){
  warning('Wave effect may be incorrect for 2D dimension')
}

# define life stages in different models
life_stages_in_PTM<-c('egg','yolksac_larvae','early_larvae')
life_stages_in_Kinesis<-c('late_larvae','juvenile','adult')
life_stages_in_DEB<-c('early_larvae','late_larvae','juvenile','adult')
life_stages_in_Stokes_drift<-c('egg','yolksac_larvae')

# set up time step for simulation-------
#【vert important】
# 1/24 = daily, 1 = hourly, 2 = 30 min, 3 = 20 min, 4 = 15 min, 5 = 12 min
timestep.set <- 1 / 24 # daily

# determine reference density--------
# 2D
Db_LL<-0.34             # ind./m3 
Db_juv<-0.0033*10       # ind./m3

# 3D
Db_LL_3D<-Db_LL*30             # ind./m2 projection to 2D
Db_juv_3D<-Db_juv*100       # ind./m2 projection to 2D

# determine mass density of particle------------------
# The mass of sea density is 1013~1018 kg/m3
# if we set up 1000 for particle, it will float at the surface
# if we set up 1020 for particle, it will go down to the sea floor
# I think it is the floating eggs, so I set up rho_particle == 1000
rho_particle <- 1000  #  kg/m3,

# set up simulation duration, all of the simulation start from start_date (00:00:00) to end_date (23:00:00)-----
start_date <- ymd_hms("1990-01-01 00:00:00")
end_date <- ymd_hms("1994-12-31 23:00:00")
year_spinup<-'1990'  # run five years for spin-up to tell model generate new eggs by model themselves

# Total running years and days (running_days did not include spin-up days)
years <- c(rep(year_spinup,times=5),year(start_date):year(end_date))
running_days<-as.numeric(difftime(as.Date(end_date), as.Date(start_date))) + 1 # formal simulation days without counting spin-up days
spinup_days_in_year<-ifelse(leap_year(as.numeric(year_spinup)), 366, 365)

# set up number of super individuals (SI) and parallel cores----------
number_individual <- 4800

# *******!!!!!!no need too much CPUs because they will cost lots of time to allocate data frame to each CPU (120 CPUs is not a good idea)
# However, we used fst library, it will use more available CPU to read in dataframe.
# When submitting jobs, it is better make the available CPUs 2 times of the number_core
# i.e., number_cores == 30, BATCH CPU == 60;
# i.e., number_cores == 60, BATCH CPU == 120
number_cores <- 30  # this value MUST be fully divided by "number_individual". 

individual_per_core <- number_individual/number_cores

split_indices<-split(1:number_individual,rep(1:number_cores,each=individual_per_core))

number_SI_for_each_region <- 5      # this is related to how to set up spawning events

# create cluster
cl <- makeCluster(number_cores)

# set up file path-------
path_read_in <- '/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation'
#path_save <- '/gpfs/home/acad/ulg-mast/yuhaolin/rtest'
path_save <- '/gpfs/scratch/acad/bsmfc/Haolin_temp'

output_path_name<-'/outputs_population_model'
path_postfix<-'_5years_1990_spawning_678_2D_interpolation_daily_PTM_EL_update_DD'

current_test_name<-paste0('/SI_',number_individual,path_postfix) # this should be changed as different purposes
dir.create(paste0(path_save,output_path_name,current_test_name)) # if there is already has a same-name file folder there, it will do nothing
setwd(dir = path_read_in)

# set up initial number of each age class----------------
# run the first time to get the equilibrium values
# this could be set up by stock assessment report
# SSB in 2010: 500,000 t
# SSB in 1990: 200,000 t
# So, in 1990, the initial worth = 2/5 * 2010

pop_num_age_1<-3.8e+09
pop_num_age_2<-1.28e+09
pop_num_age_3<-4.25e+08
pop_num_age_4<-1.4e+08

#############################################
####
###     read in basic functions (DEB, movement, mortality) and basic information of black sea
####
########################################

# read in the final basic point information of Black Sea
blacksea_points <- readRDS('blacksea_points_PTM.rds')
blacksea_points_uo<-readRDS('blacksea_points_uo.rds')
blacksea_points_vo<-readRDS('blacksea_points_vo.rds')
basic.grid <- readRDS('basic.grid.rds')
cell.number <- readRDS('cell_number.rds')
empty.cell.number <- readRDS('empty_cell_number.rds')
max_depth_level<-readRDS('max_depth_level.rds')
all_depth_levels_t<-readRDS('depth_levels_from_BAMHBI.rds')
all_depth_levels_t<-as.vector(all_depth_levels_t)
all_depth_levels_w<-readRDS('depthw.rds')
all_depth_levels_w<-as.vector(all_depth_levels_w)

# set up parameters for the whole simulation------
source('set_up_all_parameters.R')
source('extract_data_in_2D_and_3D_HPC.R')
source('functional_response_HPC.R')
source('functions_to_transfer_ordinate_with_cell_grid.R')
source('interpolation_functions_2D_and_3D_HPC.R')
source('get_Temp_food_function_2D_HPC.R')
source('get_Temp_food_function_3D_HPC.R')
source('determine_spawn_day_end_31_Aug.R') 
source('DEB_function_HPC_spawning_678.R')
source('Kinesis_movement_function_HPC.R')
source('PTM_functions_2D_HPC.R')
source('PTM_functions_3D_HPC.R')
source('running_func_age_0_4_2D_seperate.R')
source('running_func_age_0_4_3D_seperate.R')
source('running_func_all_ages_2D_combine.R')
source('running_func_all_ages_3D_combine.R')
source('mortality_functions_2D.R')
source('mortality_functions_3D.R')

#############################################
####
###     allocate related packages, basic functions, and variables to each CPU
####
########################################

clusterEvalQ(cl, {
  library(lubridate)  
  library(sf)
  library(tidyverse)
  library(dplyr)
  library(raster)
  library(terra)
  library(sp)
  library(fst)
  
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/set_up_all_parameters.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/extract_data_in_2D_and_3D_HPC.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/functional_response_HPC.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/functions_to_transfer_ordinate_with_cell_grid.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/interpolation_functions_2D_and_3D_HPC.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/get_Temp_food_function_2D_HPC.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/get_Temp_food_function_3D_HPC.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/determine_spawn_day_end_31_Aug.R') 
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/DEB_function_HPC_spawning_678.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/Kinesis_movement_function_HPC.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/PTM_functions_2D_HPC.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/PTM_functions_3D_HPC.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/running_func_age_0_4_2D_seperate.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/running_func_age_0_4_3D_seperate.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/running_func_all_ages_2D_combine.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/running_func_all_ages_3D_combine.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/mortality_functions_2D.R')
  source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/Functions_variables_for_formal_simulation/mortality_functions_3D.R')
})

clusterExport(cl, c('split_indices',
                    'input.func',
                    'input.dimension',
                    'input.data.location',
                    'input.stoke.drift',
                    'timestep.set',
                    'parms.total',
                    'rho_particle',
                    'individual_per_core',
                    'blacksea_points',
                    'blacksea_points_uo',
                    'blacksea_points_vo',
                    'basic.grid',
                    'cell.number',
                    'empty.cell.number',
                    'max_depth_level',
                    'all_depth_levels_t',
                    'all_depth_levels_w',
                    'life_stages_in_PTM',
                    'life_stages_in_Kinesis',
                    'life_stages_in_DEB',
                    'life_stages_in_Stokes_drift'))

#############################################
####
###     define initial data frames, outputs, and set up parallel  cores
####
########################################

# 1. define initial individuals and initial state variables for each age class------------

# different initial points - 4800 individuals for each age class (Age-class 1 ~ 4)
# Age_class 0 should wait the new born individuals

set.seed(101)
potential_points_num_age_1 <- sample(1:nrow(blacksea_points),
                                     size = number_individual,
                                     replace = F)

set.seed(102)
potential_points_num_age_2 <- sample(1:nrow(blacksea_points),
                                     size = number_individual,
                                     replace = F)

set.seed(103)
potential_points_num_age_3 <- sample(1:nrow(blacksea_points),
                                     size = number_individual,
                                     replace = F)

set.seed(104)
potential_points_num_age_4 <- sample(1:nrow(blacksea_points),
                                     size = number_individual,
                                     replace = F)

dif_points_age_1 <- data.frame(Lon = blacksea_points$Longitude[potential_points_num_age_1], Lat = blacksea_points$Latitude[potential_points_num_age_1])
dif_points_age_2 <- data.frame(Lon = blacksea_points$Longitude[potential_points_num_age_2], Lat = blacksea_points$Latitude[potential_points_num_age_2])
dif_points_age_3 <- data.frame(Lon = blacksea_points$Longitude[potential_points_num_age_3], Lat = blacksea_points$Latitude[potential_points_num_age_3])
dif_points_age_4 <- data.frame(Lon = blacksea_points$Longitude[potential_points_num_age_4], Lat = blacksea_points$Latitude[potential_points_num_age_4])

cell_number_initial_age_1 <- NULL
for (i in 1:number_individual) {
  aaa <- trans_lonlat_to_cellnumber(lon = dif_points_age_1[i, 1], lat = dif_points_age_1[i, 2])
  cell_number_initial_age_1 <- c(cell_number_initial_age_1, aaa)
}
cell_number_initial_age_2 <- NULL
for (i in 1:number_individual) {
  aaa <- trans_lonlat_to_cellnumber(lon = dif_points_age_2[i, 1], lat = dif_points_age_2[i, 2])
  cell_number_initial_age_2 <- c(cell_number_initial_age_2, aaa)
}
cell_number_initial_age_3 <- NULL
for (i in 1:number_individual) {
  aaa <- trans_lonlat_to_cellnumber(lon = dif_points_age_3[i, 1], lat = dif_points_age_3[i, 2])
  cell_number_initial_age_3 <- c(cell_number_initial_age_3, aaa)
}
cell_number_initial_age_4 <- NULL
for (i in 1:number_individual) {
  aaa <- trans_lonlat_to_cellnumber(lon = dif_points_age_4[i, 1], lat = dif_points_age_4[i, 2])
  cell_number_initial_age_4 <- c(cell_number_initial_age_4, aaa)
}

# determine dist to shoreline
dis_to_coast_initial_age_1<-blacksea_points[match(cell_number_initial_age_1,blacksea_points$Cell_num),'dis_to_coast']
dis_to_coast_initial_age_2<-blacksea_points[match(cell_number_initial_age_2,blacksea_points$Cell_num),'dis_to_coast']
dis_to_coast_initial_age_3<-blacksea_points[match(cell_number_initial_age_3,blacksea_points$Cell_num),'dis_to_coast']
dis_to_coast_initial_age_4<-blacksea_points[match(cell_number_initial_age_4,blacksea_points$Cell_num),'dis_to_coast']

# determine initial region
region_initial_age_1<-blacksea_points[match(cell_number_initial_age_1,blacksea_points$Cell_num),'EEZs_adj']
region_initial_age_2<-blacksea_points[match(cell_number_initial_age_2,blacksea_points$Cell_num),'EEZs_adj']
region_initial_age_3<-blacksea_points[match(cell_number_initial_age_3,blacksea_points$Cell_num),'EEZs_adj']
region_initial_age_4<-blacksea_points[match(cell_number_initial_age_4,blacksea_points$Cell_num),'EEZs_adj']

# determine initial depth
max_depth_initial_age_1<- max_depth_level[match(cell_number_initial_age_1,max_depth_level$Cell_num),'max_depth_level']
max_depth_initial_age_2<- max_depth_level[match(cell_number_initial_age_2,max_depth_level$Cell_num),'max_depth_level']
max_depth_initial_age_3<- max_depth_level[match(cell_number_initial_age_3,max_depth_level$Cell_num),'max_depth_level']
max_depth_initial_age_4<- max_depth_level[match(cell_number_initial_age_4,max_depth_level$Cell_num),'max_depth_level']

depth_level_initial_age_1<-NULL
depth_initial_age_1<-NULL
set.seed(101)
for (i in 1:number_individual) {
  limited_depth_level<-max_depth_initial_age_1[i]
  if(max_depth_initial_age_1[i]>30){ # initial depth level <=30, so depth <= 100m
    limited_depth_level<-30
  }
  aaa <- sample(1:limited_depth_level,size = 1)
  depth_level_initial_age_1 <- c(depth_level_initial_age_1, aaa)
  bbb <- all_depth_levels_t[aaa]
  depth_initial_age_1<-c(depth_initial_age_1,bbb)
}
depth_level_initial_age_2<-NULL
depth_initial_age_2<-NULL
set.seed(102)
for (i in 1:number_individual) {
  limited_depth_level<-max_depth_initial_age_2[i]
  if(max_depth_initial_age_2[i]>30){ # initial depth level <=30, so depth <= 100m
    limited_depth_level<-30
  }
  aaa <- sample(1:limited_depth_level,size = 1)
  depth_level_initial_age_2 <- c(depth_level_initial_age_2, aaa)
  bbb <- all_depth_levels_t[aaa]
  depth_initial_age_2<-c(depth_initial_age_2,bbb)
}
depth_level_initial_age_3<-NULL
depth_initial_age_3<-NULL
set.seed(103)
for (i in 1:number_individual) {
  limited_depth_level<-max_depth_initial_age_3[i]
  if(max_depth_initial_age_3[i]>30){ # initial depth level <=30, so depth <= 100m
    limited_depth_level<-30
  }
  aaa <- sample(1:limited_depth_level,size = 1)
  depth_level_initial_age_3 <- c(depth_level_initial_age_3, aaa)
  bbb <- all_depth_levels_t[aaa]
  depth_initial_age_3<-c(depth_initial_age_3,bbb)
}
depth_level_initial_age_4<-NULL
depth_initial_age_4<-NULL
set.seed(104)
for (i in 1:number_individual) {
  limited_depth_level<-max_depth_initial_age_4[i]
  if(max_depth_initial_age_4[i]>30){ # initial depth level <=30, so depth <= 100m
    limited_depth_level<-30
  }
  aaa <- sample(1:limited_depth_level,size = 1)
  depth_level_initial_age_4 <- c(depth_level_initial_age_4, aaa)
  bbb <- all_depth_levels_t[aaa]
  depth_initial_age_4<-c(depth_initial_age_4,bbb)
}

# define initial date variable
date_vector = format(start_date, "%Y-%m-%d %H:%M:%S")
initial.Date <- c(rep(date_vector, number_individual))

# read in initial state variable for different age-class
state_var_initial_movement_age_1<-readRDS('state_var_move_age_1.rds')
state_var_initial_DEB_age_1<-readRDS('state_var_DEB_age_1.rds')

state_var_initial_movement_age_2<-readRDS('state_var_move_age_2.rds')
state_var_initial_DEB_age_2<-readRDS('state_var_DEB_age_2.rds')

state_var_initial_movement_age_3<-readRDS('state_var_move_age_3.rds')
state_var_initial_DEB_age_3<-readRDS('state_var_DEB_age_3.rds')

state_var_initial_movement_age_4<-readRDS('state_var_move_age_4.rds')
state_var_initial_DEB_age_4<-readRDS('state_var_DEB_age_4.rds')

# set up initial Temp and Food as NA for initial condition
initial_Temp_Food_age_1<-data.frame(Temp = NA_real_,
                                    Temp_spawn = NA_real_,
                                    Z_micro_zoo = NA_real_,
                                    Z_meso_zoo = NA_real_,
                                    f_micro_zoo = NA_real_,
                                    f_meso_zoo = NA_real_,
                                    f = NA_real_)
initial_Temp_Food_age_2<-data.frame(Temp = NA_real_,
                                    Temp_spawn = NA_real_,
                                    Z_micro_zoo = NA_real_,
                                    Z_meso_zoo = NA_real_,
                                    f_micro_zoo = NA_real_,
                                    f_meso_zoo = NA_real_,
                                    f = NA_real_)
initial_Temp_Food_age_3<-data.frame(Temp = NA_real_,
                                    Temp_spawn = NA_real_,
                                    Z_micro_zoo = NA_real_,
                                    Z_meso_zoo = NA_real_,
                                    f_micro_zoo = NA_real_,
                                    f_meso_zoo = NA_real_,
                                    f = NA_real_)
initial_Temp_Food_age_4<-data.frame(Temp = NA_real_,
                                    Temp_spawn = NA_real_,
                                    Z_micro_zoo = NA_real_,
                                    Z_meso_zoo = NA_real_,
                                    f_micro_zoo = NA_real_,
                                    f_meso_zoo = NA_real_,
                                    f = NA_real_)

# initialize PTM
initial_PTM_age_1<-data.frame(PTM_X = NA_real_,
                              PTM_Y = NA_real_,
                              PTM_Z = NA_real_)
initial_PTM_age_2<-data.frame(PTM_X = NA_real_,
                              PTM_Y = NA_real_,
                              PTM_Z = NA_real_)
initial_PTM_age_3<-data.frame(PTM_X = NA_real_,
                              PTM_Y = NA_real_,
                              PTM_Z = NA_real_)
initial_PTM_age_4<-data.frame(PTM_X = NA_real_,
                              PTM_Y = NA_real_,
                              PTM_Z = NA_real_)

# combine each other to generate initial condition for each age class

combined_initial_location_age_1<-cbind(No_ind=1:number_individual,
                                       Date=initial.Date,
                                       dif_points_age_1,
                                       Cell_num=cell_number_initial_age_1,
                                       depth=depth_initial_age_1,
                                       depth_level=depth_level_initial_age_1,
                                       dis_to_coast=dis_to_coast_initial_age_1,
                                       EEZs_region=region_initial_age_1,
                                       initial_Temp_Food_age_1,
                                       state_var_initial_DEB_age_1,
                                       state_var_initial_movement_age_1,
                                       initial_PTM_age_1,
                                       running_status='alive',
                                       death_causes=NA_character_,
                                       running_pattern='spin_up',
                                       running_year=1,
                                       age=1,
                                       worth=pop_num_age_1/number_individual,
                                       density=pop_num_age_1/number_individual/(2780*(2780*cos(dif_points_age_1$Lat*pi/180)))/100,
                                       norm_density=NA_real_, # all adults
                                       natural_mortality=0,
                                       ZFAC=0,
                                       final_mortality=0,
                                       jump_egg_days=1,  # jump developing days from eggs to yolksac larvae, 3 days in total
                                       jump_yolksac_days=1,  # jump developing days from yolksac larvae to early_larvae, 3 days in total
                                       survival_duration=0  # survival duration used to delete slowly grow-up individuals
) 

combined_initial_location_age_2<-cbind(No_ind=1:number_individual,
                                       Date=initial.Date,
                                       dif_points_age_2,
                                       Cell_num=cell_number_initial_age_2,
                                       depth=depth_initial_age_2,
                                       depth_level=depth_level_initial_age_2,
                                       dis_to_coast=dis_to_coast_initial_age_2,
                                       EEZs_region=region_initial_age_2,
                                       initial_Temp_Food_age_2,
                                       state_var_initial_DEB_age_2,
                                       state_var_initial_movement_age_2,
                                       initial_PTM_age_2,
                                       running_status='alive',
                                       death_causes=NA_character_,
                                       running_pattern='spin_up',
                                       running_year=1,
                                       age=2,
                                       worth=pop_num_age_2/number_individual,
                                       density=pop_num_age_2/number_individual/(2780*(2780*cos(dif_points_age_2$Lat*pi/180)))/100,
                                       norm_density=NA_real_, # all adults
                                       natural_mortality=0,
                                       ZFAC=0,
                                       final_mortality=0,
                                       jump_egg_days=1,  # jump developing days from eggs to yolksac larvae, 3 days in total
                                       jump_yolksac_days=1,  # jump developing days from yolksac larvae to early_larvae, 3 days in total
                                       survival_duration=0)

combined_initial_location_age_3<-cbind(No_ind=1:number_individual,
                                       Date=initial.Date,
                                       dif_points_age_3,
                                       Cell_num=cell_number_initial_age_3,
                                       depth=depth_initial_age_3,
                                       depth_level=depth_level_initial_age_3,
                                       dis_to_coast=dis_to_coast_initial_age_3,
                                       EEZs_region=region_initial_age_3,
                                       initial_Temp_Food_age_3,
                                       state_var_initial_DEB_age_3,
                                       state_var_initial_movement_age_3,
                                       initial_PTM_age_3,
                                       running_status='alive',
                                       death_causes=NA_character_,
                                       running_pattern='spin_up',
                                       running_year=1,
                                       age=3,
                                       worth=pop_num_age_3/number_individual,
                                       density=pop_num_age_3/number_individual/(2780*(2780*cos(dif_points_age_3$Lat*pi/180)))/100,
                                       norm_density=NA_real_, # all adults
                                       natural_mortality=0,
                                       ZFAC=0,
                                       final_mortality=0,
                                       jump_egg_days=1,  # jump developing days from eggs to yolksac larvae, 3 days in total
                                       jump_yolksac_days=1,  # jump developing days from yolksac larvae to early_larvae, 3 days in total
                                       survival_duration=0)

combined_initial_location_age_4<-cbind(No_ind=1:number_individual,
                                       Date=initial.Date,
                                       dif_points_age_4,
                                       Cell_num=cell_number_initial_age_4,
                                       depth=depth_initial_age_4,
                                       depth_level=depth_level_initial_age_4,
                                       dis_to_coast=dis_to_coast_initial_age_4,
                                       EEZs_region=region_initial_age_4,
                                       initial_Temp_Food_age_4,
                                       state_var_initial_DEB_age_4,
                                       state_var_initial_movement_age_4,
                                       initial_PTM_age_4,
                                       running_status='alive',
                                       death_causes=NA_character_,
                                       running_pattern='spin_up',
                                       running_year=1,
                                       age=4,
                                       worth=pop_num_age_4/number_individual,
                                       density=pop_num_age_4/number_individual/(2780*(2780*cos(dif_points_age_4$Lat*pi/180)))/100,
                                       norm_density=NA_real_, # all adults
                                       natural_mortality=0,
                                       ZFAC=0,
                                       final_mortality=0,
                                       jump_egg_days=1,  # jump developing days from eggs to yolksac larvae, 3 days in total
                                       jump_yolksac_days=1,  # jump developing days from yolksac larvae to early_larvae, 3 days in total
                                       survival_duration=0)

# Age-0 is special because they are waiting for the initialization
# transfer every element into 0, and only set up life_stage as none, and running_status as silence

combined_initial_location_age_0<-combined_initial_location_age_1
combined_initial_location_age_0[,c(3:ncol(combined_initial_location_age_0))]<-0  # age was set up 0 here
combined_initial_location_age_0[,'life_stage']<-NA_character_
combined_initial_location_age_0[,'EEZs_region']<-NA_character_
combined_initial_location_age_0[,'running_status']<-'silence'
combined_initial_location_age_0[,'death_causes']<-NA_character_
combined_initial_location_age_0[,'running_pattern']<-'spin_up'
combined_initial_location_age_0[,'running_year']<-1

# 2. define dynamic (spawning) and final outputs (Global monitor and Age Lists) ----------
# dynamic spawning data frame 
spawning_dynamic_monitor<-data.frame(Date=NA_character_,
                                     Ukraine_1 = 0,
                                     Ukraine_2 = 0,
                                     Ukraine_3 = 0,
                                     Romania = 0,
                                     Bulgaria = 0,
                                     Turkey_1 = 0,
                                     Turkey_2 = 0,
                                     Turkey_3 = 0,
                                     Georgia = 0,
                                     Russia = 0,
                                     running_pattern = 'spin_up',
                                     running_year = 1,
                                     number_of_spawner_x_worth = 0,
                                     total_eggs_x_worth = 0)

# global data frame to monitor status of SIs in each age class on each time step 
state_var_initial_combine<-rbind(combined_initial_location_age_0,
                                 combined_initial_location_age_1,
                                 combined_initial_location_age_2,
                                 combined_initial_location_age_3,
                                 combined_initial_location_age_4)

global_monitor_eachtimestep<-state_var_initial_combine %>%
  mutate(biomass = ifelse(running_status == 'alive', worth * Www, NA_real_)) %>%
  group_by(age)%>%
  summarise(Date=first(Date),
            Dist.mean=mean(ifelse(running_status == 'alive', Dist, NA_real_), na.rm = TRUE),
            depth.mean=mean(ifelse(running_status == 'alive', depth, NA_real_), na.rm = TRUE),
            dis.to.coast.mean=mean(ifelse(running_status == 'alive', dis_to_coast, NA_real_), na.rm = TRUE),
            Temp.mean=mean(ifelse(running_status == 'alive', Temp, NA_real_), na.rm = TRUE),
            Temp.spawn.mean=mean(ifelse(running_status == 'alive', Temp_spawn, NA_real_), na.rm = TRUE),
            Z.micro.zoo.mean=mean(ifelse(running_status == 'alive', Z_micro_zoo, NA_real_), na.rm = TRUE),
            Z.meso.zoo.mean=mean(ifelse(running_status == 'alive', Z_meso_zoo, NA_real_), na.rm = TRUE),
            f.micro.zoo.mean=mean(ifelse(running_status == 'alive', f_micro_zoo, NA_real_), na.rm = TRUE),
            f.meso.zoo.mean=mean(ifelse(running_status == 'alive', f_meso_zoo, NA_real_), na.rm = TRUE),
            Food.f.mean=mean(ifelse(running_status == 'alive', f, NA_real_), na.rm = TRUE),
            TL.mean=mean(ifelse(running_status == 'alive', TL, NA_real_), na.rm = TRUE),
            Www.mean=mean(ifelse(running_status == 'alive', Www, NA_real_), na.rm = TRUE),
            eggs.mean=mean(ifelse(running_status == 'alive', eggs, NA_real_), na.rm = TRUE),
            worth.mean=mean(ifelse(running_status == 'alive', worth, NA_real_), na.rm = TRUE),
            total.biomass=sum(biomass, na.rm = TRUE),
            total.worth=sum(ifelse(running_status == 'alive', worth, NA_real_), na.rm = TRUE),
            density.mean=mean(ifelse(running_status == 'alive', density, NA_real_), na.rm = TRUE),
            norm.density.mean=mean(ifelse(running_status == 'alive', norm_density, NA_real_), na.rm = TRUE),
            natural.mortality.mean=mean(ifelse(running_status == 'alive', natural_mortality, NA_real_), na.rm = TRUE),
            ZFAC.mean=mean(ifelse(running_status == 'alive', ZFAC, NA_real_), na.rm = TRUE),
            final.mortality.mean=mean(ifelse(running_status == 'alive', final_mortality, NA_real_), na.rm = TRUE),
            total.eggs.worth=sum(ifelse(running_status == 'alive', eggs, NA_real_), na.rm = TRUE),
            running.pattern=running_pattern[1],
            running.year=running_year[1],
            number.adults=sum(ifelse(life_stage == 'adult', worth, NA_real_), na.rm = TRUE),
            Num_alive_SI=sum(running_status == 'alive', na.rm = TRUE),
            remainder_rate_SI=(sum(running_status == 'alive', na.rm = TRUE)/number_individual),
            remainder_rate_SI_worth=(sum(ifelse(running_status == "alive", worth, 0), na.rm = TRUE)/sum(worth, na.rm = TRUE)),
            Num_Killed_by_stunted=sum(running_status == 'Killed_by_stunted', na.rm = TRUE),
            Num_Killed_by_low_worth=sum(running_status == 'Killed_by_low_worth', na.rm = TRUE))%>%
  as.data.frame()

global_monitor_eachday<-state_var_initial_combine %>%
  mutate(biomass = ifelse(running_status == 'alive', worth * Www, NA_real_)) %>%
  group_by(age)%>%
  summarise(Date=first(Date),
            Dist.mean=mean(ifelse(running_status == 'alive', Dist, NA_real_), na.rm = TRUE),
            depth.mean=mean(ifelse(running_status == 'alive', depth, NA_real_), na.rm = TRUE),
            dis.to.coast.mean=mean(ifelse(running_status == 'alive', dis_to_coast, NA_real_), na.rm = TRUE),
            Temp.mean=mean(ifelse(running_status == 'alive', Temp, NA_real_), na.rm = TRUE),
            Temp.spawn.mean=mean(ifelse(running_status == 'alive', Temp_spawn, NA_real_), na.rm = TRUE),
            Z.micro.zoo.mean=mean(ifelse(running_status == 'alive', Z_micro_zoo, NA_real_), na.rm = TRUE),
            Z.meso.zoo.mean=mean(ifelse(running_status == 'alive', Z_meso_zoo, NA_real_), na.rm = TRUE),
            f.micro.zoo.mean=mean(ifelse(running_status == 'alive', f_micro_zoo, NA_real_), na.rm = TRUE),
            f.meso.zoo.mean=mean(ifelse(running_status == 'alive', f_meso_zoo, NA_real_), na.rm = TRUE),
            Food.f.mean=mean(ifelse(running_status == 'alive', f, NA_real_), na.rm = TRUE),
            TL.mean=mean(ifelse(running_status == 'alive', TL, NA_real_), na.rm = TRUE),
            Www.mean=mean(ifelse(running_status == 'alive', Www, NA_real_), na.rm = TRUE),
            eggs.mean=mean(ifelse(running_status == 'alive', eggs, NA_real_), na.rm = TRUE),
            worth.mean=mean(ifelse(running_status == 'alive', worth, NA_real_), na.rm = TRUE),
            total.biomass=sum(biomass, na.rm = TRUE),
            total.worth=sum(ifelse(running_status == 'alive', worth, NA_real_), na.rm = TRUE),
            density.mean=mean(ifelse(running_status == 'alive', density, NA_real_), na.rm = TRUE),
            norm.density.mean=mean(ifelse(running_status == 'alive', norm_density, NA_real_), na.rm = TRUE),
            natural.mortality.mean=mean(ifelse(running_status == 'alive', natural_mortality, NA_real_), na.rm = TRUE),
            ZFAC.mean=mean(ifelse(running_status == 'alive', ZFAC, NA_real_), na.rm = TRUE),
            final.mortality.mean=mean(ifelse(running_status == 'alive', final_mortality, NA_real_), na.rm = TRUE),
            total.eggs.worth=sum(ifelse(running_status == 'alive', eggs, NA_real_), na.rm = TRUE),
            running.pattern=running_pattern[1],
            running.year=running_year[1],
            number.adults=sum(ifelse(life_stage == 'adult', worth, NA_real_), na.rm = TRUE),
            Num_alive_SI=sum(running_status == 'alive', na.rm = TRUE),
            remainder_rate_SI=(sum(running_status == 'alive', na.rm = TRUE)/number_individual),
            remainder_rate_SI_worth=(sum(ifelse(running_status == "alive", worth, 0), na.rm = TRUE)/sum(worth, na.rm = TRUE)),
            Num_Killed_by_stunted=sum(running_status == 'Killed_by_stunted', na.rm = TRUE),
            Num_Killed_by_low_worth=sum(running_status == 'Killed_by_low_worth', na.rm = TRUE))%>%
  as.data.frame()

# 4. define a list to record ZFAC for each day
mortality_age_monitor<-list()

# 5. define a data frame to accumulate record the duration of life stage, spawning batches, mean distance with shoreline, and 
#    experienced days in each regions at each year
annual_record_SI<-data.frame()

# define final results with the big list
results_total_age_0 <- list(combined_initial_location_age_0)
results_total_age_1 <- list(combined_initial_location_age_1)
results_total_age_2 <- list(combined_initial_location_age_2)
results_total_age_3 <- list(combined_initial_location_age_3)
results_total_age_4 <- list(combined_initial_location_age_4)


# 3. define other output data frames --------------------
# for seasonal spatial map of abundance (eggs) and biomass (early-, late-larvae, juvenile, and adult) 
egg_seasonal_matrix_abundance<-matrix(0, nrow = running_days, ncol = length(cell.number))
yolksac_larvae_seasonal_matrix_abundance<-matrix(0, nrow = running_days, ncol = length(cell.number))

early_larvae_seasonal_matrix_biomass<-matrix(0, nrow = running_days, ncol = length(cell.number))
late_larvae_seasonal_matrix_biomass<-matrix(0, nrow = running_days, ncol = length(cell.number))
juvenile_seasonal_matrix_biomass<-matrix(0, nrow = running_days, ncol = length(cell.number))
adult_seasonal_matrix_biomass<-matrix(0, nrow = running_days, ncol = length(cell.number))

early_larvae_seasonal_matrix_abundance<-matrix(0, nrow = running_days, ncol = length(cell.number))
late_larvae_seasonal_matrix_abundance<-matrix(0, nrow = running_days, ncol = length(cell.number))
juvenile_seasonal_matrix_abundance<-matrix(0, nrow = running_days, ncol = length(cell.number))
adult_seasonal_matrix_abundance<-matrix(0, nrow = running_days, ncol = length(cell.number))

total_seasonal_matrix_abundance<-matrix(0, nrow = running_days, ncol = length(cell.number))
total_seasonal_matrix_biomass<-matrix(0, nrow = running_days, ncol = length(cell.number))

# for S-R relationship
SR_relationship<-data.frame()

# annual weight and length for calibration
annual_weight_age<-data.frame()
annual_length_age<-data.frame()
annual_Lon_Lat_age<-data.frame()

# summary_eggs_per_batch_per_age_per_year
eggs_batch_num_age_0<-matrix(0,nrow = running_days,ncol = number_individual)
eggs_batch_num_age_1<-matrix(0,nrow = running_days,ncol = number_individual)
eggs_batch_num_age_2<-matrix(0,nrow = running_days,ncol = number_individual)
eggs_batch_num_age_3<-matrix(0,nrow = running_days,ncol = number_individual)
eggs_batch_num_age_4<-matrix(0,nrow = running_days,ncol = number_individual)

eggs_total_num_age_0<-matrix(0,nrow = running_days,ncol = number_individual)
eggs_total_num_age_1<-matrix(0,nrow = running_days,ncol = number_individual)
eggs_total_num_age_2<-matrix(0,nrow = running_days,ncol = number_individual)
eggs_total_num_age_3<-matrix(0,nrow = running_days,ncol = number_individual)
eggs_total_num_age_4<-matrix(0,nrow = running_days,ncol = number_individual)

# randomly find each 3 typical individuals born in June, July, and August
information_selected_individuals<-data.frame()
target_ind_row_num<-c(c(51,
                        501,
                        1001),                                        # born in June
                      c((number_SI_for_each_region*10*(30)+51),
                        (number_SI_for_each_region*10*(30)+501),
                        (number_SI_for_each_region*10*(30)+1001)),    # born in July
                      c((number_SI_for_each_region*10*(30+31)+51),
                        (number_SI_for_each_region*10*(30+31)+501),
                        (number_SI_for_each_region*10*(30+31)+1001))) # born in Aug

# orther monitor and outputs
all_potential_spawn_day<-data.frame()
spawning_eachday_output<-data.frame()
day_counter<-0

# specific date that do not have values for environmental variables (NEMO, BAMHBI)
unusable.date<-c(as.Date('1999-02-17'),as.Date('2000-01-21'),as.Date('2000-06-02'),as.Date('2001-08-02'),
                 as.Date('2003-06-16'),as.Date('2004-11-24'),as.Date('2006-11-16'),as.Date('2007-04-07'),
                 as.Date('2007-10-19'),as.Date('2009-04-18'),as.Date('2010-01-12'),as.Date('2014-10-04'),
                 as.Date('2015-04-05'),as.Date('2016-12-25'),as.Date('2017-02-23'),as.Date('2018-09-15'),
                 as.Date('2020-01-12'),as.Date('2021-07-14'),as.Date('2021-07-27'),as.Date('2021-08-29'),
                 as.Date('2021-09-06'))

# define data frame of mortality for all days
overlap.occur.all.days<-data.frame()

#############################################
####
###     start simulation for each age class
####
########################################

for (i in 1:length(years)) {
  # >>>>>>> start loop for years
  # determine earliest and last possible day for spawn
  spawn_day<-determine.spawn.day(current.year = years[i])
  all_potential_spawn_day<-rbind(all_potential_spawn_day,spawn_day)
  
  # determine annual record of SI
  accumulate_daily_record_SI<-data.frame(Year = years[i],
                                         age = rep(0:4,each = number_individual),
                                         No_ind = rep(1:number_individual,times = 5),
                                         egg_days = 0,
                                         yolksac_larvae_days = 0,
                                         early_larvae_days = 0,
                                         late_larvae_days = 0,
                                         juvenile_days = 0,
                                         adult_days = 0,
                                         batch_num_no_worth = 0,
                                         batch_num_with_worth = 0,
                                         eggs_total_num_no_worth = 0,
                                         eggs_total_num_with_worth = 0,
                                         sum_dis_to_coast = 0,
                                         days_in_Ukraine_1 = 0,
                                         days_in_Ukraine_2 = 0,
                                         days_in_Ukraine_3 = 0,
                                         days_in_Romania = 0,
                                         days_in_Bulgaria = 0,
                                         days_in_Turkey_1 = 0,
                                         days_in_Turkey_2 = 0,
                                         days_in_Turkey_3 = 0,
                                         days_in_Georgia = 0,
                                         days_in_Russia = 0)
  
  # find how many days in each year
  # if initial days is not from the first day of one year, so we need this "if" argument
  if (length(years) > 1) {
    # more than one year
    if (i == 1) {
      # the first year
      days_in_year <- as.numeric(difftime(as.Date(paste0(years[i], '-12-31')), as.Date(start_date), units = 'days')) + 1
    } else if (i == length(years)) {
      # the last year
      days_in_year <- as.numeric(difftime(as.Date(end_date), as.Date(paste0(years[i], '-01-01')), units = 'days')) + 1
    } else{
      # the normal year (i.e., not the first or the last year)
      days_in_year <- ifelse(leap_year(as.numeric(years[i])), 366, 365)
    }
  } else{
    # only one year
    days_in_year <- as.numeric(difftime(as.Date(end_date), as.Date(start_date), units = 'days')) + 1
  }
  
  for (j in 1:days_in_year) {
    # >>>>>>>> start loop for days
    
    # read in the latest state variable as input for each day ==========================
    if(i>=2 & j == 1){ # update age-class, age, and killed age-4 class on the first day of each year
      input.state.var.age.0 <- combined_initial_location_age_0
      input.state.var.age.1 <- results_total_age_0[[length(results_total_age_0)]]
      input.state.var.age.2 <- results_total_age_1[[length(results_total_age_1)]]
      input.state.var.age.3 <- results_total_age_2[[length(results_total_age_2)]]
      input.state.var.age.4 <- results_total_age_3[[length(results_total_age_3)]]
      
      input.state.var.age.1$age<-input.state.var.age.1$age+1
      input.state.var.age.2$age<-input.state.var.age.2$age+1
      input.state.var.age.3$age<-input.state.var.age.3$age+1
      input.state.var.age.4$age<-input.state.var.age.4$age+1
      
      input.state.var.age.0[,'running_year']<-i
      input.state.var.age.1[,'running_year']<-i
      input.state.var.age.2[,'running_year']<-i
      input.state.var.age.3[,'running_year']<-i
      input.state.var.age.4[,'running_year']<-i
      
    }else{
      input.state.var.age.0 <- results_total_age_0[[length(results_total_age_0)]]
      input.state.var.age.1 <- results_total_age_1[[length(results_total_age_1)]]
      input.state.var.age.2 <- results_total_age_2[[length(results_total_age_2)]]
      input.state.var.age.3 <- results_total_age_3[[length(results_total_age_3)]]
      input.state.var.age.4 <- results_total_age_4[[length(results_total_age_4)]]
    }
    
    # Formal simulation: change running_pattern into 'formal_simulation'
    if(i>=6){ 
      input.state.var.age.0[,'running_pattern']<-'formal_simulation'
      input.state.var.age.1[,'running_pattern']<-'formal_simulation'
      input.state.var.age.2[,'running_pattern']<-'formal_simulation'
      input.state.var.age.3[,'running_pattern']<-'formal_simulation'
      input.state.var.age.4[,'running_pattern']<-'formal_simulation'
    }
    
    # update current date of the day ===========================
    if(i<=5){ # during the spin-up ### if only test one same year for several times, this needs to be changed 
      
      current.Date.in.day <- start_date + days((j-1)%%spinup_days_in_year)
      
    }else{ # in the formal simulation
      # current date in formal simulation has different logic with spin-up
      day_counter <- day_counter+1
      current.Date.in.day <- start_date + days(day_counter-1)
    }
    
    # read in potential environmental data we need ===========================
    # some dates do not have values for 2D and 3D, so we only used the data from the previous day
    if(current.Date.in.day %in% unusable.date){
      current.Date.in.day <- current.Date.in.day - days(1)
    }
    
    if(input.dimension == 1){ # 2D
      
      if(any(input.state.var.age.0$running_status == 'alive') & any(input.state.var.age.0$life_stage %in% life_stages_in_PTM)){ # during the spawning season that we need PTM data frames
        
        All.Temp.grid.current.in.day<-extract_current_Temp_2D(updated_date = current.Date.in.day)
        All.Food.grid.current.in.day<-extract_current_Food_2D(updated_date = current.Date.in.day)
        # whatever daily or hourly, PTM still needs hourly
        All.uo.grid.current.in.day<-extract_current_uo_2D(updated_date = current.Date.in.day)
        All.vo.grid.current.in.day<-extract_current_vo_2D(updated_date = current.Date.in.day)
        All.uo.grid.next.in.day<-extract_next_uo_2D(updated_date = current.Date.in.day)
        All.vo.grid.next.in.day<-extract_next_vo_2D(updated_date = current.Date.in.day)
        
        if(input.stoke.drift == 1){ # considering stokes drift - wave effects
          All.sd.grid.current.in.day<-extract_current_sd_2D(updated_date = current.Date.in.day)
          All.sd.grid.next.in.day<-extract_next_sd_2D(updated_date = current.Date.in.day)
        }else{
          All.sd.grid.current.in.day<-NULL
          All.sd.grid.next.in.day<-NULL
        }
        
        if(timestep.set >= 1){ # hourly: we NEED temporal interpolation for current and NEXT day
          All.Temp.grid.next.in.day<-extract_next_Temp_2D(updated_date = current.Date.in.day)
          All.Food.grid.next.in.day<-extract_next_Food_2D(updated_date = current.Date.in.day)
        }else{
          All.Temp.grid.next.in.day<-NULL
          All.Food.grid.next.in.day<-NULL
        }
        
      }else{ # if this is not the spawning season, we don't need PTM data frame
        All.Temp.grid.current.in.day<-extract_current_Temp_2D(updated_date = current.Date.in.day)
        All.Food.grid.current.in.day<-extract_current_Food_2D(updated_date = current.Date.in.day)
        
        All.uo.grid.current.in.day<-NULL
        All.vo.grid.current.in.day<-NULL
        All.uo.grid.next.in.day<-NULL
        All.vo.grid.next.in.day<-NULL
        All.sd.grid.current.in.day<-NULL
        All.sd.grid.next.in.day<-NULL
        
        if(timestep.set >= 1){ # hourly: we NEED temporal interpolation for current and NEXT day
          All.Temp.grid.next.in.day<-extract_next_Temp_2D(updated_date = current.Date.in.day)
          All.Food.grid.next.in.day<-extract_next_Food_2D(updated_date = current.Date.in.day)
        }else{
          All.Temp.grid.next.in.day<-NULL
          All.Food.grid.next.in.day<-NULL
        }
      }
      
    }else{ # 3D
      
      if(any(input.state.var.age.0$running_status == 'alive') & any(input.state.var.age.0$life_stage %in% life_stages_in_PTM)){ # during the spawning season that we need PTM data frames
        
        All.Temp.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                         var_name = 'Temp')
        All.Temp.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                   var_name = 'Temp')   # PTM need hourly Temp for 3D
        All.MES.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                        var_name = 'MES')
        All.MIC.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                        var_name = 'MIC')
        # whatever daily or hourly, PTM still needs hourly
        All.uo.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                       var_name = 'uo')
        All.vo.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                       var_name = 'vo')
        All.wo.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                       var_name = 'wo')
        All.rho.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                        var_name = 'rho')
        All.Sal.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                        var_name = 'Sal')
        
        All.uo.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                 var_name = 'uo')
        All.vo.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                 var_name = 'vo')
        All.wo.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                 var_name = 'wo')
        All.rho.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                  var_name = 'rho')
        All.Sal.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                  var_name = 'Sal')
        
        if(input.stoke.drift == 1){ # considering stokes drift - wave effects
          All.sd.grid.current.in.day<-extract_current_sd_2D(updated_date = current.Date.in.day)
          All.sd.grid.next.in.day<-extract_next_sd_2D(updated_date = current.Date.in.day)
        }else{
          All.sd.grid.current.in.day<-NULL
          All.sd.grid.next.in.day<-NULL
        }
        
        if(timestep.set >= 1){ # hourly: we NEED temporal interpolation for current and NEXT day
          All.MES.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                    var_name = 'MES')
          All.MIC.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                    var_name = 'MIC')
        }else{
          All.MES.grid.next.in.day<-NULL
          All.MIC.grid.next.in.day<-NULL
        }
        
      }else{ # if this is not the spawning season, we don't need PTM data frame
        All.Temp.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                         var_name = 'Temp')
        All.MES.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                        var_name = 'MES')
        All.MIC.grid.current.in.day<-extract_current_3D(updated_date = current.Date.in.day,
                                                        var_name = 'MIC')
        
        All.uo.grid.current.in.day<-NULL
        All.vo.grid.current.in.day<-NULL
        All.uo.grid.next.in.day<-NULL
        All.vo.grid.next.in.day<-NULL
        All.sd.grid.current.in.day<-NULL
        All.sd.grid.next.in.day<-NULL
        All.wo.grid.current.in.day<-NULL
        All.wo.grid.next.in.day<-NULL
        All.rho.grid.current.in.day<-NULL
        All.rho.grid.next.in.day<-NULL
        All.Sal.grid.current.in.day<-NULL
        All.Sal.grid.next.in.day<-NULL
        
        if(timestep.set >= 1){ # hourly: we NEED temporal interpolation for current and NEXT day
          All.Temp.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                     var_name = 'Temp')
          All.MES.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                    var_name = 'MES')
          All.MIC.grid.next.in.day<-extract_next_3D(updated_date = current.Date.in.day,
                                                    var_name = 'MIC')
        }else{
          All.Temp.grid.next.in.day<-NULL
          All.MES.grid.next.in.day<-NULL
          All.MIC.grid.next.in.day<-NULL
        }
      } # End IF for non spawning season (don't need PTM) 
    } # End IF for 3D
    
    # define number of hours in each day ========================
    if (timestep.set >= 1) {
      # hourly or smaller time step
      hours_in_day <- 24
    } else{
      # daily
      hours_in_day <- 1
    }
    
    # allocate environmental data each day and do not need every timestep
    if (input.dimension == 1) { # 2D
      
      # export necessary variables to each CPU because those will be changed at different time point
      clusterExport(cl, varlist = c("All.Temp.grid.current.in.day",
                                    "All.Food.grid.current.in.day",
                                    "All.Temp.grid.next.in.day",
                                    "All.Food.grid.next.in.day",
                                    "All.uo.grid.current.in.day",
                                    "All.vo.grid.current.in.day",
                                    "All.uo.grid.next.in.day",
                                    "All.vo.grid.next.in.day",
                                    "All.sd.grid.current.in.day",
                                    "All.sd.grid.next.in.day"))
    }else{ # 3D
      
      # export necessary variables to each CPU because those will be changed at different time point
      clusterExport(cl, varlist = c("All.Temp.grid.current.in.day",
                                    "All.Temp.grid.next.in.day",
                                    "All.MIC.grid.current.in.day",
                                    "All.MIC.grid.next.in.day",
                                    "All.MES.grid.current.in.day",
                                    "All.MES.grid.next.in.day",
                                    "All.uo.grid.current.in.day",
                                    "All.uo.grid.next.in.day",
                                    "All.vo.grid.current.in.day",
                                    "All.vo.grid.next.in.day",
                                    "All.wo.grid.current.in.day",
                                    "All.wo.grid.next.in.day",
                                    "All.rho.grid.current.in.day",
                                    "All.rho.grid.next.in.day",
                                    "All.Sal.grid.current.in.day",
                                    "All.Sal.grid.next.in.day",
                                    "All.sd.grid.current.in.day",
                                    "All.sd.grid.next.in.day"))
    }
    
    spawning_eachhour<-data.frame()
    for (k in 1:hours_in_day) {
      # define timestep in 1 / 24 = daily, 1 = hourly, 2 = 30 min, 3 = 20 min, 4 = 15 min, 5 = 12 min ..., 60 = 1 min
      
      if(timestep.set==1){ # hourly
        input.state.var.age.0 <- results_total_age_0[[length(results_total_age_0)]]
        input.state.var.age.1 <- results_total_age_1[[length(results_total_age_1)]]
        input.state.var.age.2 <- results_total_age_2[[length(results_total_age_2)]]
        input.state.var.age.3 <- results_total_age_3[[length(results_total_age_3)]]
        input.state.var.age.4 <- results_total_age_4[[length(results_total_age_4)]]
      }
      
      spawning_eachtimestep<-data.frame()
      for (m in 1:timestep.set) {  # this step is used to provide more tinny time step like 30min or 1 min
        
        if(timestep.set>1){ # if time step in 30 min, 20 min, or 1min
          input.state.var.age.0 <- results_total_age_0[[length(results_total_age_0)]]
          input.state.var.age.1 <- results_total_age_1[[length(results_total_age_1)]]
          input.state.var.age.2 <- results_total_age_2[[length(results_total_age_2)]]
          input.state.var.age.3 <- results_total_age_3[[length(results_total_age_3)]]
          input.state.var.age.4 <- results_total_age_4[[length(results_total_age_4)]]
        }
        
        # define current time point/length ==============================
        timepoint <- length(results_total_age_1)-1
        
        # update the current time point for each time step
        
        if(i<=5){ # during the spin-up ### if only test one same year for several times, this needs to be changed 
          
          if (timestep.set > 1) {
            current.Date <- start_date + minutes(timepoint%%(spinup_days_in_year*24*(60 / timestep.set)))
          } else if (timestep.set == 1) {
            current.Date <- start_date + hours(timepoint%%(spinup_days_in_year*24*timestep.set))
          } else{
            current.Date <- start_date + days(timepoint%%spinup_days_in_year)
          }
          
        }else{ # in the formal simulation
          # current date in formal simulation has different logic with spin-up
          # the current date in formal is continuous adding 
          if (timestep.set > 1) {
            new_timepoint <- timepoint-5*spinup_days_in_year*24*timestep.set
            current.Date <- start_date + minutes(new_timepoint* (60 / timestep.set))
          } else if (timestep.set == 1) {
            new_timepoint<-timepoint-5*spinup_days_in_year*24
            current.Date <- start_date + hours(new_timepoint)
          } else{
            new_timepoint<-timepoint-5*spinup_days_in_year
            current.Date <- start_date + days(new_timepoint)
          }
        }
        
        date_current <- format(current.Date, "%Y-%m-%d %H:%M:%S")
        current_month <- as.numeric(format(current.Date, "%m"))
        
        # update date for each SI ============================================
        input.state.var.age.0[,'Date']<-date_current
        input.state.var.age.1[,'Date']<-date_current
        input.state.var.age.2[,'Date']<-date_current
        input.state.var.age.3[,'Date']<-date_current
        input.state.var.age.4[,'Date']<-date_current
        
        # construct a list to contain all input state variables ==================================
        input.state.var.list <- list(input.state.var.age.0,
                                     input.state.var.age.1,
                                     input.state.var.age.2,
                                     input.state.var.age.3,
                                     input.state.var.age.4)
        
        # redefine row number which has eggs and yolksac to allocate to each CPU====================
        # this is aimed to let each CPU run PTM faster
        # no need to wait other CPUs
        temporal_split_indices<-split_indices
        
        if(any(input.state.var.age.0$life_stage %in% life_stages_in_PTM)){
          
          potential_row_num<-which(input.state.var.age.0$life_stage %in% life_stages_in_PTM)
          
          # 1. delete duplicated elements from split_indices and potential_row_num 
          filtered_elements <- setdiff(1:number_individual, potential_row_num)  # maintain elements not in the potential_row_num
          
          # 2. re-allocate filtered_elements into the number of used CPU
          num_groups <- length(split_indices)  # number of CPU
          num_filtered <- length(filtered_elements)  # 4800 - length(potential_row_num)
          new_split_indices <- split(filtered_elements, rep(1:num_groups, length.out=num_filtered))
          
          # 3. re-allocate potential_row_num into the number of used CPU
          num_potential <- length(potential_row_num)
          group_assignments <- split(potential_row_num, rep(1:num_groups, length.out=num_potential))
          
          # **supplement length of group_assignments as empty ** to avoid less than the number of CPUs
          while (length(group_assignments) < num_groups) {
            group_assignments[[length(group_assignments) + 1]] <- integer(0)  # supplement empty subset list
          }
          
          # 4. combine the new-allocated split_indices and potential_row_num
          temporal_split_indices <- lapply(seq_len(num_groups), function(i) {
            c(new_split_indices[[i]], group_assignments[[i]])  # combine directly by index
          })
          
        }
        
        # running in parallel in one time step ===============================
        
        clusterExport(cl, varlist = c("input.state.var.list","current_month","temporal_split_indices"))  # this should be updated in each time step
        
        if (input.dimension == 1) { # 2D
          
          results <- parLapply(cl, 1:number_cores, function(num_cpu) {
            set.seed(num_cpu)
            running_func_all_ages_2D(batch.of.No.individual = temporal_split_indices[[num_cpu]],
                                     input.state.var.list = input.state.var.list,
                                     All.Temp.grid.current.in.day = All.Temp.grid.current.in.day,
                                     All.Food.grid.current.in.day = All.Food.grid.current.in.day,
                                     All.Temp.grid.next.in.day = All.Temp.grid.next.in.day,
                                     All.Food.grid.next.in.day = All.Food.grid.next.in.day,
                                     All.uo.grid.current.in.day = All.uo.grid.current.in.day,
                                     All.vo.grid.current.in.day = All.vo.grid.current.in.day,
                                     All.uo.grid.next.in.day = All.uo.grid.next.in.day,
                                     All.vo.grid.next.in.day = All.vo.grid.next.in.day,
                                     All.sd.grid.current.in.day = All.sd.grid.current.in.day,
                                     All.sd.grid.next.in.day = All.sd.grid.next.in.day,
                                     input.func = input.func,
                                     input.data.location = input.data.location,
                                     input.stoke.drift = input.stoke.drift,
                                     timestep = timestep.set,
                                     parms.total = parms.total,
                                     individual_per_core = individual_per_core,
                                     current_month = current_month)
          })
        }else{ # 3D
          
          results <- parLapply(cl, 1:number_cores, function(num_cpu) {
            set.seed(num_cpu)
            running_func_all_ages_3D(batch.of.No.individual = temporal_split_indices[[num_cpu]],
                                     input.state.var.list = input.state.var.list,
                                     All.Temp.grid.current.in.day = All.Temp.grid.current.in.day,
                                     All.Temp.grid.next.in.day = All.Temp.grid.next.in.day,
                                     All.MIC.grid.current.in.day = All.MIC.grid.current.in.day,
                                     All.MIC.grid.next.in.day = All.MIC.grid.next.in.day,
                                     All.MES.grid.current.in.day = All.MES.grid.current.in.day,
                                     All.MES.grid.next.in.day = All.MES.grid.next.in.day,
                                     All.uo.grid.current.in.day = All.uo.grid.current.in.day,
                                     All.uo.grid.next.in.day = All.uo.grid.next.in.day,
                                     All.vo.grid.current.in.day = All.vo.grid.current.in.day,
                                     All.vo.grid.next.in.day = All.vo.grid.next.in.day,
                                     All.wo.grid.current.in.day = All.wo.grid.current.in.day,
                                     All.wo.grid.next.in.day = All.wo.grid.next.in.day,
                                     All.rho.grid.current.in.day = All.rho.grid.current.in.day,
                                     All.rho.grid.next.in.day = All.rho.grid.next.in.day,
                                     All.Sal.grid.current.in.day = All.Sal.grid.current.in.day,
                                     All.Sal.grid.next.in.day = All.Sal.grid.next.in.day,
                                     All.sd.grid.current.in.day = All.sd.grid.current.in.day,
                                     All.sd.grid.next.in.day = All.sd.grid.next.in.day,
                                     input.func = input.func,
                                     input.data.location = input.data.location,
                                     input.stoke.drift = input.stoke.drift,
                                     timestep = timestep.set,
                                     parms.total = parms.total,
                                     individual_per_core = individual_per_core,
                                     current_month = current_month)
          })
        }
        
        # combine information from all cores to generate a dataframe in each time step
        results.all.cores.eachtimestep.age.0 <- as.data.frame(do.call(rbind, lapply(1:number_cores, function(kkk) {
          results[[kkk]][[1]]
        })))
        results.all.cores.eachtimestep.age.0 <- results.all.cores.eachtimestep.age.0[order(results.all.cores.eachtimestep.age.0$No_ind), ]
        
        results.all.cores.eachtimestep.age.1 <- as.data.frame(do.call(rbind, lapply(1:number_cores, function(kkk) {
          results[[kkk]][[2]]
        })))
        results.all.cores.eachtimestep.age.1 <- results.all.cores.eachtimestep.age.1[order(results.all.cores.eachtimestep.age.1$No_ind), ]
        
        results.all.cores.eachtimestep.age.2 <- as.data.frame(do.call(rbind, lapply(1:number_cores, function(kkk) {
          results[[kkk]][[3]]
        })))
        results.all.cores.eachtimestep.age.2 <- results.all.cores.eachtimestep.age.2[order(results.all.cores.eachtimestep.age.2$No_ind), ]
        
        results.all.cores.eachtimestep.age.3 <- as.data.frame(do.call(rbind, lapply(1:number_cores, function(kkk) {
          results[[kkk]][[4]]
        })))
        results.all.cores.eachtimestep.age.3 <- results.all.cores.eachtimestep.age.3[order(results.all.cores.eachtimestep.age.3$No_ind), ]
        
        results.all.cores.eachtimestep.age.4 <- as.data.frame(do.call(rbind, lapply(1:number_cores, function(kkk) {
          results[[kkk]][[5]]
        })))
        results.all.cores.eachtimestep.age.4 <- results.all.cores.eachtimestep.age.4[order(results.all.cores.eachtimestep.age.4$No_ind), ]
        
        # combine all SIs including all age classes in each time step （only used to judge spawning events）
        results_total_eachtimestep_all_ages<-rbind(results.all.cores.eachtimestep.age.0,
                                                   results.all.cores.eachtimestep.age.1,
                                                   results.all.cores.eachtimestep.age.2,
                                                   results.all.cores.eachtimestep.age.3,
                                                   results.all.cores.eachtimestep.age.4)
        
        if(results_total_eachtimestep_all_ages %>%
           filter(running_status == "alive") %>%
           summarise(total_eggs = sum(eggs, na.rm = TRUE)) %>%
           pull(total_eggs) > 0){        # this is important to avoid zombia SIs till spawning eggs
          # spawn occur
          
          # spawned eggs should be adjusted by worth!!! we did not do it in DEB_function, and we will do it here!
          results.all.cores.eachtimestep.age.0$eggs<-results.all.cores.eachtimestep.age.0$eggs*results.all.cores.eachtimestep.age.0$worth
          results.all.cores.eachtimestep.age.1$eggs<-results.all.cores.eachtimestep.age.1$eggs*results.all.cores.eachtimestep.age.1$worth
          results.all.cores.eachtimestep.age.2$eggs<-results.all.cores.eachtimestep.age.2$eggs*results.all.cores.eachtimestep.age.2$worth
          results.all.cores.eachtimestep.age.3$eggs<-results.all.cores.eachtimestep.age.3$eggs*results.all.cores.eachtimestep.age.3$worth
          results.all.cores.eachtimestep.age.4$eggs<-results.all.cores.eachtimestep.age.4$eggs*results.all.cores.eachtimestep.age.4$worth
          
          results_total_eachtimestep_all_ages<-rbind(results.all.cores.eachtimestep.age.0,
                                                     results.all.cores.eachtimestep.age.1,
                                                     results.all.cores.eachtimestep.age.2,
                                                     results.all.cores.eachtimestep.age.3,
                                                     results.all.cores.eachtimestep.age.4)
          
          # summary how many eggs spawned in each region at each time step
          spawn_eachtimestep<-results_total_eachtimestep_all_ages%>%
            filter(Cell_num!=0 & running_status == 'alive')%>%   # this step is important to defend SI dead but frozen by always spawning eggs
            left_join(blacksea_points[,c('Cell_num','EEZs_adj')],by = 'Cell_num')%>%
            group_by(EEZs_adj)%>%
            summarise(eggs_sum=sum(eggs,na.rm = TRUE))%>%
            as.data.frame()
          
          new_colnames<-as.character(spawn_eachtimestep$EEZs_adj)
          new_row<-t(spawn_eachtimestep$eggs_sum)
          aaa<-as.data.frame(new_row)
          temporal_spawning_dynamic_monitor<-cbind(date_current,aaa)
          colnames(temporal_spawning_dynamic_monitor)<-c('Date',new_colnames)
          temporal_spawning_dynamic_monitor$running_year<-results_total_eachtimestep_all_ages[1,'running_year']
          temporal_spawning_dynamic_monitor$running_pattern<-results_total_eachtimestep_all_ages[1,'running_pattern']
          temporal_spawning_dynamic_monitor$number_of_spawner_x_worth<-sum(results_total_eachtimestep_all_ages[which(results_total_eachtimestep_all_ages$eggs>0),'worth'])
          temporal_spawning_dynamic_monitor$total_eggs_x_worth<-sum(results_total_eachtimestep_all_ages[which(results_total_eachtimestep_all_ages$eggs>0),'eggs'])
          temporal_spawning_dynamic_monitor<-temporal_spawning_dynamic_monitor[colnames(spawning_dynamic_monitor)]
          
          # 【output_1】update and accumulate spawning dynamic monitor
          spawning_dynamic_monitor<-rbind(spawning_dynamic_monitor,temporal_spawning_dynamic_monitor)
          saveRDS(spawning_dynamic_monitor,paste0(path_save,output_path_name,current_test_name,'/spawning_dynamic_monitor.rds'))
          
          spawning_eachtimestep<-rbind(spawning_eachtimestep,temporal_spawning_dynamic_monitor)
        }else{
          # no spawning occur
          # nothing happened
        }
        
        monitor_eachtimestep<-as.data.frame(results_total_eachtimestep_all_ages %>%
                                              mutate(biomass = ifelse(running_status == 'alive', worth * Www, NA_real_)) %>%
                                              group_by(age)%>%
                                              summarise(Date=first(Date),
                                                        Dist.mean=mean(ifelse(running_status == 'alive', Dist, NA_real_), na.rm = TRUE),
                                                        depth.mean=mean(ifelse(running_status == 'alive', depth, NA_real_), na.rm = TRUE),
                                                        dis.to.coast.mean=mean(ifelse(running_status == 'alive', dis_to_coast, NA_real_), na.rm = TRUE),
                                                        Temp.mean=mean(ifelse(running_status == 'alive', Temp, NA_real_), na.rm = TRUE),
                                                        Temp.spawn.mean=mean(ifelse(running_status == 'alive', Temp_spawn, NA_real_), na.rm = TRUE),
                                                        Z.micro.zoo.mean=mean(ifelse(running_status == 'alive', Z_micro_zoo, NA_real_), na.rm = TRUE),
                                                        Z.meso.zoo.mean=mean(ifelse(running_status == 'alive', Z_meso_zoo, NA_real_), na.rm = TRUE),
                                                        f.micro.zoo.mean=mean(ifelse(running_status == 'alive', f_micro_zoo, NA_real_), na.rm = TRUE),
                                                        f.meso.zoo.mean=mean(ifelse(running_status == 'alive', f_meso_zoo, NA_real_), na.rm = TRUE),
                                                        Food.f.mean=mean(ifelse(running_status == 'alive', f, NA_real_), na.rm = TRUE),TL.mean=mean(ifelse(running_status == 'alive', TL, NA), na.rm = TRUE),
                                                        Www.mean=mean(ifelse(running_status == 'alive', Www, NA_real_), na.rm = TRUE),
                                                        eggs.mean=mean(ifelse(running_status == 'alive', eggs, NA_real_), na.rm = TRUE),
                                                        worth.mean=mean(ifelse(running_status == 'alive', worth, NA_real_), na.rm = TRUE),
                                                        total.biomass=sum(biomass, na.rm = TRUE),
                                                        total.worth=sum(ifelse(running_status == 'alive', worth, NA_real_), na.rm = TRUE),
                                                        density.mean=mean(ifelse(running_status == 'alive', density, NA_real_), na.rm = TRUE),
                                                        norm.density.mean=mean(ifelse(running_status == 'alive', norm_density, NA_real_), na.rm = TRUE),
                                                        natural.mortality.mean=mean(ifelse(running_status == 'alive', natural_mortality, NA_real_), na.rm = TRUE),
                                                        ZFAC.mean=mean(ifelse(running_status == 'alive', ZFAC, NA_real_), na.rm = TRUE),
                                                        final.mortality.mean=mean(ifelse(running_status == 'alive', final_mortality, NA_real_), na.rm = TRUE),
                                                        total.eggs.worth=sum(ifelse(running_status == 'alive', eggs, NA_real_), na.rm = TRUE),
                                                        running.pattern=running_pattern[1],
                                                        running.year=running_year[1],
                                                        number.adults=sum(ifelse(life_stage == 'adult', worth, NA_real_), na.rm = TRUE),
                                                        Num_alive_SI=sum(running_status == 'alive', na.rm = TRUE),
                                                        remainder_rate_SI=(sum(running_status == 'alive', na.rm = TRUE)/number_individual),
                                                        remainder_rate_SI_worth=(sum(ifelse(running_status == "alive", worth, 0), na.rm = TRUE)/sum(worth, na.rm = TRUE)),
                                                        Num_Killed_by_stunted=sum(running_status == 'Killed_by_stunted', na.rm = TRUE),
                                                        Num_Killed_by_low_worth=sum(running_status == 'Killed_by_low_worth', na.rm = TRUE)))
        
        # 【Output_2】monitor data frame
        # update global_monitor for each timestep
        global_monitor_eachtimestep<-rbind(global_monitor_eachtimestep,monitor_eachtimestep)
        
        saveRDS(global_monitor_eachtimestep,paste0(path_save,output_path_name,current_test_name,'/global_monitor_eachtimestep.rds'))
        
        # randomly find typical individuals in a certain year June, July, and August 
        if(i==6){
          information_selected_individuals<-rbind(information_selected_individuals,results.all.cores.eachtimestep.age.0[target_ind_row_num,])
        }else if(i==7){
          information_selected_individuals<-rbind(information_selected_individuals,results.all.cores.eachtimestep.age.1[target_ind_row_num,])
        }else if(i==8){
          information_selected_individuals<-rbind(information_selected_individuals,results.all.cores.eachtimestep.age.2[target_ind_row_num,])
        }else if(i==9){
          information_selected_individuals<-rbind(information_selected_individuals,results.all.cores.eachtimestep.age.3[target_ind_row_num,])
        }else if(i==10){
          information_selected_individuals<-rbind(information_selected_individuals,results.all.cores.eachtimestep.age.4[target_ind_row_num,])
        }
        
        if(timestep.set>1 & m != timestep.set){ # save output in each timestep (30 min, 20 min, or 10 min)
          # collect all data frame to the final list and can be used as input state variable for the next time step
          results_total_age_0<-c(results_total_age_0,list(results.all.cores.eachtimestep.age.0))
          results_total_age_1<-c(results_total_age_1,list(results.all.cores.eachtimestep.age.1))
          results_total_age_2<-c(results_total_age_2,list(results.all.cores.eachtimestep.age.2))
          results_total_age_3<-c(results_total_age_3,list(results.all.cores.eachtimestep.age.3))
          results_total_age_4<-c(results_total_age_4,list(results.all.cores.eachtimestep.age.4))
        }
        
      } # end loop for all time steps
      
      spawning_eachhour<-rbind(spawning_eachhour,spawning_eachtimestep)
      
      if(timestep.set == 1 & k != 24){ # hourly to save as total results
        # collect all data frame to the final list and can be used as input state variable for the next time step
        results_total_age_0<-c(results_total_age_0,list(results.all.cores.eachtimestep.age.0))
        results_total_age_1<-c(results_total_age_1,list(results.all.cores.eachtimestep.age.1))
        results_total_age_2<-c(results_total_age_2,list(results.all.cores.eachtimestep.age.2))
        results_total_age_3<-c(results_total_age_3,list(results.all.cores.eachtimestep.age.3))
        results_total_age_4<-c(results_total_age_4,list(results.all.cores.eachtimestep.age.4))
      }
      
    } # end loop for all hours
    
    #### ONE DAY START
    
    # examine the duration of egg and yolksac larvae =======================================
    # According to observations, the pre-larval period lasts for about 2.0 to 2.5 d at water temperatures above 20°C (Lisovenko and Andrianov 1996).
    # Therefore, when temperature = 20 degree, we set developing period of egg is 24 hours, and developing period of yolksac larvae is 36 hours
    # So, we used 4.25*(10^-5)*(20^2.3)*24 ~ 1; and 4.2*(10^-5)*(20^2.17)*36 ~ 1 to calibrate the function
    results.all.cores.eachtimestep.age.0<-results.all.cores.eachtimestep.age.0%>%
      mutate(jump_egg_days=case_when(
        running_status=='alive' & life_stage =='egg' ~ jump_egg_days + 4.25*(10^-5)*(Temp^2.3)*24,
        TRUE ~ jump_egg_days
      ))%>%
      mutate(jump_yolksac_days=case_when(
        running_status=='alive' & life_stage =='yolksac_larvae' ~ jump_yolksac_days + 4.25*(10^-5)*(Temp^2.17)*24,
        TRUE ~ jump_yolksac_days
      ))%>%
      as.data.frame()
    
    # switch life_stage - egg to yolksac_larvae to early_larvae SIs
    results.all.cores.eachtimestep.age.0 <- results.all.cores.eachtimestep.age.0 %>%
      mutate(life_stage = case_when(
        running_status=='alive' & life_stage == 'egg' & jump_egg_days >= 1 ~ 'yolksac_larvae',
        running_status=='alive' & life_stage == 'yolksac_larvae' & jump_yolksac_days >= 1 ~ 'early_larvae',
        TRUE ~ life_stage
      ))
    
    # update running survival duration for each individual =========================
    results.all.cores.eachtimestep.age.0[which(results.all.cores.eachtimestep.age.0$running_status%in%'alive'),'survival_duration']<-
      results.all.cores.eachtimestep.age.0[which(results.all.cores.eachtimestep.age.0$running_status%in%'alive'),'survival_duration']+1
    
    results.all.cores.eachtimestep.age.1[which(results.all.cores.eachtimestep.age.1$running_status%in%'alive'),'survival_duration']<-
      results.all.cores.eachtimestep.age.1[which(results.all.cores.eachtimestep.age.1$running_status%in%'alive'),'survival_duration']+1
    
    results.all.cores.eachtimestep.age.2[which(results.all.cores.eachtimestep.age.2$running_status%in%'alive'),'survival_duration']<-
      results.all.cores.eachtimestep.age.2[which(results.all.cores.eachtimestep.age.2$running_status%in%'alive'),'survival_duration']+1
    
    results.all.cores.eachtimestep.age.3[which(results.all.cores.eachtimestep.age.3$running_status%in%'alive'),'survival_duration']<-
      results.all.cores.eachtimestep.age.3[which(results.all.cores.eachtimestep.age.3$running_status%in%'alive'),'survival_duration']+1
    
    results.all.cores.eachtimestep.age.4[which(results.all.cores.eachtimestep.age.4$running_status%in%'alive'),'survival_duration']<-
      results.all.cores.eachtimestep.age.4[which(results.all.cores.eachtimestep.age.4$running_status%in%'alive'),'survival_duration']+1
    
    # spawning events occur =============================
    if(nrow(spawning_eachhour)>0){ 
      
      # delete unnecessary columns
      spawning_eachhour_slim<-spawning_eachhour%>%
        dplyr::select(-running_year, -running_pattern, -number_of_spawner_x_worth, -total_eggs_x_worth)%>%
        summarise(across(-Date, sum)) %>%                           # ignore first col, calculate sum of the other cols
        mutate(Date = last(spawning_eachhour$Date)) %>%         
        dplyr::select(Date, everything()) 
      
      # allocate SIs for each region
      spawn_with_region<-gather(spawning_eachhour_slim,key = 'EEZs_adj',value = 'eggs_sum',
                                -'Date')          # transfer from horizontal to vertical
      spawn_with_region<-spawn_with_region[which(spawn_with_region$eggs_sum>0),]
      spawn_with_region$allocate_worth<-round(spawn_with_region$eggs_sum / number_SI_for_each_region,digits = 0) # allocated worth to each SI 
      
      # judge whether we have enough available places in different regions
      potential_points<-blacksea_points%>%
        filter(EEZs_adj%in%unique(spawn_with_region$EEZs_adj))%>%
        filter(Depth <= 200)%>%
        as.data.frame()
      
      # update Temp for new-born eggs
      # because new eggs only located at the cell center without spatial interpolation
      Temp.egg<-get_Temp_new_born_eggs_2D_cell_center(updated_date = spawn_with_region$Date[1], 
                                                      cell_number = potential_points$Cell_num)
      potential_points<-cbind(potential_points,Temp.egg)
      
      # start to judge
      judge_avail_points<-potential_points%>%
        mutate(EEZs_adj = as.factor(EEZs_adj)) %>% 
        filter(Temp_0_10m >= 16)%>%
        group_by(EEZs_adj)%>%
        summarize(avail_num_cell=n())%>%
        right_join(potential_points %>% distinct(EEZs_adj), by = "EEZs_adj") %>%
        replace_na(list(avail_num_cell = 0)) %>%  # 将缺失值替换为 0
        mutate(judge_avail = case_when(
          avail_num_cell < number_SI_for_each_region ~ FALSE,
          TRUE ~ TRUE
        ))
      
      # combine the judging results with spawning activities 
      combine_judge_spawn_with_region<-spawn_with_region%>%
        left_join(judge_avail_points,by = 'EEZs_adj')%>%
        filter(judge_avail)%>% # we have already write TRUE or FALSE in the previous step
        as.data.frame()
      
      avail_SI_age_0<-results.all.cores.eachtimestep.age.0[which(results.all.cores.eachtimestep.age.0$running_status%in%c('silence')),]
      num_allocate_SI<-nrow(combine_judge_spawn_with_region) * number_SI_for_each_region
      avail_SI_age_0<-avail_SI_age_0[1:num_allocate_SI,]
      
      # here, we need to allocate SI randomly distribute places within the regions but all should higher than 16 degree
      # otherwise, SI will died quickly because cannot grow up to early_larvae stage
      # extract Temp from BAMHBI dataset and combine with potential points
      
      set.seed(timepoint)
      random_location_spawn_region<-potential_points%>%
        filter(Temp_0_10m >= 16 & Depth <= 200)%>%
        group_by(EEZs_adj)%>%
        slice_sample(n = number_SI_for_each_region,replace = FALSE)%>%
        left_join(spawn_with_region,by = 'EEZs_adj')%>%       # do not need replicate, they will generate replicated rows automatically
        as.data.frame()
      
      # update the selected rows in results.all.cores.eachtimestep.age.0 
      # and assign initial values to the available SIs
      avail_SI_age_0[c('Lon','Lat','Cell_num','worth','EEZs_region','Temp')]<-
        random_location_spawn_region[c('Longitude','Latitude','Cell_num','allocate_worth','EEZs_adj','Temp_0_10m')]
      avail_SI_age_0[,'dis_to_coast'] <- blacksea_points$dis_to_coast[
        match(avail_SI_age_0[,'Cell_num'], blacksea_points$Cell_num)
      ]
      avail_SI_age_0[,'E']<-parms.total['E0']
      avail_SI_age_0[,'V']<-parms.total['V0']
      avail_SI_age_0[,'TL']<-0.36
      avail_SI_age_0[,'running_status']<-'alive'
      avail_SI_age_0[,'life_stage']<-'egg'
      
      # generate depth in a random position, this may be updated in 3D !!!!!!!!!!!!!!!!!!!!!!###############$$$$$$$$$$$$$$$
      avail_SI_age_0[,'depth']<-sample(seq(0,10,by=0.01),size = nrow(avail_SI_age_0),replace = FALSE)
      avail_SI_age_0[,'depth_level']<-define_depth_level(depth = avail_SI_age_0[,'depth'])
      # update the total age-0 data frame
      results.all.cores.eachtimestep.age.0[results.all.cores.eachtimestep.age.0$No_ind %in%
                                             avail_SI_age_0$No_ind,] <- avail_SI_age_0
      
      # 【Output_3】summed spawning egg information on each day
      spawning_eachday <- spawning_eachhour %>%
        summarise(across(-c(Date,running_year,running_pattern), sum)) %>%                           # ignore first col, calculate sum of the other cols
        mutate(Date = last(spawning_eachhour$Date),           # get the final row of Date
               running_year = last(spawning_eachhour$running_year),
               running_pattern = last(spawning_eachhour$running_pattern)) %>%         
        select(Date, everything()) 
      
      spawning_eachday_output<-rbind(spawning_eachday_output,spawning_eachday)
      
      saveRDS(spawning_eachday_output,paste0(path_save,output_path_name,current_test_name,'/spawning_eachday.rds'))
      
    }
    
    # special event for the last day of spawning =================================
    if(format(current.Date, "%Y-%m-%d") %in% as.character(spawn_day$latest_day)){
      # the last day to spawn that all available SIs will be used by divided worth
      
      # allocate SIs for each region
      
      avail_SI_age_0<-results.all.cores.eachtimestep.age.0[which(results.all.cores.eachtimestep.age.0$running_status%in%c('silence')),]
      potential.name<-colnames(avail_SI_age_0)
      potential.name<-potential.name[potential.name != 'No_ind']
      for (avai_SI in 1:nrow(avail_SI_age_0)) {
        aaa <- results.all.cores.eachtimestep.age.0 %>%
          filter(running_status == 'alive') %>%
          slice(which.max(worth))
        
        # do not overwrite the No_ind var
        avail_SI_age_0[avai_SI, potential.name] <- aaa[potential.name]
        max_worth_index <- aaa[1,'No_ind'] # make sure only one max value
        results.all.cores.eachtimestep.age.0[which(results.all.cores.eachtimestep.age.0$No_ind==max_worth_index), 'worth'] <- 
          round(results.all.cores.eachtimestep.age.0[which(results.all.cores.eachtimestep.age.0$No_ind==max_worth_index), 'worth'] / 2,digits = 0)
        avail_SI_age_0[avai_SI, 'worth'] <- round(avail_SI_age_0[avai_SI, 'worth'] / 2,digits = 0)
        
        # allocate cell number near the target SI 
        if (avail_SI_age_0[avai_SI, 'Cell_num'] == 141101) { # Attention! 141101 is the only one without neighbor points
          allocate.cell.num <- 140523
        }else{
          allocate.cell.num <- avail_SI_age_0[avai_SI, 'Cell_num'] + 1
          if (!(allocate.cell.num %in% blacksea_points$Cell_num)) {
            allocate.cell.num <- avail_SI_age_0[avai_SI, 'Cell_num'] - 1
          }
        }
        avail_SI_age_0[avai_SI, 'Cell_num'] <- allocate.cell.num
        
        avail_SI_age_0[avai_SI,'Lat']<-blacksea_points[which(blacksea_points$Cell_num%in%allocate.cell.num),'Latitude']
        avail_SI_age_0[avai_SI,'Lon']<-blacksea_points[which(blacksea_points$Cell_num%in%allocate.cell.num),'Longitude']
        results.all.cores.eachtimestep.age.0[which(results.all.cores.eachtimestep.age.0$No_ind %in%
                                                     avail_SI_age_0[avai_SI,'No_ind']),] <- avail_SI_age_0[avai_SI,]
      }
      
    }else{
      # not the last day of spawning window
      # nothing happened
    }
    
    
    # output weight- and length-at-age, moving latitude and longitude-at-age and spawned eggs ============================
    # with all age combined on each day before adding mortality 
    if(i>=6){ # formal historical simulation
      
      results_total_eachtimestep_all_ages<-rbind(results.all.cores.eachtimestep.age.0,
                                                 results.all.cores.eachtimestep.age.1,
                                                 results.all.cores.eachtimestep.age.2,
                                                 results.all.cores.eachtimestep.age.3,
                                                 results.all.cores.eachtimestep.age.4)
      
      # 1. calculate weight-at-age, length-at-age, and moving latitude and longitude-at-age
      # weight and length
      temporal_weight_age <- results_total_eachtimestep_all_ages %>%
        filter(running_status=='alive')%>%
        group_by(Date,running_year,age) %>%
        summarize(
          min = min(Www, na.rm = TRUE),
          q1 = quantile(Www, 0.25, na.rm = TRUE),
          median = median(Www, na.rm = TRUE),
          q3 = quantile(Www, 0.75, na.rm = TRUE),
          max = max(Www, na.rm = TRUE),
          mean = mean(Www,na.rm = TRUE),
          sd = sd(Www,na.rm = TRUE)
        )
      
      temporal_length_age <- results_total_eachtimestep_all_ages %>%
        filter(running_status=='alive')%>%
        group_by(Date,running_year,age) %>%
        summarize(
          min = min(TL, na.rm = TRUE),
          q1 = quantile(TL, 0.25, na.rm = TRUE),
          median = median(TL, na.rm = TRUE),
          q3 = quantile(TL, 0.75, na.rm = TRUE),
          max = max(TL, na.rm = TRUE),
          mean = mean(TL,na.rm = TRUE),
          sd = sd(TL,na.rm = TRUE)
        )
      
      annual_weight_age<-rbind(annual_weight_age,temporal_weight_age)
      annual_length_age<-rbind(annual_length_age,temporal_length_age)
      
      # moving latitude and longitude
      aaa_Lon_Lat_age <- results_total_eachtimestep_all_ages %>%
        filter(running_status=='alive')%>%
        filter(age!=0)%>%
        group_by(Date,running_year,age) %>%
        summarize(
          Lon_mean = mean(Lon, na.rm = TRUE),
          Lat_mean = mean(Lat, na.rm = TRUE))%>%
        mutate(age=as.character(age))
      
      aaa_total_Lon_Lat_age<-results_total_eachtimestep_all_ages%>%
        filter(running_status=='alive')%>%
        filter(age!=0)%>%
        group_by(Date,running_year)%>%
        summarize(
          age='Total',
          Lon_mean = mean(Lon, na.rm = TRUE),
          Lat_mean = mean(Lat, na.rm = TRUE))
      
      temporal_Lon_Lat_age<-bind_rows(aaa_Lon_Lat_age,aaa_total_Lon_Lat_age)%>%
        arrange(Date)
      
      annual_Lon_Lat_age<-rbind(annual_Lon_Lat_age,temporal_Lon_Lat_age)
      # the files will be saved at the end of the year
      
      # 2. mean spawned eggs per female
      summary_eggs_aaa<-results_total_eachtimestep_all_ages%>%
        filter(running_status=='alive')%>%
        group_by(age, No_ind)%>%
        summarise(all_batches=sum(eggs>0,na.rm = TRUE),
                  all_eggs = sum(eggs, na.rm = TRUE))%>%
        as.data.frame()
      
      eggs_batch_num_age_0[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==0),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==0),'all_batches']
      eggs_batch_num_age_1[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==1),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==1),'all_batches']
      eggs_batch_num_age_2[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==2),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==2),'all_batches']
      eggs_batch_num_age_3[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==3),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==3),'all_batches']
      eggs_batch_num_age_4[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==4),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==4),'all_batches']
      
      eggs_total_num_age_0[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==0),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==0),'all_eggs']
      eggs_total_num_age_1[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==1),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==1),'all_eggs']
      eggs_total_num_age_2[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==2),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==2),'all_eggs']
      eggs_total_num_age_3[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==3),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==3),'all_eggs']
      eggs_total_num_age_4[day_counter,summary_eggs_aaa[which(summary_eggs_aaa$age==4),'No_ind']]<-summary_eggs_aaa[which(summary_eggs_aaa$age==4),'all_eggs']
      
    }
    
    # ADD mortality =========================================
    # Why here? because eggs per female cannot be adjusted by the worth for this part, 
    # otherwise, you never can statistic the real number of the spawned eggs
    # if we add here, we use the initial worth to spawn eggs, and later we summary the information of spatial map of biomass and abundance adjusted by worth
    
    # ADD mortality to the worth of each SI
    # determine Density-Dependent mortality for different life stage
    
    if(input.dimension == 1){ # 2D
      total.mortality.age.0<-natural_and_DD_mortality_2D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.0)
      total.mortality.age.1<-natural_and_DD_mortality_2D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.1)
      total.mortality.age.2<-natural_and_DD_mortality_2D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.2)
      total.mortality.age.3<-natural_and_DD_mortality_2D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.3)
      total.mortality.age.4<-natural_and_DD_mortality_2D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.4)
      
      results.all.cores.eachtimestep.age.0<-total.mortality.age.0[[1]]
      results.all.cores.eachtimestep.age.1<-total.mortality.age.1[[1]]
      results.all.cores.eachtimestep.age.2<-total.mortality.age.2[[1]]
      results.all.cores.eachtimestep.age.3<-total.mortality.age.3[[1]]
      results.all.cores.eachtimestep.age.4<-total.mortality.age.4[[1]]
      
      overlap.occur.one.day<-rbind(total.mortality.age.0[[2]],
                                   total.mortality.age.1[[2]],
                                   total.mortality.age.2[[2]],
                                   total.mortality.age.3[[2]],
                                   total.mortality.age.4[[2]])
      
      if(!is.null(overlap.occur.one.day)){
        overlap.occur.all.days<-rbind(overlap.occur.all.days,
                                      overlap.occur.one.day)
        
      }
      
    }else if (input.dimension == 2){ # 3D
      total.mortality.age.0<-natural_and_DD_mortality_3D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.0)
      total.mortality.age.1<-natural_and_DD_mortality_3D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.1)
      total.mortality.age.2<-natural_and_DD_mortality_3D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.2)
      total.mortality.age.3<-natural_and_DD_mortality_3D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.3)
      total.mortality.age.4<-natural_and_DD_mortality_3D(input.DD.mortality,
                                                         results.all.cores.eachtimestep.age.4)
      
      results.all.cores.eachtimestep.age.0<-total.mortality.age.0[[1]]
      results.all.cores.eachtimestep.age.1<-total.mortality.age.1[[1]]
      results.all.cores.eachtimestep.age.2<-total.mortality.age.2[[1]]
      results.all.cores.eachtimestep.age.3<-total.mortality.age.3[[1]]
      results.all.cores.eachtimestep.age.4<-total.mortality.age.4[[1]]
      
      overlap.occur.one.day<-rbind(total.mortality.age.0[[2]],
                                   total.mortality.age.1[[2]],
                                   total.mortality.age.2[[2]],
                                   total.mortality.age.3[[2]],
                                   total.mortality.age.4[[2]])
      
      if(!is.null(overlap.occur.one.day)){
        overlap.occur.all.days<-rbind(overlap.occur.all.days,
                                      overlap.occur.one.day)
        
      }
    }
    
    # output monitor for D-D mortality daily
    DD_age_combine <- rbind(
      results.all.cores.eachtimestep.age.0[c('No_ind', 'age', 'running_status', 'life_stage', 'TL',
                                             'density', 'norm_density', 'natural_mortality', 'ZFAC', 'final_mortality')],
      results.all.cores.eachtimestep.age.1[c('No_ind', 'age', 'running_status', 'life_stage', 'TL',
                                             'density', 'norm_density', 'natural_mortality', 'ZFAC', 'final_mortality')],
      results.all.cores.eachtimestep.age.2[c('No_ind', 'age', 'running_status', 'life_stage', 'TL',
                                             'density', 'norm_density', 'natural_mortality', 'ZFAC', 'final_mortality')],
      results.all.cores.eachtimestep.age.3[c('No_ind', 'age', 'running_status', 'life_stage', 'TL',
                                             'density', 'norm_density', 'natural_mortality', 'ZFAC', 'final_mortality')],
      results.all.cores.eachtimestep.age.4[c('No_ind', 'age', 'running_status', 'life_stage', 'TL',
                                             'density', 'norm_density', 'natural_mortality', 'ZFAC', 'final_mortality')]
    )
    # generate mornitor for ZFAC and density and output once a year
    # IMPORTANT!!! Here density and normalized density are already adjusted by mortality
    mortality_age_monitor<-c(mortality_age_monitor,list(DD_age_combine))
    
    # Kill low-growing individuals, determined by maximum survival days
    kill_low_growing_ind <- function(df) {
      df <- df %>%
        mutate(running_status = case_when(
          running_status=='alive' & life_stage %in% c("egg",'yolksac_larvae') & survival_duration > parms.total['max_day_egg_yolksac'] ~ "Killed_by_stunted",
          running_status=='alive' & life_stage == "early_larvae" & survival_duration > parms.total['max_day_early_larvae'] ~ "Killed_by_stunted",
          running_status=='alive' & life_stage == "late_larvae" & survival_duration > parms.total['max_day_late_larvae'] ~ "Killed_by_stunted",
          running_status=='alive' & life_stage == "juvenile" & survival_duration > parms.total['max_day_juvenile'] ~ "Killed_by_stunted",
          TRUE ~ running_status
        ))%>%
        mutate(death_causes = case_when(
          running_status=='Killed_by_stunted' & life_stage == "egg" ~ "Stunted_in_egg",
          running_status=='Killed_by_stunted' & life_stage == 'yolksac_larvae' ~ "Stunted_in_yolksac_larvae",
          running_status=='Killed_by_stunted' & life_stage == "early_larvae"  ~ "Stunted_in_early_larvae",
          running_status=='Killed_by_stunted' & life_stage == "late_larvae"  ~ "Stunted_in_late_larvae",
          running_status=='Killed_by_stunted' & life_stage == "juvenile"  ~ "Stunted_in_juvenile",
          TRUE ~ death_causes
        ))
      return(df)
    }
    
    results.all.cores.eachtimestep.age.0<-kill_low_growing_ind(results.all.cores.eachtimestep.age.0)
    results.all.cores.eachtimestep.age.1<-kill_low_growing_ind(results.all.cores.eachtimestep.age.1)
    results.all.cores.eachtimestep.age.2<-kill_low_growing_ind(results.all.cores.eachtimestep.age.2)
    results.all.cores.eachtimestep.age.3<-kill_low_growing_ind(results.all.cores.eachtimestep.age.3)
    results.all.cores.eachtimestep.age.4<-kill_low_growing_ind(results.all.cores.eachtimestep.age.4)
    
    # Kill very low values of worth (Kill individuals with worth < 1)
    kill_low_worth_ind <- function(df) {
      df <- df %>%
        mutate(running_status = case_when(
          running_status=='alive' & worth < 1 ~ "Killed_by_low_worth",
          TRUE ~ running_status
        ))%>%
        mutate(death_causes = case_when(
          running_status=='Killed_by_low_worth' & life_stage == "egg" ~ "Low_worth_in_egg",
          running_status=='Killed_by_low_worth' & life_stage == 'yolksac_larvae' ~ "Low_worth_in_yolksac_larvae",
          running_status=='Killed_by_low_worth' & life_stage == "early_larvae"  ~ "Low_worth_in_early_larvae",
          running_status=='Killed_by_low_worth' & life_stage == "late_larvae"  ~ "Low_worth_in_late_larvae",
          running_status=='Killed_by_low_worth' & life_stage == "juvenile"  ~ "Low_worth_in_juvenile",
          running_status=='Killed_by_low_worth' & life_stage == "adult"  ~ "Low_worth_in_adult",
          TRUE ~ death_causes
        ))%>%
        as.data.frame()
      return(df)
    }
    
    results.all.cores.eachtimestep.age.0<-kill_low_worth_ind(results.all.cores.eachtimestep.age.0)
    results.all.cores.eachtimestep.age.1<-kill_low_worth_ind(results.all.cores.eachtimestep.age.1)
    results.all.cores.eachtimestep.age.2<-kill_low_worth_ind(results.all.cores.eachtimestep.age.2)
    results.all.cores.eachtimestep.age.3<-kill_low_worth_ind(results.all.cores.eachtimestep.age.3)
    results.all.cores.eachtimestep.age.4<-kill_low_worth_ind(results.all.cores.eachtimestep.age.4)
    
    # collect all data frame to the final list and can be used as input state variable for the next time step
    results_total_age_0<-c(results_total_age_0,list(results.all.cores.eachtimestep.age.0))
    results_total_age_1<-c(results_total_age_1,list(results.all.cores.eachtimestep.age.1))
    results_total_age_2<-c(results_total_age_2,list(results.all.cores.eachtimestep.age.2))
    results_total_age_3<-c(results_total_age_3,list(results.all.cores.eachtimestep.age.3))
    results_total_age_4<-c(results_total_age_4,list(results.all.cores.eachtimestep.age.4))
    
    # 【output_10】global monitor for each day (adding the final reproduction and worth with all adjusted by mortality)
    
    # IMPORTANT!!! Variables in eachday_monitor: density and eggs etc. are all adjusted by mortality.
    
    results_total_eachtimestep_all_ages<-rbind(results.all.cores.eachtimestep.age.0,
                                               results.all.cores.eachtimestep.age.1,
                                               results.all.cores.eachtimestep.age.2,
                                               results.all.cores.eachtimestep.age.3,
                                               results.all.cores.eachtimestep.age.4)
    
    monitor_eachday<-results_total_eachtimestep_all_ages %>%
      mutate(biomass = ifelse(running_status == 'alive', worth * Www, NA_real_)) %>%
      group_by(age)%>%
      summarise(Date=first(Date),
                Dist.mean=mean(ifelse(running_status == 'alive', Dist, NA_real_), na.rm = TRUE),
                depth.mean=mean(ifelse(running_status == 'alive', depth, NA_real_), na.rm = TRUE),
                dis.to.coast.mean=mean(ifelse(running_status == 'alive', dis_to_coast, NA_real_), na.rm = TRUE),
                Temp.mean=mean(ifelse(running_status == 'alive', Temp, NA_real_), na.rm = TRUE),
                Temp.spawn.mean=mean(ifelse(running_status == 'alive', Temp_spawn, NA_real_), na.rm = TRUE),
                Z.micro.zoo.mean=mean(ifelse(running_status == 'alive', Z_micro_zoo, NA_real_), na.rm = TRUE),
                Z.meso.zoo.mean=mean(ifelse(running_status == 'alive', Z_meso_zoo, NA_real_), na.rm = TRUE),
                f.micro.zoo.mean=mean(ifelse(running_status == 'alive', f_micro_zoo, NA_real_), na.rm = TRUE),
                f.meso.zoo.mean=mean(ifelse(running_status == 'alive', f_meso_zoo, NA_real_), na.rm = TRUE),
                Food.f.mean=mean(ifelse(running_status == 'alive', f, NA_real_), na.rm = TRUE),
                TL.mean=mean(ifelse(running_status == 'alive', TL, NA_real_), na.rm = TRUE),
                Www.mean=mean(ifelse(running_status == 'alive', Www, NA_real_), na.rm = TRUE),
                eggs.mean=mean(ifelse(running_status == 'alive', eggs, NA_real_), na.rm = TRUE),
                worth.mean=mean(ifelse(running_status == 'alive', worth, NA_real_), na.rm = TRUE),
                total.biomass=sum(biomass, na.rm = TRUE),
                total.worth=sum(ifelse(running_status == 'alive', worth, NA_real_), na.rm = TRUE),
                density.mean=mean(ifelse(running_status == 'alive', density, NA_real_), na.rm = TRUE),
                norm.density.mean=mean(ifelse(running_status == 'alive', norm_density, NA_real_), na.rm = TRUE),
                natural.mortality.mean=mean(ifelse(running_status == 'alive', natural_mortality, NA_real_), na.rm = TRUE),
                ZFAC.mean=mean(ifelse(running_status == 'alive', ZFAC, NA_real_), na.rm = TRUE),
                final.mortality.mean=mean(ifelse(running_status == 'alive', final_mortality, NA_real_), na.rm = TRUE),
                total.eggs.worth=sum(ifelse(running_status == 'alive', eggs, NA_real_), na.rm = TRUE),
                running.pattern=running_pattern[1],
                running.year=running_year[1],
                number.adults=sum(ifelse(life_stage == 'adult', worth, NA_real_), na.rm = TRUE),
                Num_alive_SI=sum(running_status == 'alive', na.rm = TRUE),
                remainder_rate_SI=(sum(running_status == 'alive', na.rm = TRUE)/number_individual),
                remainder_rate_SI_worth=(sum(ifelse(running_status == "alive", worth, 0), na.rm = TRUE)/sum(worth, na.rm = TRUE)),
                Num_Killed_by_stunted=sum(running_status == 'Killed_by_stunted', na.rm = TRUE),
                Num_Killed_by_low_worth=sum(running_status == 'Killed_by_low_worth', na.rm = TRUE)) %>%
      as.data.frame()
    
    global_monitor_eachday<-rbind(global_monitor_eachday,monitor_eachday)
    saveRDS(global_monitor_eachday,paste0(path_save,output_path_name,current_test_name,'/global_monitor_eachday.rds'))
    
    # 3. For the spatial abundance for egg and biomass for early-, late-larvae, juvenile, and adult
    if(i>=6){ 
      summary_cell_num<-results_total_eachtimestep_all_ages%>%
        filter(running_status=='alive')%>%
        mutate(biomass=Www*worth)%>%
        group_by(life_stage,Cell_num)%>%
        summarise(sum_abundance = sum(worth,na.rm = TRUE),
                  sum_biomass = sum(biomass, na.rm = TRUE))%>%
        mutate(Cell_num_ord = which(cell.number%in%Cell_num))%>%
        as.data.frame()
      
      total_cell_num<-summary_cell_num%>%
        group_by(Cell_num)%>%
        summarise(total_abundance = sum(sum_abundance,na.rm = TRUE),
                  total_biomass = sum(sum_biomass,na.rm = TRUE))%>%
        mutate(Cell_num_ord = which(cell.number%in%Cell_num))%>%
        as.data.frame()
      
      egg_summary_cell_num<-summary_cell_num[which(summary_cell_num$life_stage=='egg'),]
      yolksac_larvae_summary_cell_num<-summary_cell_num[which(summary_cell_num$life_stage=='yolksac_larvae'),]
      early_larvae_summary_cell_num<-summary_cell_num[which(summary_cell_num$life_stage=='early_larvae'),]
      late_larvae_summary_cell_num<-summary_cell_num[which(summary_cell_num$life_stage=='late_larvae'),]
      juvenile_summary_cell_num<-summary_cell_num[which(summary_cell_num$life_stage=='juvenile'),]
      adult_summary_cell_num<-summary_cell_num[which(summary_cell_num$life_stage=='adult'),]
      
      egg_seasonal_matrix_abundance[day_counter,egg_summary_cell_num$Cell_num_ord]<-egg_summary_cell_num$sum_abundance
      yolksac_larvae_seasonal_matrix_abundance[day_counter,yolksac_larvae_summary_cell_num$Cell_num_ord]<-yolksac_larvae_summary_cell_num$sum_abundance
      
      early_larvae_seasonal_matrix_biomass[day_counter,early_larvae_summary_cell_num$Cell_num_ord]<-early_larvae_summary_cell_num$sum_biomass
      late_larvae_seasonal_matrix_biomass[day_counter,late_larvae_summary_cell_num$Cell_num_ord]<-late_larvae_summary_cell_num$sum_biomass
      juvenile_seasonal_matrix_biomass[day_counter,juvenile_summary_cell_num$Cell_num_ord]<-juvenile_summary_cell_num$sum_biomass
      adult_seasonal_matrix_biomass[day_counter,adult_summary_cell_num$Cell_num_ord]<-adult_summary_cell_num$sum_biomass
      
      early_larvae_seasonal_matrix_abundance[day_counter,early_larvae_summary_cell_num$Cell_num_ord]<-early_larvae_summary_cell_num$sum_abundance
      late_larvae_seasonal_matrix_abundance[day_counter,late_larvae_summary_cell_num$Cell_num_ord]<-late_larvae_summary_cell_num$sum_abundance
      juvenile_seasonal_matrix_abundance[day_counter,juvenile_summary_cell_num$Cell_num_ord]<-juvenile_summary_cell_num$sum_abundance
      adult_seasonal_matrix_abundance[day_counter,adult_summary_cell_num$Cell_num_ord]<-adult_summary_cell_num$sum_abundance
      
      total_seasonal_matrix_abundance[day_counter,total_cell_num$cell_num_ord]<-total_cell_num$total_abundance
      total_seasonal_matrix_biomass[day_counter,total_cell_num$cell_num_ord]<-total_cell_num$total_biomass
      # output RDS file at the end of each year
      
      # 4. generate relationship between spawner stock biomass and recruitment
      if(j==days_in_year){# this relationship will be collected on the end of one year
        SR_summary_collect<-results_total_eachtimestep_all_ages%>%
          filter(running_status=='alive')%>%
          mutate(biomass=Www*worth)%>%
          summarise(Date=first(Date),
                    running.year=first(running_year),
                    running.pattern=first(running_pattern),
                    recruitment_abundance = sum(worth[age < 1],na.rm = TRUE),
                    recruitment_biomass = sum(biomass[age < 1], na.rm = TRUE),
                    SS_abundance = sum(worth[age >= 1 & life_stage == 'adult'],na.rm = TRUE),
                    SS_biomass = sum(biomass[age >= 1 & life_stage == 'adult'],na.rm = TRUE))%>%
          as.data.frame()
        SR_relationship<-rbind(SR_relationship,SR_summary_collect)
        # 【output_12】SR-relationship at the end day of each year
        saveRDS(SR_relationship,paste0(path_save,output_path_name,current_test_name,'/SR_relationship.rds'))
      }
      
      # 5. update annual accumulate recording SI of duration of each life stage, spawning batch, mean dis to shoreline and days in each region
      # group by age、ID and life_stage, calculate accumulated days in different life stages
      
      stage_summary <- results_total_eachtimestep_all_ages %>%
        group_by(age, No_ind, life_stage) %>%
        summarize(
          stage_days = n(),
          eggs_total_num_no_worth = sum(eggs/worth, na.rm = TRUE),
          eggs_total_num_with_worth = sum(eggs, na.rm = TRUE),
          batch_num_no_worth = sum(eggs > 0),
          batch_num_with_worth = sum(eggs > 0)*worth,
          sum_dis_to_coast = sum(dis_to_coast,na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        pivot_wider(
          names_from = life_stage,
          values_from = stage_days,
          values_fill = list(stage_days = 0)
        )  %>%
        mutate(
          egg_days = if ("egg" %in% names(.)) egg else 0,
          yolksac_larvae_days = if ("yolksac_larvae" %in% names(.)) yolksac_larvae else 0,
          early_larvae_days = if ("early_larvae" %in% names(.)) early_larvae else 0,
          late_larvae_days = if ("late_larvae" %in% names(.)) late_larvae else 0,
          juvenile_days = if ("juvenile" %in% names(.)) juvenile else 0,
          adult_days = if ("adult" %in% names(.)) adult else 0
        ) %>%
        as.data.frame()
      
      # 计算每个个体在各个区域的累计次数
      region_summary <- results_total_eachtimestep_all_ages %>%
        group_by(age, No_ind, EEZs_region) %>%
        summarize(region_count = n(), 
                  .groups = 'drop') %>%
        pivot_wider(
          names_from = EEZs_region,
          values_from = region_count,
          values_fill = list(region_count = 0)
        )%>%
        mutate(
          days_in_Ukraine_1 = if ("Ukraine_1" %in% names(.)) Ukraine_1 else 0,
          days_in_Ukraine_2 = if ("Ukraine_2" %in% names(.)) Ukraine_2 else 0,
          days_in_Ukraine_3 = if ("Ukraine_3" %in% names(.)) Ukraine_3 else 0,
          days_in_Romania = if ("Romania" %in% names(.)) Romania else 0,
          days_in_Bulgaria = if ("Bulgaria" %in% names(.)) Bulgaria else 0,
          days_in_Turkey_1 = if ("Turkey_1" %in% names(.)) Turkey_1 else 0,
          days_in_Turkey_2 = if ("Turkey_2" %in% names(.)) Turkey_2 else 0,
          days_in_Turkey_3 = if ("Turkey_3" %in% names(.)) Turkey_3 else 0,
          days_in_Georgia = if ("Georgia" %in% names(.)) Georgia else 0,
          days_in_Russia = if ("Russia" %in% names(.)) Russia else 0
        ) %>%
        as.data.frame()
      
      # 合并阶段数据和区域数据
      summary_total <- stage_summary %>%
        left_join(region_summary, by = c("age", "No_ind"))
      
      # 累加每个阶段的天数、产卵量和每个区域的累计次数
      accumulate_daily_record_SI <- accumulate_daily_record_SI %>%
        mutate(
          egg_days = egg_days + summary_total$egg_days,
          yolksac_larvae_days = yolksac_larvae_days + summary_total$yolksac_larvae_days,
          early_larvae_days = early_larvae_days + summary_total$early_larvae_days,
          late_larvae_days = late_larvae_days + summary_total$late_larvae_days,
          juvenile_days = juvenile_days + summary_total$juvenile_days,
          adult_days = adult_days + summary_total$adult_days,
          
          eggs_total_num_no_worth = eggs_total_num_no_worth + summary_total$eggs_total_num_no_worth,
          eggs_total_num_with_worth = eggs_total_num_with_worth + summary_total$eggs_total_num_with_worth,
          batch_num_no_worth = batch_num_no_worth + summary_total$batch_num_no_worth,
          batch_num_with_worth = batch_num_with_worth + summary_total$batch_num_with_worth,
          
          sum_dis_to_coast = sum_dis_to_coast + summary_total$sum_dis_to_coast,
          
          days_in_Ukraine_1 = days_in_Ukraine_1 + summary_total$days_in_Ukraine_1,
          days_in_Ukraine_2 = days_in_Ukraine_2 + summary_total$days_in_Ukraine_2,
          days_in_Ukraine_3 = days_in_Ukraine_3 + summary_total$days_in_Ukraine_3,
          days_in_Romania = days_in_Romania + summary_total$days_in_Romania,
          days_in_Bulgaria = days_in_Bulgaria + summary_total$days_in_Bulgaria,
          days_in_Turkey_1 = days_in_Turkey_1 + summary_total$days_in_Turkey_1,
          days_in_Turkey_2 = days_in_Turkey_2 + summary_total$days_in_Turkey_2,
          days_in_Turkey_3 = days_in_Turkey_3 + summary_total$days_in_Turkey_3,
          days_in_Georgia = days_in_Georgia + summary_total$days_in_Georgia,
          days_in_Russia = days_in_Russia + summary_total$days_in_Russia
        )%>%
        as.data.frame()
    }
    
  } # end loop for all days in one year
  
  # 【output_6】seasonal abundance and biomass per age-class per year
  saveRDS(egg_seasonal_matrix_abundance,paste0(path_save,output_path_name,current_test_name,'/egg_seasonal_matrix_abundance.rds'))
  saveRDS(yolksac_larvae_seasonal_matrix_abundance,paste0(path_save,output_path_name,current_test_name,'/yolksac_larvae_seasonal_matrix_abundance.rds'))
  
  saveRDS(early_larvae_seasonal_matrix_biomass,paste0(path_save,output_path_name,current_test_name,'/early_larvae_seasonal_matrix_biomass.rds'))
  saveRDS(late_larvae_seasonal_matrix_biomass,paste0(path_save,output_path_name,current_test_name,'/late_larvae_seasonal_matrix_biomass.rds'))
  saveRDS(juvenile_seasonal_matrix_biomass,paste0(path_save,output_path_name,current_test_name,'/juvenile_seasonal_matrix_biomass.rds'))
  saveRDS(adult_seasonal_matrix_biomass,paste0(path_save,output_path_name,current_test_name,'/adult_seasonal_matrix_biomass.rds'))
  
  saveRDS(early_larvae_seasonal_matrix_abundance,paste0(path_save,output_path_name,current_test_name,'/early_larvae_seasonal_matrix_abundance.rds'))
  saveRDS(late_larvae_seasonal_matrix_abundance,paste0(path_save,output_path_name,current_test_name,'/late_larvae_seasonal_matrix_abundance.rds'))
  saveRDS(juvenile_seasonal_matrix_abundance,paste0(path_save,output_path_name,current_test_name,'/juvenile_seasonal_matrix_abundance.rds'))
  saveRDS(adult_seasonal_matrix_abundance,paste0(path_save,output_path_name,current_test_name,'/adult_seasonal_matrix_abundance.rds'))
  
  saveRDS(total_seasonal_matrix_abundance,paste0(path_save,output_path_name,current_test_name,'/total_seasonal_matrix_abundance.rds'))
  saveRDS(total_seasonal_matrix_biomass,paste0(path_save,output_path_name,current_test_name,'/total_seasonal_matrix_biomass.rds'))
  
  # 【output_7】annual length and weight-at-age
  saveRDS(annual_weight_age,paste0(path_save,output_path_name,current_test_name,'/annual_weight_age.rds'))
  saveRDS(annual_length_age,paste0(path_save,output_path_name,current_test_name,'/annual_length_age.rds'))
  saveRDS(annual_Lon_Lat_age,paste0(path_save,output_path_name,current_test_name,'/annual_Lon_Lat_age.rds'))
  
  # 【output_8】mean spawned eggs and batches per female
  saveRDS(eggs_batch_num_age_0,paste0(path_save,output_path_name,current_test_name,'/eggs_batch_num_age_0.rds'))
  saveRDS(eggs_batch_num_age_1,paste0(path_save,output_path_name,current_test_name,'/eggs_batch_num_age_1.rds'))
  saveRDS(eggs_batch_num_age_2,paste0(path_save,output_path_name,current_test_name,'/eggs_batch_num_age_2.rds'))
  saveRDS(eggs_batch_num_age_3,paste0(path_save,output_path_name,current_test_name,'/eggs_batch_num_age_3.rds'))
  saveRDS(eggs_batch_num_age_4,paste0(path_save,output_path_name,current_test_name,'/eggs_batch_num_age_4.rds'))
  
  saveRDS(eggs_total_num_age_0,paste0(path_save,output_path_name,current_test_name,'/eggs_total_num_age_0.rds'))
  saveRDS(eggs_total_num_age_1,paste0(path_save,output_path_name,current_test_name,'/eggs_total_num_age_1.rds'))
  saveRDS(eggs_total_num_age_2,paste0(path_save,output_path_name,current_test_name,'/eggs_total_num_age_2.rds'))
  saveRDS(eggs_total_num_age_3,paste0(path_save,output_path_name,current_test_name,'/eggs_total_num_age_3.rds'))
  saveRDS(eggs_total_num_age_4,paste0(path_save,output_path_name,current_test_name,'/eggs_total_num_age_4.rds'))
  
  # 【output_9】Lagrangian continue recording for the 8 selected individuals
  saveRDS(information_selected_individuals,paste0(path_save,output_path_name,current_test_name,'/information_selected_individuals.rds'))
  
  # 【output_11】monitor of all kinds of mortality (ZFAC and natural mortality) daily but output once a year
  saveRDS(mortality_age_monitor,paste0(path_save,output_path_name,current_test_name,'/High_ZFAC.rds'))
  
  # 【output_12】Annual accumulate record of SIs with duration of life stage, spawning batches, dis to shoreline, and days in regions
  if(i>=6){
    annual_record_SI<-rbind(annual_record_SI, accumulate_daily_record_SI)
    saveRDS(annual_record_SI,paste0(path_save,output_path_name,current_test_name,'/annual_record_SI.rds'))
  }
  
  # 【output_13】overlap with multi-SIs in the same cells
  saveRDS(overlap.occur.all.days,paste0(path_save,output_path_name,current_test_name,'/overlap_cell_SIs_in_DD_mortality.rds'))
  
} # end loop for all years

# 【ouput_4】 all information of each time step for each state variable
saveRDS(results_total_age_0,paste0(path_save,output_path_name,current_test_name,'/results_total_age_0.rds'))
saveRDS(results_total_age_1,paste0(path_save,output_path_name,current_test_name,'/results_total_age_1.rds'))
saveRDS(results_total_age_2,paste0(path_save,output_path_name,current_test_name,'/results_total_age_2.rds'))
saveRDS(results_total_age_3,paste0(path_save,output_path_name,current_test_name,'/results_total_age_3.rds'))
saveRDS(results_total_age_4,paste0(path_save,output_path_name,current_test_name,'/results_total_age_4.rds'))

# 【output_5】determined possible earliest and latest spawning day for each year
saveRDS(all_potential_spawn_day,paste0(path_save,output_path_name,current_test_name,'/all_potential_spawn_day.rds'))

stopCluster(cl)
