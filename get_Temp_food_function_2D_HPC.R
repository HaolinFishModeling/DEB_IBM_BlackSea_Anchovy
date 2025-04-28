# reconstruct extracting data frame from temperautre and food 20250120
# prey width versus fish size and proportional prey categories for anchovy
food_allocation_2D_3D <- function(first_feeding_length = 0.36, # first_feeding length only be changed when adjust model basic parameter
                            total_length) {
  min_w <- 0.04 + 0.005 * (first_feeding_length * 10 - 4)
  max_w <- pmin(0.035 * total_length * 10, 3)
  percent_MIC <- pmin(pmax((0.2 - min_w) / (max_w - min_w), 0),1)  
  percent_MES <- pmin(pmax((max_w - 0.2) / (max_w - min_w), 0),1)
  
  final <- data.frame(
    total_length = total_length,
    min_w = min_w,
    max_w = max_w,
    percent_MIC = percent_MIC,
    percent_MES = percent_MES,
    V_micro_zoo = 1,
    V_meso_zoo = 1
  )
  return(final)
}

get_Temp_2D_cell_center<-function(grid_data, life_stage, cell_number){
  # this function is used for a single individual
  # output Temp for the selected cell_number
  Temp.data<-grid_data[which(grid_data$Cell_number==cell_number),]
  
  # Use ifelse to handle vectorized life_stage
  Temp <- case_when(
    life_stage %in% c('egg', 'yolksac_larvae') ~ Temp.data[,'Temp_0_10m'],  # only for judging new-born eggs location
    life_stage %in% c('early_larvae', 'late_larvae') ~ Temp.data[,'Temp_0_30m'],
    life_stage %in% c('juvenile', 'adult') ~ Temp.data[,'Temp_0_100m'],
    TRUE ~ NA_real_
  )
  
  Temp_spawn <- Temp.data[,'Temp_0_10m']
  
  Temp_final<-data.frame(Temp=Temp,
                         Temp_spawn=Temp_spawn)
  return(Temp_final)
  
}

get_Temp_2D_interpolation<-function(grid_data, life_stage, cell_number, lon, lat){
  # this function is used for a single individual and cannot be transferred by vertor
  
  # The function Temp_2D_spatial_interpolation considered different life stages !!!
  Temp_and_Temp_spawn <- Temp_2D_spatial_interpolation(var_data = grid_data,
                                                       life_stage = life_stage,
                                                       cell_number = cell_number,
                                                       lon = lon,
                                                       lat = lat)
  
  Temp_final<-data.frame(Temp=Temp_and_Temp_spawn[1],
                         Temp_spawn=Temp_and_Temp_spawn[2])
  
  return(Temp_final)
  
}

get_Food_2D_cell_center<-function(grid_data, TL, life_stage, cell_number, parms.total){
  # this is used for a single individual
  
  if(life_stage %in% c('early_larvae','late_larvae','juvenile','adult')){
    Food.data<-grid_data[which(grid_data$Cell_number==cell_number),]
    
    # Calculate food allocation
    food.allocation <- food_allocation_2D_3D(total_length = TL)
    
    # Vectorized food calculation based on life_stage
    Food <- functional_response_cell_center_2D(Food_data = Food.data, 
                                               life_stage = life_stage, 
                                               parms.total = parms.total,
                                               V_micro_zoo = food.allocation$V_micro_zoo,
                                               V_meso_zoo = food.allocation$V_meso_zoo)
  }else{ # egg and yolksac larvae
    # this data frame is same as functional_response_cell_center_2D
    Food <- data.frame(Z_micro_zoo = NA_real_,
                       Z_meso_zoo = NA_real_,
                       f_micro_zoo = NA_real_,
                       f_meso_zoo = NA_real_,
                       f = NA_real_)
  }
  
  return(Food)
}

get_Food_2D_interpolation<-function(grid_data, TL, life_stage, cell_number, parms.total, lon, lat){
  # this is used for a single individual
  
  if(life_stage %in% c('early_larvae','late_larvae','juvenile','adult')){
    Food.data <- Food_2D_spatial_interpolation(var_data = grid_data,
                                               life_stage = life_stage,
                                               cell_number = cell_number,
                                               lon = lon,
                                               lat = lat)
    
    # Calculate food allocation
    food.allocation <- food_allocation_2D_3D(total_length = TL)
    
    # Vectorized food calculation based on life_stage
    Food <- functional_response_interpolation_2D_3D(Food_data = Food.data, 
                                                    life_stage = life_stage, 
                                                    parms.total = parms.total,
                                                    V_micro_zoo = food.allocation$V_micro_zoo,
                                                    V_meso_zoo = food.allocation$V_meso_zoo)
    
  }else{ # egg and yolksac larvae
    # this data frame is same as functional_response_cell_center_2D
    Food <- data.frame(Z_micro_zoo = NA_real_,
                       Z_meso_zoo = NA_real_,
                       f_micro_zoo = NA_real_,
                       f_meso_zoo = NA_real_,
                       f = NA_real_)
  }
  return(Food)
}

get_Temp_new_born_eggs_2D_cell_center<-function(updated_date, cell_number){
  updated_date <- as.Date(updated_date)
  path1<-"/gpfs/projects/acad/bsmfc/Fish_model/preprocessed_data/2D/Temp/"
  aaa <- readRDS(paste0(path1,'Temp_',updated_date,'.rds'))
  
  # generate dataframe to save order of cell_number
  selected_df <- data.frame(Cell_number = cell_number, order = seq_along(cell_number))
  
  # output Temp for the selected cell_number
  Temp.data<-selected_df%>%
    left_join(aaa[,c('Cell_number','Temp_0_10m')],by='Cell_number')%>%   # delete Date to suit Age_0 with silence individual
    arrange(order)%>%
    as.data.frame()
  
  return(Temp.data)
  
}
