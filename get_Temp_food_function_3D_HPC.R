# reconstruct extracting data frame from temperautre and food 20250125
# For 3D

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

get_Temp_3D_cell_center<-function(grid_data, life_stage, depth_level, cell_number){
  # this function is used for a single individual
  
  Temp<-grid_data[cell_number,depth_level]
  Temp <- ifelse(is.na(Temp), 0, Temp)
  
  Temp_spawn<-Temp
  
  Temp_final<-data.frame(Temp=Temp,
                         Temp_spawn=Temp_spawn) # Temp_spawn is similar with Temp
  return(Temp_final)
  
}

get_Temp_3D_interpolation<-function(grid_data, 
                                    life_stage, 
                                    cell_number, 
                                    lon, lat, 
                                    depth, 
                                    depth_level){
  # this function is used for a single individual and cannot be transferred by vertor
  f0_level_Temp <- var_3D_one_level_horizontal_spatial_interpolation(lon = lon, 
                                                                     lat = lat, 
                                                                     var_data = grid_data,
                                                                     depth_level = depth_level)
  
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell_number),'max_depth_level']
  if(depth<=all_depth_levels_t[1] | depth>=all_depth_levels_t[max_depth_level_new] | depth == all_depth_levels_t[depth_level]){ 
    # upper the top cell or below the bottom cell, no need for vertical interpolation
    
    Temp<-f0_level_Temp
    
  }else{
    
    if(depth < all_depth_levels_t[depth_level]){
      next_depth_level<-depth_level-1
    }else{
      next_depth_level<-depth_level+1
    }
    f1_level_Temp <- var_3D_one_level_horizontal_spatial_interpolation(lon = lon, 
                                                                       lat = lat, 
                                                                       var_data = grid_data,
                                                                       depth_level = next_depth_level)
  }
  
  
  Temp <- var_3D_vertical_spatial_interpolation(f0_level_var = f0_level_Temp,
                                                f1_level_var = f1_level_Temp,
                                                cell_number = cell_number,
                                                depth = depth,
                                                depth_level = depth_level)
  
  Temp_spawn <- Temp
  
  Temp_final<-data.frame(Temp=Temp,
                         Temp_spawn=Temp_spawn)
  
  return(Temp_final)
  
}

get_Food_3D_cell_center<-function(MIC.data,
                                  MES.data, 
                                  TL, depth_level, life_stage, cell_number, parms.total){
  # this is used for a single individual
  if(life_stage %in% c('early_larvae','late_larvae','juvenile','adult')){
    MIC.value<-MIC.data[cell_number,depth_level]
    MES.value<-MES.data[cell_number,depth_level]
    MIC.value <- ifelse(is.na(MIC.value), 0, MIC.value)
    MES.value <- ifelse(is.na(MES.value), 0, MES.value)
    
    Food.data<-data.frame(MIC = MIC.value,
                          MES = MES.value)
    # Calculate food allocation
    food.allocation <- food_allocation_2D_3D(total_length = TL)
    
    # Vectorized food calculation based on life_stage
    Food <- functional_response_cell_center_3D(Food_data = Food.data, 
                                               life_stage = life_stage, 
                                               parms.total = parms.total,
                                               V_micro_zoo = food.allocation$V_micro_zoo,
                                               V_meso_zoo = food.allocation$V_meso_zoo)
  }else{
    # this data frame is same as functional_response_cell_center_3D
    Food <- data.frame(Z_micro_zoo = NA_real_,
                       Z_meso_zoo = NA_real_,
                       f_micro_zoo = NA_real_,
                       f_meso_zoo = NA_real_,
                       f = NA_real_)
  }
  
  # return a data frame
  return(Food)
}

get_Food_3D_interpolation<-function(MIC.data,
                                    MES.data,
                                    TL,
                                    life_stage, 
                                    cell_number, 
                                    lon, lat, 
                                    depth, 
                                    depth_level){
  
  # this function is used for a single individual and cannot be transferred by vertor
  
  if(life_stage %in% c('early_larvae','late_larvae','juvenile','adult')){
    f0_level_MIC <- var_3D_one_level_horizontal_spatial_interpolation(lon = lon, 
                                                                      lat = lat, 
                                                                      var_data = MIC.data,
                                                                      depth_level = depth_level)
    f0_level_MES <- var_3D_one_level_horizontal_spatial_interpolation(lon = lon, 
                                                                      lat = lat, 
                                                                      var_data = MES.data,
                                                                      depth_level = depth_level)
    
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell_number),'max_depth_level']
    if(depth<=all_depth_levels_t[1] | depth>=all_depth_levels_t[max_depth_level_new] | depth == all_depth_levels_t[depth_level]){ 
      # upper the top cell or below the bottom cell, no need for vertical interpolation
      
      MIC.value<-f0_level_MIC
      MES.value<-f0_level_MES
      
    }else{
      
      if(depth < all_depth_levels_t[depth_level]){
        next_depth_level<-depth_level-1
      }else{
        next_depth_level<-depth_level+1
      }
      
      f1_level_MIC <- var_3D_one_level_horizontal_spatial_interpolation(lon = lon, 
                                                                        lat = lat, 
                                                                        var_data = MIC.data,
                                                                        depth_level = next_depth_level)
      f1_level_MES <- var_3D_one_level_horizontal_spatial_interpolation(lon = lon, 
                                                                        lat = lat, 
                                                                        var_data = MES.data,
                                                                        depth_level = next_depth_level)
      
      MIC.value <- var_3D_vertical_spatial_interpolation(f0_level_var = f0_level_MIC,
                                                         f1_level_var = f1_level_MIC,
                                                         cell_number = cell_number,
                                                         depth = depth,
                                                         depth_level = depth_level)
      MES.value <- var_3D_vertical_spatial_interpolation(f0_level_var = f0_level_MES,
                                                         f1_level_var = f1_level_MES,
                                                         cell_number = cell_number,
                                                         depth = depth,
                                                         depth_level = depth_level)
    }
    # Because BAMHBI outputs negative values for MIC and MES
    # if they are negative, they will output 0
    MIC.value <- pmax(MIC.value, 0)
    MES.value <- pmax(MES.value, 0)
    
    Food.data<-data.frame(MIC = MIC.value,
                          MES = MES.value)
    
    # Calculate food allocation
    food.allocation <- food_allocation_2D_3D(total_length = TL)
    
    # Vectorized food calculation based on life_stage
    Food <- functional_response_interpolation_2D_3D(Food_data = Food.data, 
                                                    life_stage = life_stage, 
                                                    parms.total = parms.total,
                                                    V_micro_zoo = food.allocation$V_micro_zoo,
                                                    V_meso_zoo = food.allocation$V_meso_zoo)
  }else{ # egg and yolksac larvae
    Food <- data.frame(Z_micro_zoo = NA_real_,
                       Z_meso_zoo = NA_real_,
                       f_micro_zoo = NA_real_,
                       f_meso_zoo = NA_real_,
                       f = NA_real_)
  }
  
  # return a data frame
  return(Food)
  
}





