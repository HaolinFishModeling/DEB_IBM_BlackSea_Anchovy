running_func_age_0_4_2D<-function(batch.of.No.individual,
                                  input.state.var.age.dif,
                                  All.Temp.grid.current.in.day,
                                  All.Food.grid.current.in.day,
                                  All.Temp.grid.next.in.day,
                                  All.Food.grid.next.in.day,
                                  All.uo.grid.current.in.day,
                                  All.vo.grid.current.in.day,
                                  All.uo.grid.next.in.day,
                                  All.vo.grid.next.in.day,
                                  All.sd.grid.current.in.day,
                                  All.sd.grid.next.in.day,
                                  input.func,
                                  input.data.location,
                                  input.stoke.drift,
                                  timestep,
                                  parms.total,
                                  individual_per_core,
                                  current_month
){
  
  # below source should be deleted
  #source('/gpfs/home/acad/ulg-mast/yuhaolin/rtest/functions_to_transfer_ordinate_with_cell_grid.R')
  
  input.state.var<-input.state.var.age.dif[batch.of.No.individual,]           # age_dif
  
  for (n in 1:length(batch.of.No.individual)) {
    aaa <- input.state.var[n,]
    
    # if there is no alive SI, do not extract Temp and Food and do not execute DEB, Kinesis, and PTM
    if(aaa$running_status == 'alive'){
      
      # extract and update Temp and food ==========================
      if(timestep < 1){ # daily without temporal interpolation
        
        if(input.data.location==1){ # extract data in cell center
          Temp.total<-get_Temp_2D_cell_center(grid_data = All.Temp.grid.current.in.day, 
                                              life_stage = aaa$life_stage, 
                                              cell_number = aaa$Cell_num)
          Food.total<-get_Food_2D_cell_center(grid_data = All.Food.grid.current.in.day, 
                                              TL = aaa$TL,
                                              life_stage = aaa$life_stage, 
                                              cell_number = aaa$Cell_num,
                                              parms.total = parms.total)
        }else{ # extract data from spatial interpolation
          Temp.total<-get_Temp_2D_interpolation(grid_data = All.Temp.grid.current.in.day, 
                                                life_stage = aaa$life_stage, 
                                                cell_number = aaa$Cell_num,
                                                lon = aaa$Lon, 
                                                lat = aaa$Lat)
          Food.total<-get_Food_2D_interpolation(grid_data = All.Food.grid.current.in.day, 
                                                TL = aaa$TL,
                                                life_stage = aaa$life_stage, 
                                                cell_number = aaa$Cell_num,
                                                parms.total = parms.total,
                                                lon = aaa$Lon,
                                                lat = aaa$Lat)
        }
        
      }else{ # consider hourly here with temporal interpolation
        if(input.data.location==1){ # extract data in cell center
          Temp.total.current<-get_Temp_2D_cell_center(grid_data = All.Temp.grid.current.in.day, 
                                                      life_stage = aaa$life_stage, 
                                                      cell_number = aaa$Cell_num)
          Food.total.current<-get_Food_2D_cell_center(grid_data = All.Food.grid.current.in.day, 
                                                      TL = aaa$TL,
                                                      life_stage = aaa$life_stage, 
                                                      cell_number = aaa$Cell_num,
                                                      parms.total = parms.total)
          Temp.total.next<-get_Temp_2D_cell_center(grid_data = All.Temp.grid.next.in.day, 
                                                   life_stage = aaa$life_stage, 
                                                   cell_number = aaa$Cell_num)
          Food.total.next<-get_Food_2D_cell_center(grid_data = All.Food.grid.next.in.day, 
                                                   TL = aaa$TL,
                                                   life_stage = aaa$life_stage, 
                                                   cell_number = aaa$Cell_num,
                                                   parms.total = parms.total)
        }else{ # extract data from spatial interpolation
          Temp.total.current<-get_Temp_2D_interpolation(grid_data = All.Temp.grid.current.in.day, 
                                                        life_stage = aaa$life_stage, 
                                                        cell_number = aaa$Cell_num,
                                                        lon = aaa$Lon, 
                                                        lat = aaa$Lat)
          Food.total.current<-get_Food_2D_interpolation(grid_data = All.Food.grid.current.in.day, 
                                                        TL = aaa$TL,
                                                        life_stage = aaa$life_stage, 
                                                        cell_number = aaa$Cell_num,
                                                        parms.total = parms.total,
                                                        lon = aaa$Lon,
                                                        lat = aaa$Lat)
          Temp.total.next<-get_Temp_2D_interpolation(grid_data = All.Temp.grid.next.in.day, 
                                                     life_stage = aaa$life_stage, 
                                                     cell_number = aaa$Cell_num,
                                                     lon = aaa$Lon, 
                                                     lat = aaa$Lat)
          Food.total.next<-get_Food_2D_interpolation(grid_data = All.Food.grid.next.in.day, 
                                                     TL = aaa$TL,
                                                     life_stage = aaa$life_stage, 
                                                     cell_number = aaa$Cell_num,
                                                     parms.total = parms.total,
                                                     lon = aaa$Lon,
                                                     lat = aaa$Lat)
          
        }
        # get the specific hour number of one day
        hourly_step_num<-hour(aaa$Date)+1
        
        # Temperature
        Temp.total<-data.frame(Temp=temporal_interpolation_2D_3D(current_value = Temp.total.current$Temp,
                                                                 next_value = Temp.total.next$Temp,
                                                                 hourly_step_num = hourly_step_num),
                               Temp_spawn=temporal_interpolation_2D_3D(current_value = Temp.total.current$Temp_spawn,
                                                                       next_value = Temp.total.next$Temp_spawn,
                                                                       hourly_step_num = hourly_step_num))
        
        # Food
        if(aaa$life_stage %in% c('early_larvae','late_larvae','juvenile','adult')){
          Food.total<-data.frame(Z_micro_zoo=temporal_interpolation_2D_3D(current_value = Food.total.current$Z_micro_zoo,
                                                                          next_value = Food.total.next$Z_micro_zoo,
                                                                          hourly_step_num = hourly_step_num),
                                 Z_meso_zoo=temporal_interpolation_2D_3D(current_value = Food.total.current$Z_meso_zoo,
                                                                         next_value = Food.total.next$Z_meso_zoo,
                                                                         hourly_step_num = hourly_step_num),
                                 f_micro_zoo = temporal_interpolation_2D_3D(current_value = Food.total.current$f_micro_zoo,
                                                                            next_value = Food.total.next$f_micro_zoo,
                                                                            hourly_step_num = hourly_step_num),
                                 f_meso_zoo = temporal_interpolation_2D_3D(current_value = Food.total.current$f_meso_zoo,
                                                                           next_value = Food.total.next$f_meso_zoo,
                                                                           hourly_step_num = hourly_step_num),
                                 f = temporal_interpolation_2D_3D(current_value = Food.total.current$f,
                                                                  next_value = Food.total.next$f,
                                                                  hourly_step_num = hourly_step_num))
        }else{
          Food.total<-data.frame(Z_micro_zoo=NA_real_,
                                 Z_meso_zoo=NA_real_,
                                 f_micro_zoo = NA_real_,
                                 f_meso_zoo = NA_real_,
                                 f = NA_real_)
        }
      } # End IF hourly
      
      # update state variables for each SI
      input.state.var[n,colnames(Temp.total)]<-Temp.total
      input.state.var[n,colnames(Food.total)]<-Food.total
      
      # call different models to finish different work ==============================
      if(input.func == 1){
        # DEB and Movement and PTM
        
        # 1. call DEB model ----------------
        results_DEB <- DEB_func(
          h = 1 / (24 * timestep),
          state_variables = input.state.var[n,],
          parms = parms.total,
          current_month = current_month
        )
        
        # update state variables of DEB for a single SI
        input.state.var[n,colnames(results_DEB)]<-results_DEB
        
        # 2. call Kinesis model ----------------
        results_movement <- move_func_kinesis(
          delta_t = (60 / timestep) * 60,
          state_variables = input.state.var[n,],
          parms = parms.total
        )
        
        # update state variables of movement for a single SI
        input.state.var[n,colnames(results_movement)]<-results_movement
        
        # 3. call PTM model --------------------
        if(timestep<1){ # daily
          # Daily: only save values at the last hour of the day
          
          results_PTM<-input.state.var[n,]
          for (hour_i in 1:24) {
            aaa<-PTM_function_2D(input.stoke.drift = input.stoke.drift,
                                 state_variables = results_PTM,
                                 dt = 60*60, # one hour
                                 current_uo_data = All.uo.grid.current.in.day,
                                 next_uo_data = All.uo.grid.next.in.day,
                                 current_vo_data = All.vo.grid.current.in.day,
                                 next_vo_data = All.vo.grid.next.in.day,
                                 current_sd_data = All.sd.grid.current.in.day,
                                 next_sd_data = All.sd.grid.next.in.day,
                                 hourly_step_num = hour_i)
            
            results_PTM[1,colnames(aaa)]<-aaa
          }
          input.state.var[n,colnames(results_PTM)]<-results_PTM
          
        }else{ # hourly
          results_PTM<-PTM_function_2D(input.stoke.drift = input.stoke.drift,
                                       state_variables = input.state.var[n,],
                                       dt = 60*60, # one hour
                                       current_uo_data = All.uo.grid.current.in.day,
                                       next_uo_data = All.uo.grid.next.in.day,
                                       current_vo_data = All.vo.grid.current.in.day,
                                       next_vo_data = All.vo.grid.next.in.day,
                                       current_sd_data = All.sd.grid.current.in.day,
                                       next_sd_data = All.sd.grid.next.in.day,
                                       hourly_step_num = hour(input.state.var[n,'Date'])+1)
          
          input.state.var[n,colnames(results_PTM)]<-results_PTM
        }
        
        # End IF of input.fun == 1
        
      }else if(input.func == 2){
        # only DEB function
        # call DEB model ----------------
        results_DEB <- DEB_func(
          h = 1 / (24 * timestep),
          state_variables = input.state.var[n,],
          parms = parms.total,
          current_month = current_month
        )
        
        # update state variables of DEB for a single SI
        input.state.var[n,colnames(results_DEB)]<-results_DEB
        
        # End IF of input.fun == 2
        
      }else if(input.func == 3){
        # only Kinesis function without DEB and PTM
        
        # call Kinesis model --------------
        results_movement <- move_func_kinesis(
          delta_t = (60 / timestep) * 60,
          state_variables = input.state.var[n,],
          parms = parms.total
        )
        
        # update state variables of movement for a single SI
        input.state.var[n,colnames(results_movement)]<-results_movement
        
        
      } # End IF of input.fun == 3
      
    } # End IF of judge whether alive
    
  } # End final loop for each SI
  
  return(input.state.var)
}
