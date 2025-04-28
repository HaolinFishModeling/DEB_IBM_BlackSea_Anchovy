functional_response_cell_center_2D<-function(Food_data, 
                                          life_stage, 
                                          parms.total,
                                          V_micro_zoo,
                                          V_meso_zoo){
  
  if(life_stage %in% 'early_larvae'){
    Z_micro_zoo<-Food_data[,'Food_MIC_0_30m']
    Z_meso_zoo<-Food_data[,'Food_MES_0_30m']
    kf_micro_zoo<-parms.total['kf_early_larvae_micro_zoo']
    kf_meso_zoo<-parms.total['kf_early_larvae_meso_zoo']
  }else if(life_stage %in% 'late_larvae'){
    Z_micro_zoo<-Food_data[,'Food_MIC_0_30m']
    Z_meso_zoo<-Food_data[,'Food_MES_0_30m']
    kf_micro_zoo<-parms.total['kf_late_larvae_micro_zoo']
    kf_meso_zoo<-parms.total['kf_late_larvae_meso_zoo']
  }else if(life_stage %in% 'juvenile'){
    Z_micro_zoo<-Food_data[,'Food_MIC_0_100m']
    Z_meso_zoo<-Food_data[,'Food_MES_0_100m']
    kf_micro_zoo<-parms.total['kf_juvenile_micro_zoo']
    kf_meso_zoo<-parms.total['kf_juvenile_meso_zoo']
  }else if(life_stage %in% 'adult'){
    Z_micro_zoo<-Food_data[,'Food_MIC_0_100m']
    Z_meso_zoo<-Food_data[,'Food_MES_0_100m']
    kf_micro_zoo<-parms.total['kf_adult_micro_zoo']
    kf_meso_zoo<-parms.total['kf_adult_meso_zoo']
  }else{
    Z_micro_zoo<-NA_real_
    Z_meso_zoo<-NA_real_
    kf_micro_zoo<-NA_real_
    kf_meso_zoo<-NA_real_
  }
  
  # Calculate functional response
  f_micro_zoo <- (Z_micro_zoo * V_micro_zoo / kf_micro_zoo) /
    (1 + ((Z_micro_zoo * V_micro_zoo / kf_micro_zoo) + 
            (Z_meso_zoo * V_micro_zoo / kf_meso_zoo)))
  
  f_meso_zoo <- (Z_meso_zoo * V_meso_zoo / kf_meso_zoo) /
    (1 + ((Z_micro_zoo * V_meso_zoo / kf_micro_zoo) + 
            (Z_meso_zoo * V_meso_zoo / kf_meso_zoo)))
  
  f <- f_micro_zoo + f_meso_zoo
  
  # 返回结果
  final.functional.response <- data.frame(
    Z_micro_zoo = Z_micro_zoo,
    Z_meso_zoo = Z_meso_zoo,
    f_micro_zoo = f_micro_zoo,
    f_meso_zoo = f_meso_zoo,
    f = f
  )
  return(final.functional.response)
}


functional_response_cell_center_3D<-function(Food_data, 
                                             life_stage, 
                                             parms.total,
                                             V_micro_zoo,
                                             V_meso_zoo){
  
  Z_micro_zoo<-Food_data[,'MIC']
  Z_meso_zoo<-Food_data[,'MES']
  
  if(life_stage %in% 'early_larvae'){
    kf_micro_zoo<-parms.total['kf_early_larvae_micro_zoo']
    kf_meso_zoo<-parms.total['kf_early_larvae_meso_zoo']
  }else if(life_stage %in% 'late_larvae'){
    kf_micro_zoo<-parms.total['kf_late_larvae_micro_zoo']
    kf_meso_zoo<-parms.total['kf_late_larvae_meso_zoo']
  }else if(life_stage %in% 'juvenile'){
    kf_micro_zoo<-parms.total['kf_juvenile_micro_zoo']
    kf_meso_zoo<-parms.total['kf_juvenile_meso_zoo']
  }else if(life_stage %in% 'adult'){
    kf_micro_zoo<-parms.total['kf_adult_micro_zoo']
    kf_meso_zoo<-parms.total['kf_adult_meso_zoo']
  }else{
    kf_micro_zoo<-NA_real_
    kf_meso_zoo<-NA_real_
  }
  
  # Calculate functional response
  f_micro_zoo <- (Z_micro_zoo * V_micro_zoo / kf_micro_zoo) /
    (1 + ((Z_micro_zoo * V_micro_zoo / kf_micro_zoo) + 
            (Z_meso_zoo * V_micro_zoo / kf_meso_zoo)))
  
  f_meso_zoo <- (Z_meso_zoo * V_meso_zoo / kf_meso_zoo) /
    (1 + ((Z_micro_zoo * V_meso_zoo / kf_micro_zoo) + 
            (Z_meso_zoo * V_meso_zoo / kf_meso_zoo)))
  
  f <- f_micro_zoo + f_meso_zoo
  
  # 返回结果
  final.functional.response <- data.frame(
    Z_micro_zoo = Z_micro_zoo,
    Z_meso_zoo = Z_meso_zoo,
    f_micro_zoo = f_micro_zoo,
    f_meso_zoo = f_meso_zoo,
    f = f
  )
  return(final.functional.response)
}



functional_response_interpolation_2D_3D<-function(Food_data, 
                                                  life_stage, 
                                                  parms.total,
                                                  V_micro_zoo,
                                                  V_meso_zoo){
  
  # Food has already done filter for MIC and MES 0_30 or 0_100m in the interpolation step
  # no matter what 2D or 3D, this is only used for the specific location of the super individual
  # interpolation has already done!
  
  Z_micro_zoo<-Food_data[,'MIC']
  Z_meso_zoo<-Food_data[,'MES']
  
  if(life_stage %in% 'early_larvae'){
    kf_micro_zoo<-parms.total['kf_early_larvae_micro_zoo']
    kf_meso_zoo<-parms.total['kf_early_larvae_meso_zoo']
  }else if(life_stage %in% 'late_larvae'){
    kf_micro_zoo<-parms.total['kf_late_larvae_micro_zoo']
    kf_meso_zoo<-parms.total['kf_late_larvae_meso_zoo']
  }else if(life_stage %in% 'juvenile'){
    kf_micro_zoo<-parms.total['kf_juvenile_micro_zoo']
    kf_meso_zoo<-parms.total['kf_juvenile_meso_zoo']
  }else if(life_stage %in% 'adult'){
    kf_micro_zoo<-parms.total['kf_adult_micro_zoo']
    kf_meso_zoo<-parms.total['kf_adult_meso_zoo']
  }else{
    kf_micro_zoo<-NA_real_
    kf_meso_zoo<-NA_real_
  }
  
  # Calculate functional response
  f_micro_zoo <- (Z_micro_zoo * V_micro_zoo / kf_micro_zoo) /
    (1 + ((Z_micro_zoo * V_micro_zoo / kf_micro_zoo) + 
            (Z_meso_zoo * V_micro_zoo / kf_meso_zoo)))
  
  f_meso_zoo <- (Z_meso_zoo * V_meso_zoo / kf_meso_zoo) /
    (1 + ((Z_micro_zoo * V_meso_zoo / kf_micro_zoo) + 
            (Z_meso_zoo * V_meso_zoo / kf_meso_zoo)))
  
  f <- f_micro_zoo + f_meso_zoo
  
  # 返回结果
  final.functional.response <- data.frame(
    Z_micro_zoo = Z_micro_zoo,
    Z_meso_zoo = Z_meso_zoo,
    f_micro_zoo = f_micro_zoo,
    f_meso_zoo = f_meso_zoo,
    f = f
  )
  return(final.functional.response)
}
