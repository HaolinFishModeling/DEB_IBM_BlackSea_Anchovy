# 20250305 D-D mortality for 3D
natural_and_DD_mortality_3D<-function(input.DD.mortality,
                                      results.each.age){
  
  if(input.DD.mortality==1){  # add D-D mortality only for Late-larvae
    
    Cal_ZFAC_LL<-function(Dt,Db=Db_LL_3D,L){
      ZFAC<--1+2*Dt/Db+0.4*L-0.4*(Dt/Db)*L
      ZFAC <- ifelse(ZFAC < 0, 0, ZFAC)
      return(ZFAC)
    }
    
    Cal_ZFAC_juv<-function(Dt,Db=Db_juv_3D,L){
      ZFAC<-1
      return(ZFAC)
    }
    
    # add mortality and update state variables
    results.each.age<-results.each.age %>%
      mutate(temporal_worth = worth) %>%
      group_by(Cell_num,life_stage) %>%
      mutate(temporal_worth = case_when(
        running_status=='alive' & life_stage %in% c('late_larvae') ~ sum(worth),
        TRUE ~ temporal_worth
      )) %>%
      ungroup() %>%
      mutate(natural_mortality = case_when(
        running_status=='alive' & life_stage == 'egg' ~ parms.total['mort_egg_yolksac'],
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ parms.total['mort_egg_yolksac'],
        running_status=='alive' & life_stage == 'early_larvae' ~ parms.total['mort_early_larvae'],
        running_status=='alive' & life_stage == 'late_larvae' ~ parms.total['mort_late_larvae'],
        running_status=='alive' & life_stage == 'juvenile' ~ parms.total['mort_juvenile'],
        running_status=='alive' & life_stage == 'adult' ~ parms.total['mort_adult'],
        TRUE ~ 0
      ))%>%
      mutate(ZFAC = case_when(
        running_status=='alive' & life_stage == 'egg' ~ 1,
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ 1,
        running_status=='alive' & life_stage == 'early_larvae' ~ 1,
        running_status=='alive' & life_stage == 'late_larvae' ~ Cal_ZFAC_LL(Dt=temporal_worth/(2780*(2780*cos(Lat*pi/180))), L=TL),
        running_status=='alive' & life_stage == 'juvenile' ~ Cal_ZFAC_juv(Dt=temporal_worth/(2780*(2780*cos(Lat*pi/180))), L=TL),
        running_status=='alive' & life_stage == 'adult' ~ 1,
        TRUE ~ 0
      ))%>%
      mutate(final_mortality = case_when(
        running_status=='alive'  ~ natural_mortality*ZFAC,
        TRUE ~ final_mortality
      ))%>%
      mutate(worth = case_when(
        running_status=='alive'  ~ worth*exp(-final_mortality),
        TRUE ~ worth
      ))%>%
      mutate(density = case_when(
        running_status=='alive' & life_stage == 'egg' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'early_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'late_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'juvenile' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'adult' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        TRUE ~ density
      ))%>%
      mutate(norm_density = case_when(
        running_status=='alive' & life_stage == 'egg' ~ NA_real_,
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ NA_real_,
        running_status=='alive' & life_stage == 'early_larvae' ~ NA_real_,
        running_status=='alive' & life_stage == 'late_larvae' ~ density/Db_LL_3D,
        running_status=='alive' & life_stage == 'juvenile' ~ NA_real_,
        running_status=='alive' & life_stage == 'adult' ~ NA_real_,
        TRUE ~ norm_density
      ))%>%
      select(-temporal_worth)%>%
      as.data.frame()
    
    # examine how many SIs were overlapped when they are late_larvae
    aaa<-results.each.age %>%
      filter(life_stage %in% c('late_larvae'))%>%
      group_by(Cell_num) %>%
      mutate(count = n()) %>%
      ungroup()%>%
      filter(count>1)%>%
      as.data.frame()
    
  }else if(input.DD.mortality==2){  # add D-D mortality only for juvenile
    
    Cal_ZFAC_LL<-function(Dt,Db=Db_LL_3D,L){
      ZFAC<-1
      return(ZFAC)
    }
    
    Cal_ZFAC_juv<-function(Dt,Db=Db_juv_3D,L){
      ZFAC<--0.7+1.7*Dt/Db+0.145*L-0.145*(Dt/Db)*L
      ZFAC <- ifelse(ZFAC < 0, 0, ZFAC)
      return(ZFAC)
    }
    
    # add mortality and update state variables
    results.each.age<-results.each.age %>%
      mutate(temporal_worth = worth) %>%
      group_by(Cell_num,life_stage) %>%
      mutate(temporal_worth = case_when(
        running_status=='alive' & life_stage %in% c('juvenile') ~ sum(worth),
        TRUE ~ temporal_worth
      )) %>%
      ungroup() %>%
      mutate(natural_mortality = case_when(
        running_status=='alive' & life_stage == 'egg' ~ parms.total['mort_egg_yolksac'],
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ parms.total['mort_egg_yolksac'],
        running_status=='alive' & life_stage == 'early_larvae' ~ parms.total['mort_early_larvae'],
        running_status=='alive' & life_stage == 'late_larvae' ~ parms.total['mort_late_larvae'],
        running_status=='alive' & life_stage == 'juvenile' ~ parms.total['mort_juvenile'],
        running_status=='alive' & life_stage == 'adult' ~ parms.total['mort_adult'],
        TRUE ~ 0
      ))%>%
      mutate(ZFAC = case_when(
        running_status=='alive' & life_stage == 'egg' ~ 1,
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ 1,
        running_status=='alive' & life_stage == 'early_larvae' ~ 1,
        running_status=='alive' & life_stage == 'late_larvae' ~ Cal_ZFAC_LL(Dt=temporal_worth/(2780*(2780*cos(Lat*pi/180))), L=TL),
        running_status=='alive' & life_stage == 'juvenile' ~ Cal_ZFAC_juv(Dt=temporal_worth/(2780*(2780*cos(Lat*pi/180))), L=TL),
        running_status=='alive' & life_stage == 'adult' ~ 1,
        TRUE ~ 0
      ))%>%
      mutate(final_mortality = case_when(
        running_status=='alive'  ~ natural_mortality*ZFAC,
        TRUE ~ final_mortality
      ))%>%
      mutate(worth = case_when(
        running_status=='alive'  ~ worth*exp(-final_mortality),
        TRUE ~ worth
      ))%>%
      mutate(density = case_when(
        running_status=='alive' & life_stage == 'egg' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'early_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'late_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'juvenile' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'adult' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        TRUE ~ density
      ))%>%
      mutate(norm_density = case_when(
        running_status=='alive' & life_stage == 'egg' ~ NA_real_,
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ NA_real_,
        running_status=='alive' & life_stage == 'early_larvae' ~ NA_real_,
        running_status=='alive' & life_stage == 'late_larvae' ~ NA_real_,
        running_status=='alive' & life_stage == 'juvenile' ~ density/Db_juv_3D,
        running_status=='alive' & life_stage == 'adult' ~ NA_real_,
        TRUE ~ norm_density
      ))%>%
      select(-temporal_worth)%>%
      as.data.frame()
    
    # examine how many SIs were overlapped when they are juvenile
    aaa<-results.each.age %>%
      filter(life_stage %in% c('juvenile'))%>%
      group_by(Cell_num) %>%
      mutate(count = n()) %>%
      ungroup()%>%
      filter(count>1)%>%
      as.data.frame()
    
  }else if(input.DD.mortality==3){  # add D-D mortality both for Late-larvae and juvenile
    
    Cal_ZFAC_LL<-function(Dt,Db=Db_LL_3D,L){
      ZFAC<--1+2*Dt/Db+0.4*L-0.4*(Dt/Db)*L
      ZFAC <- ifelse(ZFAC < 0, 0, ZFAC)
      return(ZFAC)
    }
    
    Cal_ZFAC_juv<-function(Dt,Db=Db_juv_3D,L){
      ZFAC<--0.7+1.7*Dt/Db+0.145*L-0.145*(Dt/Db)*L
      ZFAC <- ifelse(ZFAC < 0, 0, ZFAC)
      return(ZFAC)
    }
    
    # add mortality and update state variables
    results.each.age<-results.each.age %>%
      mutate(temporal_worth = worth) %>%
      group_by(Cell_num,life_stage) %>%
      mutate(temporal_worth = case_when(
        running_status=='alive' & life_stage %in% c('late_larvae','juvenile') ~ sum(worth),
        TRUE ~ temporal_worth
      )) %>%
      ungroup() %>%
      mutate(natural_mortality = case_when(
        running_status=='alive' & life_stage == 'egg' ~ parms.total['mort_egg_yolksac'],
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ parms.total['mort_egg_yolksac'],
        running_status=='alive' & life_stage == 'early_larvae' ~ parms.total['mort_early_larvae'],
        running_status=='alive' & life_stage == 'late_larvae' ~ parms.total['mort_late_larvae'],
        running_status=='alive' & life_stage == 'juvenile' ~ parms.total['mort_juvenile'],
        running_status=='alive' & life_stage == 'adult' ~ parms.total['mort_adult'],
        TRUE ~ 0
      ))%>%
      mutate(ZFAC = case_when(
        running_status=='alive' & life_stage == 'egg' ~ 1,
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ 1,
        running_status=='alive' & life_stage == 'early_larvae' ~ 1,
        running_status=='alive' & life_stage == 'late_larvae' ~ Cal_ZFAC_LL(Dt=temporal_worth/(2780*(2780*cos(Lat*pi/180))), L=TL),
        running_status=='alive' & life_stage == 'juvenile' ~ Cal_ZFAC_juv(Dt=temporal_worth/(2780*(2780*cos(Lat*pi/180))), L=TL),
        running_status=='alive' & life_stage == 'adult' ~ 1,
        TRUE ~ 0
      ))%>%
      mutate(final_mortality = case_when(
        running_status=='alive'  ~ natural_mortality*ZFAC,
        TRUE ~ final_mortality
      ))%>%
      mutate(worth = case_when(
        running_status=='alive'  ~ worth*exp(-final_mortality),
        TRUE ~ worth
      ))%>%
      mutate(density = case_when(
        running_status=='alive' & life_stage == 'egg' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'early_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'late_larvae' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'juvenile' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        running_status=='alive' & life_stage == 'adult' ~ worth/(2780*(2780*cos(Lat*pi/180))),
        TRUE ~ density
      ))%>%
      mutate(norm_density = case_when(
        running_status=='alive' & life_stage == 'egg' ~ NA_real_,
        running_status=='alive' & life_stage == 'yolksac_larvae' ~ NA_real_,
        running_status=='alive' & life_stage == 'early_larvae' ~ NA_real_,
        running_status=='alive' & life_stage == 'late_larvae' ~ density/Db_LL_3D,
        running_status=='alive' & life_stage == 'juvenile' ~ density/Db_juv_3D,
        running_status=='alive' & life_stage == 'adult' ~ NA_real_,
        TRUE ~ norm_density
      ))%>%
      select(-temporal_worth)%>%
      as.data.frame()
    
    # examine how many SIs were overlapped when they are late_larvae and juvenile
    aaa_LL<-results.each.age %>%
      filter(life_stage %in% c('late_larvae'))%>%
      group_by(Cell_num) %>%
      mutate(count = n()) %>%
      ungroup()%>%
      filter(count>1)%>%
      as.data.frame()
    aaa_juv<-results.each.age %>%
      filter(life_stage %in% c('juvenile'))%>%
      group_by(Cell_num) %>%
      mutate(count = n()) %>%
      ungroup()%>%
      filter(count>1)%>%
      as.data.frame()
    aaa<-rbind(aaa_LL,aaa_juv)
  }
  
  if(nrow(aaa)<1){
    aaa <- NULL
  }
  
  return(list(results.each.age,
              aaa))
  
}


