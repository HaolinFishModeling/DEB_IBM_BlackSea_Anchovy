# construct DEB model function - 20250124----------
# only for one single SI

DEB_func <- function(h,
                     state_variables,
                     parms,
                     current_month) {
  
  # h, time step for DEB model if loop in one day, h=1; if loop in one hour, h=1/24
  # state_variables, include four elements: E, V, ER, EGAM, values for current time point
  # parms, basic parameters for DEB model
  # current_month, aim to identify whether it is belongs to April - Octorber for spawning season
  
  running_status <- as.character(state_variables[1, 'running_status'])
  life_stage <- as.character(state_variables[1, 'life_stage'])
  date_current <- state_variables[1, 'Date']   # record the spawning date
  
  if((running_status == 'alive') & (life_stage %in% life_stages_in_DEB) & !is.na(life_stage)){
    E <- as.numeric(state_variables[1, 'E'])
    V <- as.numeric(state_variables[1, 'V'])
    ER <- as.numeric(state_variables[1, 'ER'])
    EGAM <- as.numeric(state_variables[1, 'EGAM'])
    gametes_period_each <- as.numeric(state_variables[1, 'gam_dev_period'])
    
    # read in Temp and food
    Temp <- as.numeric(state_variables[1, 'Temp'])
    f    <- as.numeric(state_variables[1, 'f'])
    
    Temp_spawn  <- as.numeric(state_variables[1, 'Temp_spawn'])
    
    Temp.K <- Temp + 273.15                   # transfer to Kelvin temperature
    
    arrhenius_base = exp((parms['T_A'] / parms['T_ref']) - (parms['T_A'] /
                                                              Temp.K))
    
    arrhenius_boundaries_numerator =   1 + exp((parms['T_AL'] / parms['T_ref'])   - (parms['T_AL'] /
                                                                                       parms['T_L'])) + exp((parms['T_AH'] / parms['T_H'])   -  (parms['T_AH'] /
                                                                                                                                                   parms['T_ref']))
    
    arrhenius_boundaries_denominator = 1 + exp((parms['T_AL'] / Temp.K) - (parms['T_AL'] /
                                                                             parms['T_L'])) + exp((parms['T_AH'] / parms['T_H'])  -  (parms['T_AH'] /
                                                                                                                                        Temp.K))
    
    cT <- arrhenius_base * arrhenius_boundaries_numerator / arrhenius_boundaries_denominator
    
    
    #Fluxes---
    p_Am <- parms['p_xm'] * parms['ae'] * cT
    p_M_temp_effect <- parms['p_M']* cT
    
    PA <- p_Am * f * (V ^ (2 / 3))                                          # Assimilation
    PC <- (E / V) /
      (parms['E_G'] + parms['kap'] * (E / V)) *
      ((parms['E_G'] * p_Am * (V ^ (2 / 3))) / parms['E_m'] + p_M_temp_effect * V)     # Catabolic utilization
    PM <- p_M_temp_effect * V                                                   # Somatic maintenance
    
    PG <- max(parms['kap'] * PC - PM, 0)                                        # Growth
    
    PJ <- min(((1 - parms['kap']) / parms['kap']) * p_M_temp_effect * V, ((1 - parms['kap']) / parms['kap']) * p_M_temp_effect * parms['V_adult'])           # Maturity maintenance
    PR <- (1 - parms['kap']) * PC - PJ                                           # reproduction
    
    if (V > parms['V_adult']) {
      #judge puberty
      
      if (Temp_spawn > parms['TR'] & current_month %in% c(6:8)) {
        #judge spawning season
        PR2 <- ER *
          (p_Am / (parms['E_m'] * (V ^ (1 / 3))) + p_M_temp_effect / parms['E_G']) *
          (1 - parms['kap'] * E / (parms['E_G'] * V + parms['kap'] * E))      # Gamete (energy mobilization)
        PM2 <- min(PM - parms['kap'] * PC, PR2)                                     # Emergency somatic maintenance
        PM3 <- max((PM - (parms['kap'] * PC + PM2)) / parms['E_G'], 0)                           # Atresia (gonad resorption)
      } else{
        #not spawning season
        PR2 <- 0
        PM2 <- min(PM - parms['kap'] * PC, PR2)                                     # Emergency somatic maintenance
        PM3 <- max((PM - (parms['kap'] * PC + PM2)) / parms['E_G'], 0)                           # Atresia (gonad resorption)
        
      }
      
    } else{
      #not puberty
      PR2 <- 0
      PM2 <- min(PM - parms['kap'] * PC, PR2)                                     # Emergency somatic maintenance
      PM3 <- max((PM - (parms['kap'] * PC + PM2)) / parms['E_G'], 0)                           # Atresia (gonad resorption)
      
    }
    
    #state variable---
    dE <- h * (PA - PC)
    dV <- h * (PG / parms['E_G'])
    dER <- h * (PR - PR2)
    dEGAM.not.spawn <- h * (PR2 - PM2 - PM3)
    dEGAM.spawn <- h * (PR2 - PM2 - PM3) - 0.6 * EGAM
    
    next.E <- E + dE                                                   # energy reserves
    next.V <- V + dV                                                  # structural body volume (cm^3)
    
    # Reproduction
    if (V > parms['V_adult']) {
      # judge puberty
      
      next.ER <- ER + dER
      
      if (Temp_spawn > parms['TR'] & current_month %in% c(6:8)) {
        # judge spawning season
        
        Ebatch <- parms['Eegg'] * parms['RFbatch'] * V                             # energy per batch
        
        # judge developing period of gemates in daily or hourly
        if(h ==1){    # daily
          time_to_batch <- ceiling(3100 * (Temp_spawn ^ (-2.29)) * 1)       # time interval at each spawning batch, we delete cT, please see the reference in Pethybridge, not the P's paper itself
          gametes_period_each <- gametes_period_each + 1
          
        }else{        # hourly
          time_to_batch <- 3100 * (Temp_spawn ^ (-2.29)) * 24       # time interval at each spawning batch, we delete cT, please see the reference in Pethybridge, not the P's paper itself
          gametes_period_each <- gametes_period_each + 1
        }
        
        if (0.6 * EGAM > Ebatch) {
          # judge whether sufficient energy
          
          if (gametes_period_each >= time_to_batch) {
            # judge whether sufficient developing duration
            
            #Everything is OK to spawn
            next.EGAM <- EGAM + dEGAM.spawn  # spawning occur
            
            # record information
            gametes_period_each <- 0
            spawn_time <- date_current
            spawn_month <- current_month
            spawn_batch_energy <- Ebatch
            eggs <- spawn_batch_energy / 0.15
            
          } else{
            # continue to develop oocytes
            
            next.EGAM <- EGAM + dEGAM.not.spawn
            spawn_time <- 0
            spawn_month <- 0
            spawn_batch_energy <- 0
            eggs <- spawn_batch_energy / 0.15
          }
          
        } else{
          # continue to save energy
          
          next.EGAM <- EGAM + dEGAM.not.spawn
          spawn_time <- 0
          spawn_month <- 0
          spawn_batch_energy <- 0
          eggs <- spawn_batch_energy / 0.15
        }
        
      } else{
        # not reproduction season
        next.EGAM <- EGAM + (-EGAM)
        spawn_time <- 0
        spawn_month <- 0
        spawn_batch_energy <- 0
        eggs <- spawn_batch_energy / 0.15
        gametes_period_each <- 0
      }
      
    } else{
      # not puberty
      
      next.ER <- ER + 0
      next.EGAM <- EGAM + 0
      spawn_time <- 0
      spawn_month <- 0
      spawn_batch_energy <- 0
      eggs <- spawn_batch_energy / 0.15
    }
    
    #output from DEB
    
    if (next.V < parms['V_juvenile']) {
      TL <- next.V ^ (1 / 3) / (0.2 * next.V + parms['del_larvae'])                          # Total length (cm) Pethybridge et al., (2013) P375
    } else{
      TL <- next.V ^ (1 / 3) / (0.2 * parms['V_juvenile'] + parms['del_larvae'])                                     # Total length (cm)
    }
    
    Ev <- parms['uv'] * parms['dv'] * next.V
    Wdw <- Ev / parms['uv'] + (next.E + next.ER) / parms['uE'] + next.EGAM / parms['uG']
    
    Www <- Wdw * parms['gamma_conversion']
    
    WGON <- (0.015 * Ev / parms['uv'] + next.EGAM / parms['uG']) * parms['gamma_conversion']
    
    Ed <- ((Ev + next.E + next.ER + next.EGAM) / Www)
    
    Kful <- 100 * Www / (TL ^ 3)
    
    GSI <- WGON * 100 / Www
    
    growth_weight <- Www - as.numeric(state_variables[1, 'Www'])
    
    growth_length <- TL - as.numeric(state_variables[1, 'TL'])
    
    if (next.V < parms['V_late_larvae']) {
      #at the early larvae stage
      
      life_stage <- 'early_larvae'
      
    } else if (next.V < parms['V_juvenile']) {
      #at the late larvae stage
      
      life_stage <- 'late_larvae'
      
    } else if (next.V < parms['V_adult']) {
      #at the juvenile stage
      
      life_stage <- 'juvenile'
      
    } else{
      #at the adult stage
      
      life_stage <- 'adult'
      
    }
    
  }else{
    # silence that has not been allocated yet
    next.E <- state_variables[1, 'E']
    next.V <- state_variables[1, 'V']
    next.ER <- state_variables[1, 'ER']
    next.EGAM <- state_variables[1, 'EGAM']
    PA <- state_variables[1, 'PA']
    PG <- state_variables[1, 'PG']
    PC <- state_variables[1, 'PC']
    PM <- state_variables[1, 'PM']
    PM2 <- state_variables[1, 'PM2']
    PR2 <- state_variables[1, 'PR2']
    PM3 <- state_variables[1, 'PM3']
    PJ <- state_variables[1, 'PJ']
    PR <- state_variables[1, 'PR']
    TL <- state_variables[1, 'TL']
    Wdw <- state_variables[1, 'Wdw']
    Www <- state_variables[1, 'Www']
    WGON <- state_variables[1, 'WGON']
    Ed <- state_variables[1, 'Ed']
    Kful <- state_variables[1, 'Kful']
    GSI <- state_variables[1, 'GSI']
    gametes_period_each <- state_variables[1, 'gam_dev_period']
    spawn_time <- state_variables[1, 'spawn_time']
    spawn_month <- state_variables[1, 'spawn_month']
    spawn_batch_energy <- state_variables[1, 'spawn_batch_energy']
    eggs <- state_variables[1, 'eggs']
    life_stage <- state_variables[1, 'life_stage']
    growth_weight <- state_variables[1, 'growth_weight']
    growth_length <- state_variables[1, 'growth_length']
  }
  
  return(
    data.frame(
      E = next.E,
      V = next.V,
      ER = next.ER,
      EGAM = next.EGAM,
      PA = PA,
      PG = PG,
      PC = PC,
      PM = PM,
      PM2 = PM2,
      PR2 = PR2,
      PM3 = PM3,
      PJ = PJ,
      PR = PR,
      TL = TL,
      Wdw = Wdw,
      Www = Www,
      WGON = WGON,
      Ed = Ed,
      Kful = Kful,
      GSI = GSI,
      gam_dev_period = gametes_period_each,
      spawn_time = spawn_time,
      spawn_month = spawn_month,
      spawn_batch_energy = spawn_batch_energy,
      eggs = eggs,
      life_stage = life_stage,
      growth_weight = growth_weight,
      growth_length = growth_length
    )
  )
}


