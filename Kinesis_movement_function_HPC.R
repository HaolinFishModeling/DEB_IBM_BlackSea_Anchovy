# Kinesis algorithm movement (temperature and food - two cues) 20250124
# only run for a single SI

# construct random algorithm function------------
move_func_kinesis <- function(delta_t,
                              state_variables,
                              parms) {
  # delta_t, time step for movement model: if loop in one second, delta_t=1; if loop in an hour, delta_t=1*60*60
  # state_variables, include four elements: Lon, Lat, Cell_num, velocity values (Vx, Vy) and Temp for current time point
  # in the future, there will some state variables from DEB model
  # parms, basic parameters for Kinesis movement model
  
  running_status <- as.character(state_variables[1, 'running_status'])
  life_stage <- as.character(state_variables[1, 'life_stage'])
  
  if ((running_status == 'alive') & (life_stage %in% life_stages_in_Kinesis) & !is.na(life_stage)) {
    
    # read in updated state variables
    Lon <- as.numeric(state_variables[1, 'Lon'])
    Lat <- as.numeric(state_variables[1, 'Lat'])
    Cell_num <- as.numeric(state_variables[1, 'Cell_num'])
    Vx <- as.numeric(state_variables[1, 'Vx'])
    Vy <- as.numeric(state_variables[1, 'Vy'])
    V <- as.numeric(state_variables[1, 'V'])                     # in order to determine life stage
    total_length <- as.numeric(state_variables[1, 'TL']) * 10      # unit, transfer cm to mm
  
    # read in Temp and food
    Temp <- as.numeric(state_variables[1, 'Temp'])
    Food <- as.numeric(state_variables[1, 'f']) # f value
    
    # update of state variable
    f_T <- parms['H1'] * exp(-0.5 * (((Temp - parms['T0']) / parms['sigma_T']) ^
                                       2))
    f_P <- parms['H1'] * exp(-0.5 * (((Food - parms['P0']) / parms['sigma_P']) ^
                                       2))
    
    phi <- parms['SSmax'] * total_length / 1000
    mean_epsilon <- sqrt(0.5 * (phi) ^ 2)
    standard_deviation <- 0.5 * phi
    
    signs_x <- sample(c(-1, 1), 1, replace = T)
    signs_y <- sample(c(-1, 1), 1, replace = T)
    
    epsilon_x <- rnorm(1, mean_epsilon, standard_deviation) * signs_x # generating one normal distribution value
    epsilon_y <- rnorm(1, mean_epsilon, standard_deviation) * signs_y # generating one normal distribution value
    
    g_T <- 1 - parms['H2'] * exp(-0.5 * (((Temp - parms['T0']) / parms['sigma_T']) ^
                                           2))
    g_P <- 1 - parms['H2'] * exp(-0.5 * (((Food - parms['P0']) / parms['sigma_P']) ^
                                           2))
    
    T_alpha <- 1 + parms['tau_Temp'] * (1 - f_T)
    P_alpha <- 1 + parms['tau_Food'] * (1 - f_P)
    
    F_alpha <- (T_alpha * f_T + P_alpha * f_P) / (T_alpha + P_alpha)
    G_alpha <- (T_alpha * g_T + P_alpha * g_P) / (T_alpha + P_alpha)
    
    Ix <- Vx * F_alpha
    Iy <- Vy * F_alpha
    
    Rx <- epsilon_x * G_alpha
    Ry <- epsilon_y * G_alpha
    
    new.Vx <- Ix + Rx
    new.Vy <- Iy + Ry
    
    dLon <- new.Vx * delta_t / (1852*60*cos(Lat*pi/180)) # aim to transfer meters into decimal degree of coordinate
    dLat <- new.Vy * delta_t / (1852*60)                 # aim to transfer meters into decimal degree of coordinate
    
    dist.timestep <- sqrt((new.Vx * delta_t) ^ 2 + (new.Vy * delta_t) ^ 2)
    
    new.Lon <- Lon + dLon
    new.Lat <- Lat + dLat
    
    # define current cell number
    new.Cell_num <- trans_lonlat_to_cellnumber(lon = new.Lon, lat = new.Lat)
    
    # identify whether exceed the boundary
    # read in the total cell number
    if (!new.Cell_num %in% cell.number) {
      Theta_index <- sample(seq(1, 359, by = 0.2), replace = FALSE)  # degree
      for (i in Theta_index) {
        aaa <- rotate_clockwise(x = dLon,
                                y = dLat,
                                theta = pi / 180 * i) #transfer from degree to radian
        changed.dLon <- aaa[1]
        changed.dLat <- aaa[2]
        new.Lon <- Lon + changed.dLon
        new.Lat <- Lat + changed.dLat
        
        new.Cell_num <- trans_lonlat_to_cellnumber(lon = new.Lon, lat = new.Lat)
        if (new.Cell_num %in% cell.number) {
          break
        }
      }
    }
    
    # update distance to the coastline
    new.dis_to_coast<-blacksea_points[match(new.Cell_num,blacksea_points$Cell_num),'dis_to_coast']
    new.region<-blacksea_points[match(new.Cell_num,blacksea_points$Cell_num),'EEZs_adj']
    
  } else{
    # silence, that SI has not been allocated to Age-0
    
    Lon <- as.numeric(state_variables[1, 'Lon'])
    Lat <- as.numeric(state_variables[1, 'Lat'])
    Cell_num <- as.numeric(state_variables[1, 'Cell_num'])
    dis_to_coast <- as.numeric(state_variables[1, 'dis_to_coast'])
    
    new.Lon = Lon
    new.Lat = Lat
    new.Cell_num = Cell_num
    new.dis_to_coast = dis_to_coast
    new.region = blacksea_points[match(new.Cell_num,blacksea_points$Cell_num),'EEZs_adj']
    new.Vx = 0
    new.Vy = 0
    dist.timestep = 0
  }
  
  return(
    data.frame(
      Lon = new.Lon,
      Lat = new.Lat,
      Cell_num = new.Cell_num,
      Vx = new.Vx,
      Vy = new.Vy,
      Dist = dist.timestep,
      dis_to_coast = new.dis_to_coast,
      EEZs_region = new.region
    )
  )
}

