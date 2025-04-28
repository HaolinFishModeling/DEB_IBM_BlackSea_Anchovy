# Particle tracking model (PTM) functions 3D --- 20250125

# 1. basic functions for calculate buoyancy-----------
# calculate dynamic viscosity related by temperature and salinity
cal_miu<-function(Temp,Sal){
  A <- 1.474*(10^-3) + 1.5*(10^-5)*Temp - 3.927*(10^-8)*(Temp^2)
  B <- 1.0734*(10^-5) - 8.5*(10^-8)*Temp + 2.23*(10^-10)*(Temp^2)
  uR <- 1+A*Sal + B*(Sal^2)
  uW <- exp(-3.79418 + 604.129/(139.18+Temp))
  miu <- uW*uR*(10^-3)
  return(miu)
}

# calculate velocity of buoyancy 
cal_buoyancy<-function(g = 9.81,
                       R = 0.855/1000/2, # egg size (transfer to m) refer to (Huret et al., 2016)
                       rho_p = 1000, # kg/m3, mass density of the particle, 【default】it is floating eggs
                       rho_f, # kg/m3, mass density of the fluid
                       miu){
  w_buoyancy<-(2/9)*((rho_p-rho_f)/miu)*g*(R^2)
  return(w_buoyancy)
}

# 2.1 generate functions for interpolation of velocity and temporal (2D)------------------
interpolation_velocity_and_temporal_UVTSR_2D<-function(lon, lat, 
                                                       current_uo_data,
                                                       next_uo_data,
                                                       current_vo_data,
                                                       next_vo_data,
                                                       current_Temp_data,
                                                       next_Temp_data,
                                                       current_Sal_data,
                                                       next_Sal_data,
                                                       current_rho_data,
                                                       next_rho_data,
                                                       depth_level,
                                                       hourly_step_num){ # interpolation only for a single level
  
  # bilinear interpolation and temporal interpolation
  
  # determine cell number
  target.cell.num<-trans_lonlat_to_cellnumber(lat = lat,
                                              lon = lon)
  
  # determine row number and column number
  target.row.con.num<-trans_lonlat_to_rowcol(lat = lat,
                                             lon = lon)
  
  # determine location of cell for u0, v0, w0, Temp0, Sal0, rho0
  u0.lon<-blacksea_points_uo[which(blacksea_points_uo$Cell_num==target.cell.num),c('Longitude')]
  u0.lat<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Latitude')]
  v0.lon<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Longitude')]
  v0.lat<-blacksea_points_vo[which(blacksea_points_vo$Cell_num==target.cell.num),c('Latitude')]
  Temp0.lon<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Longitude')]
  Temp0.lat<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Latitude')]
  Sal0.lon<-Temp0.lon
  Sal0.lat<-Temp0.lat
  rho0.lon<-Temp0.lon
  rho0.lat<-Temp0.lat
  
  # determine location of other u1, u2, u3
  # locating 4 cells for u0, u1, u2, u3
  # u1
  u1.row.con.num<-as.data.frame(target.row.con.num)
  u1.row.con.num$lon_col<-u1.row.con.num$lon_col-1
  u1.cell.num<-trans_rowcol_to_cellnum(nrow_number = u1.row.con.num$lat_row,
                                       col_number = u1.row.con.num$lon_col)
  aaa<-trans_rowcol_to_lonlat(lat_row = u1.row.con.num$lat_row,
                              lon_col = u1.row.con.num$lon_col)
  u1.lon<-aaa[1]+0.0125 
  u1.lat<-aaa[2]
  
  if(lat>= u0.lat){ # on the upper panel
    # u2
    u2.row.con.num<-as.data.frame(target.row.con.num)
    u2.row.con.num$lat_row<-u2.row.con.num$lat_row-1
    u2.cell.num<-trans_rowcol_to_cellnum(nrow_number = u2.row.con.num$lat_row,
                                         col_number = u2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = u2.row.con.num$lat_row,
                                lon_col = u2.row.con.num$lon_col)
    u2.lon<-aaa[1]+0.0125
    u2.lat<-aaa[2]
    
    #u3
    u3.row.con.num<-as.data.frame(target.row.con.num)
    u3.row.con.num$lat_row<-u3.row.con.num$lat_row-1
    u3.row.con.num$lon_col<-u3.row.con.num$lon_col-1
    u3.cell.num<-trans_rowcol_to_cellnum(nrow_number = u3.row.con.num$lat_row,
                                         col_number = u3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = u3.row.con.num$lat_row,
                                lon_col = u3.row.con.num$lon_col)
    u3.lon<-aaa[1]+0.0125
    u3.lat<-aaa[2]
    
  }else{ # on the bottom panel
    # u2
    u2.row.con.num<-as.data.frame(target.row.con.num)
    u2.row.con.num$lat_row<-u2.row.con.num$lat_row+1
    u2.cell.num<-trans_rowcol_to_cellnum(nrow_number = u2.row.con.num$lat_row,
                                         col_number = u2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = u2.row.con.num$lat_row,
                                lon_col = u2.row.con.num$lon_col)
    u2.lon<-aaa[1]+0.0125
    u2.lat<-aaa[2]
    
    #u3
    u3.row.con.num<-as.data.frame(target.row.con.num)
    u3.row.con.num$lat_row<-u3.row.con.num$lat_row+1
    u3.row.con.num$lon_col<-u3.row.con.num$lon_col-1
    u3.cell.num<-trans_rowcol_to_cellnum(nrow_number = u3.row.con.num$lat_row,
                                         col_number = u3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = u3.row.con.num$lat_row,
                                lon_col = u3.row.con.num$lon_col)
    u3.lon<-aaa[1]+0.0125
    u3.lat<-aaa[2]
  }
  
  # locating 4 cells for v0, v1, v2, v3
  # v1
  v1.row.con.num<-as.data.frame(target.row.con.num)
  v1.row.con.num$lat_row<-v1.row.con.num$lat_row+1
  v1.cell.num<-trans_rowcol_to_cellnum(nrow_number = v1.row.con.num$lat_row,
                                       col_number = v1.row.con.num$lon_col)
  aaa<-trans_rowcol_to_lonlat(lat_row = v1.row.con.num$lat_row,
                              lon_col = v1.row.con.num$lon_col)
  v1.lon<-aaa[1]
  v1.lat<-aaa[2]+0.0125
  
  if(lon < v0.lon){ # on the left panel
    # v2
    v2.row.con.num<-as.data.frame(target.row.con.num)
    v2.row.con.num$lon_col<-v2.row.con.num$lon_col-1
    v2.cell.num<-trans_rowcol_to_cellnum(nrow_number = v2.row.con.num$lat_row,
                                         col_number = v2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = v2.row.con.num$lat_row,
                                lon_col = v2.row.con.num$lon_col)
    v2.lon<-aaa[1]
    v2.lat<-aaa[2]+0.0125
    
    #v3
    v3.row.con.num<-as.data.frame(target.row.con.num)
    v3.row.con.num$lat_row<-v3.row.con.num$lat_row+1
    v3.row.con.num$lon_col<-v3.row.con.num$lon_col-1
    v3.cell.num<-trans_rowcol_to_cellnum(nrow_number = v3.row.con.num$lat_row,
                                         col_number = v3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = v3.row.con.num$lat_row,
                                lon_col = v3.row.con.num$lon_col)
    v3.lon<-aaa[1]
    v3.lat<-aaa[2]+0.0125
    
  }else{ # on the right panel
    # v2
    v2.row.con.num<-as.data.frame(target.row.con.num)
    v2.row.con.num$lon_col<-v2.row.con.num$lon_col+1
    v2.cell.num<-trans_rowcol_to_cellnum(nrow_number = v2.row.con.num$lat_row,
                                         col_number = v2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = v2.row.con.num$lat_row,
                                lon_col = v2.row.con.num$lon_col)
    v2.lon<-aaa[1]
    v2.lat<-aaa[2]+0.0125
    
    #v3
    v3.row.con.num<-as.data.frame(target.row.con.num)
    v3.row.con.num$lat_row<-v3.row.con.num$lat_row+1
    v3.row.con.num$lon_col<-v3.row.con.num$lon_col+1
    v3.cell.num<-trans_rowcol_to_cellnum(nrow_number = v3.row.con.num$lat_row,
                                         col_number = v3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = v3.row.con.num$lat_row,
                                lon_col = v3.row.con.num$lon_col)
    v3.lon<-aaa[1]
    v3.lat<-aaa[2]+0.0125
  }
  
  # locating 4 cells for Temp0, Temp1, Temp2, Temp3
  if(lat>= Temp0.lat & lon<= Temp0.lon){ # on the upper left panel
    # Temp1
    Temp1.row.con.num<-as.data.frame(target.row.con.num)
    Temp1.row.con.num$lon_col<-Temp1.row.con.num$lon_col-1
    Temp1.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp1.row.con.num$lat_row,
                                            col_number = Temp1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp1.row.con.num$lat_row,
                                lon_col = Temp1.row.con.num$lon_col)
    Temp1.lon<-aaa[1] 
    Temp1.lat<-aaa[2]
    
    # Temp2
    Temp2.row.con.num<-as.data.frame(target.row.con.num)
    Temp2.row.con.num$lat_row<-Temp2.row.con.num$lat_row-1
    Temp2.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp2.row.con.num$lat_row,
                                            col_number = Temp2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp2.row.con.num$lat_row,
                                lon_col = Temp2.row.con.num$lon_col)
    Temp2.lon<-aaa[1]
    Temp2.lat<-aaa[2]
    
    #Temp3
    Temp3.row.con.num<-as.data.frame(target.row.con.num)
    Temp3.row.con.num$lat_row<-Temp3.row.con.num$lat_row-1
    Temp3.row.con.num$lon_col<-Temp3.row.con.num$lon_col-1
    Temp3.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp3.row.con.num$lat_row,
                                            col_number = Temp3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp3.row.con.num$lat_row,
                                lon_col = Temp3.row.con.num$lon_col)
    Temp3.lon<-aaa[1]
    Temp3.lat<-aaa[2]
    
  }else if(lat>= Temp0.lat & lon > Temp0.lon){ # on the upper right panel
    # Temp1
    Temp1.row.con.num<-as.data.frame(target.row.con.num)
    Temp1.row.con.num$lon_col<-Temp1.row.con.num$lon_col+1
    Temp1.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp1.row.con.num$lat_row,
                                            col_number = Temp1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp1.row.con.num$lat_row,
                                lon_col = Temp1.row.con.num$lon_col)
    Temp1.lon<-aaa[1] 
    Temp1.lat<-aaa[2]
    
    # Temp2
    Temp2.row.con.num<-as.data.frame(target.row.con.num)
    Temp2.row.con.num$lat_row<-Temp2.row.con.num$lat_row-1
    Temp2.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp2.row.con.num$lat_row,
                                            col_number = Temp2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp2.row.con.num$lat_row,
                                lon_col = Temp2.row.con.num$lon_col)
    Temp2.lon<-aaa[1]
    Temp2.lat<-aaa[2]
    
    #Temp3
    Temp3.row.con.num<-as.data.frame(target.row.con.num)
    Temp3.row.con.num$lat_row<-Temp3.row.con.num$lat_row-1
    Temp3.row.con.num$lon_col<-Temp3.row.con.num$lon_col+1
    Temp3.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp3.row.con.num$lat_row,
                                            col_number = Temp3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp3.row.con.num$lat_row,
                                lon_col = Temp3.row.con.num$lon_col)
    Temp3.lon<-aaa[1]
    Temp3.lat<-aaa[2]
  }else if(lat < Temp0.lat & lon <= Temp0.lon){ # on the bottom left panel
    # Temp1
    Temp1.row.con.num<-as.data.frame(target.row.con.num)
    Temp1.row.con.num$lon_col<-Temp1.row.con.num$lon_col-1
    Temp1.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp1.row.con.num$lat_row,
                                            col_number = Temp1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp1.row.con.num$lat_row,
                                lon_col = Temp1.row.con.num$lon_col)
    Temp1.lon<-aaa[1] 
    Temp1.lat<-aaa[2]
    
    # Temp2
    Temp2.row.con.num<-as.data.frame(target.row.con.num)
    Temp2.row.con.num$lat_row<-Temp2.row.con.num$lat_row+1
    Temp2.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp2.row.con.num$lat_row,
                                            col_number = Temp2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp2.row.con.num$lat_row,
                                lon_col = Temp2.row.con.num$lon_col)
    Temp2.lon<-aaa[1]
    Temp2.lat<-aaa[2]
    
    #Temp3
    Temp3.row.con.num<-as.data.frame(target.row.con.num)
    Temp3.row.con.num$lat_row<-Temp3.row.con.num$lat_row+1
    Temp3.row.con.num$lon_col<-Temp3.row.con.num$lon_col-1
    Temp3.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp3.row.con.num$lat_row,
                                            col_number = Temp3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp3.row.con.num$lat_row,
                                lon_col = Temp3.row.con.num$lon_col)
    Temp3.lon<-aaa[1]
    Temp3.lat<-aaa[2]
  }else{ # on the bottom right panel
    # Temp1
    Temp1.row.con.num<-as.data.frame(target.row.con.num)
    Temp1.row.con.num$lon_col<-Temp1.row.con.num$lon_col+1
    Temp1.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp1.row.con.num$lat_row,
                                            col_number = Temp1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp1.row.con.num$lat_row,
                                lon_col = Temp1.row.con.num$lon_col)
    Temp1.lon<-aaa[1] 
    Temp1.lat<-aaa[2]
    
    # Temp2
    Temp2.row.con.num<-as.data.frame(target.row.con.num)
    Temp2.row.con.num$lat_row<-Temp2.row.con.num$lat_row+1
    Temp2.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp2.row.con.num$lat_row,
                                            col_number = Temp2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp2.row.con.num$lat_row,
                                lon_col = Temp2.row.con.num$lon_col)
    Temp2.lon<-aaa[1]
    Temp2.lat<-aaa[2]
    
    #Temp3
    Temp3.row.con.num<-as.data.frame(target.row.con.num)
    Temp3.row.con.num$lat_row<-Temp3.row.con.num$lat_row+1
    Temp3.row.con.num$lon_col<-Temp3.row.con.num$lon_col+1
    Temp3.cell.num<-trans_rowcol_to_cellnum(nrow_number = Temp3.row.con.num$lat_row,
                                            col_number = Temp3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = Temp3.row.con.num$lat_row,
                                lon_col = Temp3.row.con.num$lon_col)
    Temp3.lon<-aaa[1]
    Temp3.lat<-aaa[2]
  }
  
  # Temp1,2,3 and Sal1,2,3 is the same as w1, 2,3
  Sal1.lon<-Temp1.lon
  Sal1.lat<-Temp1.lat
  rho1.lon<-Temp1.lon
  rho1.lat<-Temp1.lat
  
  Sal2.lon<-Temp2.lon
  Sal2.lat<-Temp2.lat
  rho2.lon<-Temp2.lon
  rho2.lat<-Temp2.lat
  
  Sal3.lon<-Temp3.lon
  Sal3.lat<-Temp3.lat
  rho3.lon<-Temp3.lon
  rho3.lat<-Temp3.lat
  
  # extract velocity to u0, u1, u2, u3 and v0, v1, v2, v3 and w0, w1, w2, w3
  # target.cell.num cannot be the cell in the land!!!
  original.u0<-current_uo_data[target.cell.num,depth_level]
  original.v0<-current_vo_data[target.cell.num,depth_level]
  original.Temp0<-current_Temp_data[target.cell.num,depth_level]
  original.Sal0<-current_Sal_data[target.cell.num,depth_level]
  original.rho0<-current_rho_data[target.cell.num,depth_level]
  
  next.u0<-next_uo_data[target.cell.num,depth_level]
  next.v0<-next_vo_data[target.cell.num,depth_level]
  next.Temp0<-next_Temp_data[target.cell.num,depth_level]
  next.Sal0<-next_Sal_data[target.cell.num,depth_level]
  next.rho0<-next_rho_data[target.cell.num,depth_level]
  
  if(u1.cell.num %in% blacksea_points$Cell_num){
    original.u1<-current_uo_data[u1.cell.num,depth_level]
    next.u1<-next_uo_data[u1.cell.num,depth_level]
  }else{
    original.u1<-0
    next.u1<-0
  }
  
  if(u2.cell.num %in% blacksea_points$Cell_num){
    original.u2<-current_uo_data[u2.cell.num,depth_level]
    next.u2<-next_uo_data[u2.cell.num,depth_level]
  }else{
    original.u2<-0
    next.u2<-0
  }
  
  if(u3.cell.num %in% blacksea_points$Cell_num){
    original.u3<-current_uo_data[u3.cell.num,depth_level]
    next.u3<-next_uo_data[u3.cell.num,depth_level]
  }else{
    original.u3<-0
    next.u3<-0
  }
  
  if(v1.cell.num %in% blacksea_points$Cell_num){
    original.v1<-current_vo_data[v1.cell.num,depth_level]
    next.v1<-next_vo_data[v1.cell.num,depth_level]
  }else{
    original.v1<-0
    next.v1<-0
  }
  
  if(v2.cell.num %in% blacksea_points$Cell_num){
    original.v2<-current_vo_data[v2.cell.num,depth_level]
    next.v2<-next_vo_data[v2.cell.num,depth_level]
  }else{
    original.v2<-0
    next.v2<-0
  }
  
  if(v3.cell.num %in% blacksea_points$Cell_num){
    original.v3<-current_vo_data[v3.cell.num,depth_level]
    next.v3<-next_vo_data[v3.cell.num,depth_level]
  }else{
    original.v3<-0
    next.v3<-0
  }
  
  if(Temp1.cell.num %in% blacksea_points$Cell_num){
    original.Temp1<-current_Temp_data[Temp1.cell.num,depth_level]
    next.Temp1<-next_Temp_data[Temp1.cell.num,depth_level]
    original.Sal1<-current_Sal_data[Temp1.cell.num,depth_level]
    next.Sal1<-next_Sal_data[Temp1.cell.num,depth_level]
    original.rho1<-current_rho_data[Temp1.cell.num,depth_level]
    next.rho1<-next_rho_data[Temp1.cell.num,depth_level]
    
  }else{
    original.Temp1<-original.Temp0
    original.Sal1<-original.Sal0
    original.rho1<-original.rho0
    next.Temp1<-next.Temp0
    next.Sal1<-next.Sal0
    next.rho1<-next.rho0
  }
  
  if(Temp2.cell.num %in% blacksea_points$Cell_num){
    original.Temp2<-current_Temp_data[Temp2.cell.num,depth_level]
    next.Temp2<-next_Temp_data[Temp2.cell.num,depth_level]
    original.Sal2<-current_Sal_data[Temp2.cell.num,depth_level]
    next.Sal2<-next_Sal_data[Temp2.cell.num,depth_level]
    original.rho2<-current_rho_data[Temp2.cell.num,depth_level]
    next.rho2<-next_rho_data[Temp2.cell.num,depth_level]
    
  }else{
    original.Temp2<-original.Temp0
    original.Sal2<-original.Sal0
    original.rho2<-original.rho0
    next.Temp2<-next.Temp0
    next.Sal2<-next.Sal0
    next.rho2<-next.rho0
  }
  
  if(Temp3.cell.num %in% blacksea_points$Cell_num){
    original.Temp3<-current_Temp_data[Temp3.cell.num,depth_level]
    next.Temp3<-next_Temp_data[Temp3.cell.num,depth_level]
    original.Sal3<-current_Sal_data[Temp3.cell.num,depth_level]
    next.Sal3<-next_Sal_data[Temp3.cell.num,depth_level]
    original.rho3<-current_rho_data[Temp3.cell.num,depth_level]
    next.rho3<-next_rho_data[Temp3.cell.num,depth_level]
    
  }else{
    original.w3<-0
    next.w3<-0
    
    original.Temp3<-original.Temp0
    original.Sal3<-original.Sal0
    original.rho3<-original.rho0
    next.Temp3<-next.Temp0
    next.Sal3<-next.Sal0
    next.rho3<-next.rho0
  }
  
  # start temporal interpolation
  u0<-seq(from = original.u0, to = next.u0, length.out = 24)[hourly_step_num]
  v0<-seq(from = original.v0, to = next.v0, length.out = 24)[hourly_step_num]
  Temp0<-seq(from = original.Temp0, to = next.Temp0, length.out = 24)[hourly_step_num]
  Sal0<-seq(from = original.Sal0, to = next.Sal0, length.out = 24)[hourly_step_num]
  rho0<-seq(from = original.rho0, to = next.rho0, length.out = 24)[hourly_step_num]
  
  u1<-seq(from = original.u1, to = next.u1, length.out = 24)[hourly_step_num]
  v1<-seq(from = original.v1, to = next.v1, length.out = 24)[hourly_step_num]
  Temp1<-seq(from = original.Temp1, to = next.Temp1, length.out = 24)[hourly_step_num]
  Sal1<-seq(from = original.Sal1, to = next.Sal1, length.out = 24)[hourly_step_num]
  rho1<-seq(from = original.rho1, to = next.rho1, length.out = 24)[hourly_step_num]
  
  u2<-seq(from = original.u2, to = next.u2, length.out = 24)[hourly_step_num]
  v2<-seq(from = original.v2, to = next.v2, length.out = 24)[hourly_step_num]
  Temp2<-seq(from = original.Temp2, to = next.Temp2, length.out = 24)[hourly_step_num]
  Sal2<-seq(from = original.Sal2, to = next.Sal2, length.out = 24)[hourly_step_num]
  rho2<-seq(from = original.rho2, to = next.rho2, length.out = 24)[hourly_step_num]
  
  u3<-seq(from = original.u3, to = next.u3, length.out = 24)[hourly_step_num]
  v3<-seq(from = original.v3, to = next.v3, length.out = 24)[hourly_step_num]
  Temp3<-seq(from = original.Temp3, to = next.Temp3, length.out = 24)[hourly_step_num]
  Sal3<-seq(from = original.Sal3, to = next.Sal3, length.out = 24)[hourly_step_num]
  rho3<-seq(from = original.rho3, to = next.rho3, length.out = 24)[hourly_step_num]
  
  # strat to interpolate
  xsi<-(lon-u0.lon)/(u1.lon-u0.lon)
  target.u_1<-u0 + xsi*(u1-u0)
  xsi<-(lon-u2.lon)/(u3.lon-u2.lon)
  target.u_2<-u2 + xsi*(u3-u2)
  xsi<-(lat-u1.lat)/(u3.lat-u1.lat)
  target.u<-target.u_1+xsi*(target.u_2 - target.u_1)
  
  eta<-(lat-v0.lat)/(v1.lat-v0.lat)
  target.v_1<-v0 + eta*(v1-v0)
  eta<-(lat-v2.lat)/(v3.lat-v2.lat)
  target.v_2<-v2 + eta*(v3-v2)
  eat<-(lon-v1.lon)/(v3.lon-v1.lon)
  target.v<-target.v_1+eta*(target.v_2 - target.v_1)
  
  xsi<-(lon-Temp0.lon)/(Temp1.lon-Temp0.lon)
  target.Temp_1<-Temp0 + xsi*(Temp1-Temp0)
  target.Sal_1<-Sal0 + xsi*(Sal1-Sal0)
  target.rho_1<-rho0 + xsi*(rho1-rho0)
  
  xsi<-(lon-Temp2.lon)/(Temp3.lon-Temp2.lon)
  target.Temp_2<-Temp2 + xsi*(Temp3-Temp2)
  target.Sal_2<-Sal2 + xsi*(Sal3-Sal2)
  target.rho_2<-rho2 + xsi*(rho3-rho2)
  
  xsi<-(lat-Temp1.lat)/(Temp3.lat-Temp1.lat)
  target.Temp<-target.Temp_1+xsi*(target.Temp_2 - target.Temp_1)
  target.Sal<-target.Sal_1+xsi*(target.Sal_2 - target.Sal_1)
  target.rho<-target.rho_1+xsi*(target.rho_2 - target.rho_1)
  
  output<-c(target.u,target.v,target.Temp,target.Sal,target.rho)
  return(output)
}

interpolation_velocity_and_temporal_W_2D<-function(lon, lat,
                                                   current_wo_data,
                                                   next_wo_data,
                                                   depth_level,
                                                   hourly_step_num){ # interpolation only for a single level
  
  # bilinear interpolation and temporal interpolation
  
  # determine cell number
  target.cell.num<-trans_lonlat_to_cellnumber(lat = lat,
                                              lon = lon)
  
  # determine row number and column number
  target.row.con.num<-trans_lonlat_to_rowcol(lat = lat,
                                             lon = lon)
  
  # determine location of cell for u0, v0, w0, Temp0, Sal0, rho0
  w0.lon<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Longitude')]
  w0.lat<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Latitude')]
  
  # locating 4 cells for w0, w1, w2, w3
  if(lat>= w0.lat & lon<= w0.lon){ # on the upper left panel
    # w1
    w1.row.con.num<-as.data.frame(target.row.con.num)
    w1.row.con.num$lon_col<-w1.row.con.num$lon_col-1
    w1.cell.num<-trans_rowcol_to_cellnum(nrow_number = w1.row.con.num$lat_row,
                                         col_number = w1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w1.row.con.num$lat_row,
                                lon_col = w1.row.con.num$lon_col)
    w1.lon<-aaa[1] 
    w1.lat<-aaa[2]
    
    # w2
    w2.row.con.num<-as.data.frame(target.row.con.num)
    w2.row.con.num$lat_row<-w2.row.con.num$lat_row-1
    w2.cell.num<-trans_rowcol_to_cellnum(nrow_number = w2.row.con.num$lat_row,
                                         col_number = w2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w2.row.con.num$lat_row,
                                lon_col = w2.row.con.num$lon_col)
    w2.lon<-aaa[1]
    w2.lat<-aaa[2]
    
    #w3
    w3.row.con.num<-as.data.frame(target.row.con.num)
    w3.row.con.num$lat_row<-w3.row.con.num$lat_row-1
    w3.row.con.num$lon_col<-w3.row.con.num$lon_col-1
    w3.cell.num<-trans_rowcol_to_cellnum(nrow_number = w3.row.con.num$lat_row,
                                         col_number = w3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w3.row.con.num$lat_row,
                                lon_col = w3.row.con.num$lon_col)
    w3.lon<-aaa[1]
    w3.lat<-aaa[2]
    
  }else if(lat>= w0.lat & lon > w0.lon){ # on the upper right panel
    # w1
    w1.row.con.num<-as.data.frame(target.row.con.num)
    w1.row.con.num$lon_col<-w1.row.con.num$lon_col+1
    w1.cell.num<-trans_rowcol_to_cellnum(nrow_number = w1.row.con.num$lat_row,
                                         col_number = w1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w1.row.con.num$lat_row,
                                lon_col = w1.row.con.num$lon_col)
    w1.lon<-aaa[1] 
    w1.lat<-aaa[2]
    
    # w2
    w2.row.con.num<-as.data.frame(target.row.con.num)
    w2.row.con.num$lat_row<-w2.row.con.num$lat_row-1
    w2.cell.num<-trans_rowcol_to_cellnum(nrow_number = w2.row.con.num$lat_row,
                                         col_number = w2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w2.row.con.num$lat_row,
                                lon_col = w2.row.con.num$lon_col)
    w2.lon<-aaa[1]
    w2.lat<-aaa[2]
    
    #w3
    w3.row.con.num<-as.data.frame(target.row.con.num)
    w3.row.con.num$lat_row<-w3.row.con.num$lat_row-1
    w3.row.con.num$lon_col<-w3.row.con.num$lon_col+1
    w3.cell.num<-trans_rowcol_to_cellnum(nrow_number = w3.row.con.num$lat_row,
                                         col_number = w3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w3.row.con.num$lat_row,
                                lon_col = w3.row.con.num$lon_col)
    w3.lon<-aaa[1]
    w3.lat<-aaa[2]
  }else if(lat < w0.lat & lon <= w0.lon){ # on the bottom left panel
    # w1
    w1.row.con.num<-as.data.frame(target.row.con.num)
    w1.row.con.num$lon_col<-w1.row.con.num$lon_col-1
    w1.cell.num<-trans_rowcol_to_cellnum(nrow_number = w1.row.con.num$lat_row,
                                         col_number = w1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w1.row.con.num$lat_row,
                                lon_col = w1.row.con.num$lon_col)
    w1.lon<-aaa[1] 
    w1.lat<-aaa[2]
    
    # w2
    w2.row.con.num<-as.data.frame(target.row.con.num)
    w2.row.con.num$lat_row<-w2.row.con.num$lat_row+1
    w2.cell.num<-trans_rowcol_to_cellnum(nrow_number = w2.row.con.num$lat_row,
                                         col_number = w2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w2.row.con.num$lat_row,
                                lon_col = w2.row.con.num$lon_col)
    w2.lon<-aaa[1]
    w2.lat<-aaa[2]
    
    #w3
    w3.row.con.num<-as.data.frame(target.row.con.num)
    w3.row.con.num$lat_row<-w3.row.con.num$lat_row+1
    w3.row.con.num$lon_col<-w3.row.con.num$lon_col-1
    w3.cell.num<-trans_rowcol_to_cellnum(nrow_number = w3.row.con.num$lat_row,
                                         col_number = w3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w3.row.con.num$lat_row,
                                lon_col = w3.row.con.num$lon_col)
    w3.lon<-aaa[1]
    w3.lat<-aaa[2]
  }else{ # on the bottom right panel
    # w1
    w1.row.con.num<-as.data.frame(target.row.con.num)
    w1.row.con.num$lon_col<-w1.row.con.num$lon_col+1
    w1.cell.num<-trans_rowcol_to_cellnum(nrow_number = w1.row.con.num$lat_row,
                                         col_number = w1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w1.row.con.num$lat_row,
                                lon_col = w1.row.con.num$lon_col)
    w1.lon<-aaa[1] 
    w1.lat<-aaa[2]
    
    # w2
    w2.row.con.num<-as.data.frame(target.row.con.num)
    w2.row.con.num$lat_row<-w2.row.con.num$lat_row+1
    w2.cell.num<-trans_rowcol_to_cellnum(nrow_number = w2.row.con.num$lat_row,
                                         col_number = w2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w2.row.con.num$lat_row,
                                lon_col = w2.row.con.num$lon_col)
    w2.lon<-aaa[1]
    w2.lat<-aaa[2]
    
    #w3
    w3.row.con.num<-as.data.frame(target.row.con.num)
    w3.row.con.num$lat_row<-w3.row.con.num$lat_row+1
    w3.row.con.num$lon_col<-w3.row.con.num$lon_col+1
    w3.cell.num<-trans_rowcol_to_cellnum(nrow_number = w3.row.con.num$lat_row,
                                         col_number = w3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = w3.row.con.num$lat_row,
                                lon_col = w3.row.con.num$lon_col)
    w3.lon<-aaa[1]
    w3.lat<-aaa[2]
  }
  
  # extract velocity to w0, w1, w2, w3
  # target.cell.num cannot be the cell in the land!!!
  original.w0<-current_wo_data[target.cell.num,depth_level]
  
  next.w0<-next_wo_data[target.cell.num,depth_level]
  
  
  if(w1.cell.num %in% blacksea_points$Cell_num){
    original.w1<-current_wo_data[w1.cell.num,depth_level]
    next.w1<-next_wo_data[w1.cell.num,depth_level]
    
  }else{
    original.w1<-original.w0
    next.w1<-next.w0
  }
  
  if(w2.cell.num %in% blacksea_points$Cell_num){
    original.w2<-current_wo_data[w2.cell.num,depth_level]
    next.w2<-next_wo_data[w2.cell.num,depth_level]
    
  }else{
    original.w2<-original.w0
    next.w2<-next.w0
  }
  
  if(w3.cell.num %in% blacksea_points$Cell_num){
    original.w3<-current_wo_data[w3.cell.num,depth_level]
    next.w3<-next_wo_data[w3.cell.num,depth_level]
    
  }else{
    original.w3<-original.w0
    next.w3<-next.w0
  }
  
  # start temporal interpolation
  w0<-seq(from = original.w0, to = next.w0, length.out = 24)[hourly_step_num]
  w1<-seq(from = original.w1, to = next.w1, length.out = 24)[hourly_step_num]
  w2<-seq(from = original.w2, to = next.w2, length.out = 24)[hourly_step_num]
  w3<-seq(from = original.w3, to = next.w3, length.out = 24)[hourly_step_num]
  
  # strat to interpolate
  xsi<-(lon-w0.lon)/(w1.lon-w0.lon)
  target.w_1<-w0 + xsi*(w1-w0)
  
  xsi<-(lon-w2.lon)/(w3.lon-w2.lon)
  target.w_2<-w2 + xsi*(w3-w2)
  
  xsi<-(lat-w1.lat)/(w3.lat-w1.lat)
  target.w<-target.w_1+xsi*(target.w_2 - target.w_1)
  
  output<-c(target.w)
  return(output)
}

# 2.2 define function to interpolate on depth by linear intepolation (1D) ----------
interpolation_velocity_UVTSR_vertical<-function(f0_level_velocity,
                                                f1_level_velocity,
                                                cell_number,
                                                depth,
                                                depth_level){ # including calculate w_buoyancy
  
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell_number),'max_depth_level']
  
  # for variable in T - center of the cell
  if(depth<=all_depth_levels_t[1] | depth>=all_depth_levels_t[max_depth_level_new] | depth == all_depth_levels_t[depth_level]){ # upper the top cell or below the bottom cell
    target.u.vertical<-f0_level_velocity[1]
    target.v.vertical<-f0_level_velocity[2]
    target.Temp.vertical<-f0_level_velocity[3]
    target.Sal.vertical<-f0_level_velocity[4]
    target.rho.vertical<-f0_level_velocity[5]
    
  }else{
    
    if(depth < all_depth_levels_t[depth_level]){
      next_depth_level<-depth_level-1
    }else{
      next_depth_level<-depth_level+1
    }
    f0_level_depth_t<-all_depth_levels_t[depth_level] 
    f1_level_depth_t<-all_depth_levels_t[next_depth_level]
    
    zeta_t<-(depth-f0_level_depth_t)/(f1_level_depth_t - f0_level_depth_t)
    target.u.vertical<-f0_level_velocity[1]+zeta_t*(f1_level_velocity[1] - f0_level_velocity[1])
    target.v.vertical<-f0_level_velocity[2]+zeta_t*(f1_level_velocity[2] - f0_level_velocity[2])
    target.Temp.vertical<-f0_level_velocity[3]+zeta_t*(f1_level_velocity[3] - f0_level_velocity[3])
    target.Sal.vertical<-f0_level_velocity[4]+zeta_t*(f1_level_velocity[4] - f0_level_velocity[4])
    target.rho.vertical<-f0_level_velocity[5]+zeta_t*(f1_level_velocity[5] - f0_level_velocity[5])
  }
  
  return(c(target.u.vertical,
           target.v.vertical,
           target.Temp.vertical,
           target.Sal.vertical,
           target.rho.vertical))
}

interpolation_velocity_W_vertical<-function(f0_level_velocity,
                                            f1_level_velocity,
                                            cell_number,
                                            depth,
                                            depth_level){ # including calculate w_buoyancy
  
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell_number),'max_depth_level']
  
  # for variable in w - upper edge of the cell
  if(depth==all_depth_levels_w[1] | depth>=all_depth_levels_w[max_depth_level_new]){
    target.w.vertical<-f0_level_velocity
  }else{
    f0_level_depth_w<-all_depth_levels_w[depth_level]  # actual depth value is converse with the location of "wo", so, -1
    f1_level_depth_w<-all_depth_levels_w[depth_level+1]
    
    zeta_w<-(depth-f0_level_depth_w)/(f1_level_depth_w - f0_level_depth_w)
    target.w.vertical<-f0_level_velocity+zeta_w*(f1_level_velocity - f0_level_velocity)
    
  }
  
  return(target.w.vertical)
}

update_UVW_buoyancy_specific_3D_location<-function(vertical_UVTSR_velocity,
                                                   vertical_W_velocity){
  # update w velocity with buoyancy
  target.miu<-cal_miu(Temp = vertical_UVTSR_velocity[3],
                      Sal = vertical_UVTSR_velocity[4])
  w_buoyancy<-cal_buoyancy(rho_f = vertical_UVTSR_velocity[5],
                           miu = target.miu)
  target.w<-vertical_W_velocity - w_buoyancy # here, it is a minus because in Stokes' law, <0 = upward; >0 = downward
  
  return(c(vertical_UVTSR_velocity[1], # U
           vertical_UVTSR_velocity[2], # V
           target.w))
}

# 2.3 define function to interpolate stokes drift on 2D---------
interpolation_sd_2D<-function(lon, lat, 
                              current_sd_data,
                              next_sd_data,
                              hourly_step_num){
  
  # bilineary interpolation and temporal interpolation
  
  # determine cell number
  target.cell.num<-trans_lonlat_to_cellnumber(lat = lat,
                                              lon = lon)
  # determine row number and column number
  target.row.con.num<-trans_lonlat_to_rowcol(lat = lat,
                                             lon = lon)
  
  # determine location of cell for sd0
  sd0.lon<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Longitude')]
  sd0.lat<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Latitude')]
  
  # locating 4 cells for sd0, sd1, sd2, sd3
  if(lat>= sd0.lat & lon<= sd0.lon){ # on the upper left panel
    # sd1
    sd1.row.con.num<-as.data.frame(target.row.con.num)
    sd1.row.con.num$lon_col<-sd1.row.con.num$lon_col-1
    sd1.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd1.row.con.num$lat_row,
                                          col_number = sd1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd1.row.con.num$lat_row,
                                lon_col = sd1.row.con.num$lon_col)
    sd1.lon<-aaa[1] 
    sd1.lat<-aaa[2]
    
    # sd2
    sd2.row.con.num<-as.data.frame(target.row.con.num)
    sd2.row.con.num$lat_row<-sd2.row.con.num$lat_row-1
    sd2.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd2.row.con.num$lat_row,
                                          col_number = sd2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd2.row.con.num$lat_row,
                                lon_col = sd2.row.con.num$lon_col)
    sd2.lon<-aaa[1]
    sd2.lat<-aaa[2]
    
    #sd3
    sd3.row.con.num<-as.data.frame(target.row.con.num)
    sd3.row.con.num$lat_row<-sd3.row.con.num$lat_row-1
    sd3.row.con.num$lon_col<-sd3.row.con.num$lon_col-1
    sd3.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd3.row.con.num$lat_row,
                                          col_number = sd3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd3.row.con.num$lat_row,
                                lon_col = sd3.row.con.num$lon_col)
    sd3.lon<-aaa[1]
    sd3.lat<-aaa[2]
    
  }else if(lat>= sd0.lat & lon > sd0.lon){ # on the upper right panel
    # sd1
    sd1.row.con.num<-as.data.frame(target.row.con.num)
    sd1.row.con.num$lon_col<-sd1.row.con.num$lon_col+1
    sd1.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd1.row.con.num$lat_row,
                                          col_number = sd1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd1.row.con.num$lat_row,
                                lon_col = sd1.row.con.num$lon_col)
    sd1.lon<-aaa[1] 
    sd1.lat<-aaa[2]
    
    # sd2
    sd2.row.con.num<-as.data.frame(target.row.con.num)
    sd2.row.con.num$lat_row<-sd2.row.con.num$lat_row-1
    sd2.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd2.row.con.num$lat_row,
                                          col_number = sd2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd2.row.con.num$lat_row,
                                lon_col = sd2.row.con.num$lon_col)
    sd2.lon<-aaa[1]
    sd2.lat<-aaa[2]
    
    #sd3
    sd3.row.con.num<-as.data.frame(target.row.con.num)
    sd3.row.con.num$lat_row<-sd3.row.con.num$lat_row-1
    sd3.row.con.num$lon_col<-sd3.row.con.num$lon_col+1
    sd3.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd3.row.con.num$lat_row,
                                          col_number = sd3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd3.row.con.num$lat_row,
                                lon_col = sd3.row.con.num$lon_col)
    sd3.lon<-aaa[1]
    sd3.lat<-aaa[2]
  }else if(lat < sd0.lat & lon <= sd0.lon){ # on the bottom left panel
    # sd1
    sd1.row.con.num<-as.data.frame(target.row.con.num)
    sd1.row.con.num$lon_col<-sd1.row.con.num$lon_col-1
    sd1.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd1.row.con.num$lat_row,
                                          col_number = sd1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd1.row.con.num$lat_row,
                                lon_col = sd1.row.con.num$lon_col)
    sd1.lon<-aaa[1] 
    sd1.lat<-aaa[2]
    
    # sd2
    sd2.row.con.num<-as.data.frame(target.row.con.num)
    sd2.row.con.num$lat_row<-sd2.row.con.num$lat_row+1
    sd2.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd2.row.con.num$lat_row,
                                          col_number = sd2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd2.row.con.num$lat_row,
                                lon_col = sd2.row.con.num$lon_col)
    sd2.lon<-aaa[1]
    sd2.lat<-aaa[2]
    
    #sd3
    sd3.row.con.num<-as.data.frame(target.row.con.num)
    sd3.row.con.num$lat_row<-sd3.row.con.num$lat_row+1
    sd3.row.con.num$lon_col<-sd3.row.con.num$lon_col-1
    sd3.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd3.row.con.num$lat_row,
                                          col_number = sd3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd3.row.con.num$lat_row,
                                lon_col = sd3.row.con.num$lon_col)
    sd3.lon<-aaa[1]
    sd3.lat<-aaa[2]
  }else{ # on the bottom right panel
    # sd1
    sd1.row.con.num<-as.data.frame(target.row.con.num)
    sd1.row.con.num$lon_col<-sd1.row.con.num$lon_col+1
    sd1.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd1.row.con.num$lat_row,
                                          col_number = sd1.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd1.row.con.num$lat_row,
                                lon_col = sd1.row.con.num$lon_col)
    sd1.lon<-aaa[1] 
    sd1.lat<-aaa[2]
    
    # sd2
    sd2.row.con.num<-as.data.frame(target.row.con.num)
    sd2.row.con.num$lat_row<-sd2.row.con.num$lat_row+1
    sd2.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd2.row.con.num$lat_row,
                                          col_number = sd2.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd2.row.con.num$lat_row,
                                lon_col = sd2.row.con.num$lon_col)
    sd2.lon<-aaa[1]
    sd2.lat<-aaa[2]
    
    #sd3
    sd3.row.con.num<-as.data.frame(target.row.con.num)
    sd3.row.con.num$lat_row<-sd3.row.con.num$lat_row+1
    sd3.row.con.num$lon_col<-sd3.row.con.num$lon_col+1
    sd3.cell.num<-trans_rowcol_to_cellnum(nrow_number = sd3.row.con.num$lat_row,
                                          col_number = sd3.row.con.num$lon_col)
    aaa<-trans_rowcol_to_lonlat(lat_row = sd3.row.con.num$lat_row,
                                lon_col = sd3.row.con.num$lon_col)
    sd3.lon<-aaa[1]
    sd3.lat<-aaa[2]
  }
  
  # extract velocity tosd.u0, sd.u1, sd.u2, sd.u3; and sd.v0, sd.v1, sd.v2, sd.v3
  if(target.cell.num%in%blacksea_points$Cell_num){
    original.sd.u0<-current_sd_data[which(current_sd_data$Cell_number==target.cell.num),'VSDX']
    original.sd.v0<-current_sd_data[which(current_sd_data$Cell_number==target.cell.num),'VSDY']
    original.sd.wmp0<-current_sd_data[which(current_sd_data$Cell_number==target.cell.num),'VTPK']
    original.sd.swh0<-current_sd_data[which(current_sd_data$Cell_number==target.cell.num),'VHM0']
    
    next.sd.u0<-next_sd_data[which(next_sd_data$Cell_number==target.cell.num),'VSDX']
    next.sd.v0<-next_sd_data[which(next_sd_data$Cell_number==target.cell.num),'VSDY']
    next.sd.wmp0<-next_sd_data[which(next_sd_data$Cell_number==target.cell.num),'VTPK']
    next.sd.swh0<-next_sd_data[which(next_sd_data$Cell_number==target.cell.num),'VHM0']
  }else{
    original.sd.u0<-0
    original.sd.v0<-0
    original.sd.wmp0<-0
    original.sd.swh0<-0
    
    next.sd.u0<-0
    next.sd.v0<-0
    next.sd.wmp0<-0
    next.sd.swh0<-0
  }
  
  if(sd1.cell.num %in% blacksea_points$Cell_num){
    original.sd.u1<-current_sd_data[which(current_sd_data$Cell_number==sd1.cell.num),'VSDX']
    next.sd.u1<-next_sd_data[which(next_sd_data$Cell_number==sd1.cell.num),'VSDX']
    original.sd.v1<-current_sd_data[which(current_sd_data$Cell_number==sd1.cell.num),'VSDY']
    next.sd.v1<-next_sd_data[which(next_sd_data$Cell_number==sd1.cell.num),'VSDY']
    original.sd.wmp1<-current_sd_data[which(current_sd_data$Cell_number==sd1.cell.num),'VTPK']
    next.sd.wmp1<-next_sd_data[which(next_sd_data$Cell_number==sd1.cell.num),'VTPK']
    original.sd.swh1<-current_sd_data[which(current_sd_data$Cell_number==sd1.cell.num),'VHM0']
    next.sd.swh1<-next_sd_data[which(next_sd_data$Cell_number==sd1.cell.num),'VHM0']
  }else{
    original.sd.u1<-0
    original.sd.v1<-0
    original.sd.wmp1<-0
    original.sd.swh1<-0
    
    next.sd.u1<-0
    next.sd.v1<-0
    next.sd.wmp1<-0
    next.sd.swh1<-0
  }
  
  if(sd2.cell.num %in% blacksea_points$Cell_num){
    original.sd.u2<-current_sd_data[which(current_sd_data$Cell_number==sd2.cell.num),'VSDX']
    next.sd.u2<-next_sd_data[which(next_sd_data$Cell_number==sd2.cell.num),'VSDX']
    original.sd.v2<-current_sd_data[which(current_sd_data$Cell_number==sd2.cell.num),'VSDY']
    next.sd.v2<-next_sd_data[which(next_sd_data$Cell_number==sd2.cell.num),'VSDY']
    original.sd.wmp2<-current_sd_data[which(current_sd_data$Cell_number==sd2.cell.num),'VTPK']
    next.sd.wmp2<-next_sd_data[which(next_sd_data$Cell_number==sd2.cell.num),'VTPK']
    original.sd.swh2<-current_sd_data[which(current_sd_data$Cell_number==sd2.cell.num),'VHM0']
    next.sd.swh2<-next_sd_data[which(next_sd_data$Cell_number==sd2.cell.num),'VHM0']
  }else{
    original.sd.u2<-0
    original.sd.v2<-0
    original.sd.wmp2<-0
    original.sd.swh2<-0
    next.sd.u2<-0
    next.sd.v2<-0
    next.sd.wmp2<-0
    next.sd.swh2<-0
  }
  
  if(sd3.cell.num %in% blacksea_points$Cell_num){
    original.sd.u3<-current_sd_data[which(current_sd_data$Cell_number==sd3.cell.num),'VSDX']
    next.sd.u3<-next_sd_data[which(next_sd_data$Cell_number==sd3.cell.num),'VSDX']
    original.sd.v3<-current_sd_data[which(current_sd_data$Cell_number==sd3.cell.num),'VSDY']
    next.sd.v3<-next_sd_data[which(next_sd_data$Cell_number==sd3.cell.num),'VSDY']
    original.sd.wmp3<-current_sd_data[which(current_sd_data$Cell_number==sd3.cell.num),'VTPK']
    next.sd.wmp3<-next_sd_data[which(next_sd_data$Cell_number==sd3.cell.num),'VTPK']
    original.sd.swh3<-current_sd_data[which(current_sd_data$Cell_number==sd3.cell.num),'VHM0']
    next.sd.swh3<-next_sd_data[which(next_sd_data$Cell_number==sd3.cell.num),'VHM0']
  }else{
    original.sd.u3<-0
    original.sd.v3<-0
    original.sd.wmp3<-0
    original.sd.swh3<-0
    next.sd.u3<-0
    next.sd.v3<-0
    next.sd.wmp3<-0
    next.sd.swh3<-0
  }
  
  # start temporal interpolation
  sd.u0<-seq(from = original.sd.u0, to = next.sd.u0, length.out = 24)[hourly_step_num]
  sd.v0<-seq(from = original.sd.v0, to = next.sd.v0, length.out = 24)[hourly_step_num]
  sd.u1<-seq(from = original.sd.u1, to = next.sd.u1, length.out = 24)[hourly_step_num]
  sd.v1<-seq(from = original.sd.v1, to = next.sd.v1, length.out = 24)[hourly_step_num]
  sd.u2<-seq(from = original.sd.u2, to = next.sd.u2, length.out = 24)[hourly_step_num]
  sd.v2<-seq(from = original.sd.v2, to = next.sd.v2, length.out = 24)[hourly_step_num]
  sd.u3<-seq(from = original.sd.u3, to = next.sd.u3, length.out = 24)[hourly_step_num]
  sd.v3<-seq(from = original.sd.v3, to = next.sd.v3, length.out = 24)[hourly_step_num]
  
  sd.wmp0<-seq(from = original.sd.wmp0, to = next.sd.wmp0, length.out = 24)[hourly_step_num]
  sd.swh0<-seq(from = original.sd.swh0, to = next.sd.swh0, length.out = 24)[hourly_step_num]
  sd.wmp1<-seq(from = original.sd.wmp1, to = next.sd.wmp1, length.out = 24)[hourly_step_num]
  sd.swh1<-seq(from = original.sd.swh1, to = next.sd.swh1, length.out = 24)[hourly_step_num]
  sd.wmp2<-seq(from = original.sd.wmp2, to = next.sd.wmp2, length.out = 24)[hourly_step_num]
  sd.swh2<-seq(from = original.sd.swh2, to = next.sd.swh2, length.out = 24)[hourly_step_num]
  sd.wmp3<-seq(from = original.sd.wmp3, to = next.sd.wmp3, length.out = 24)[hourly_step_num]
  sd.swh3<-seq(from = original.sd.swh3, to = next.sd.swh3, length.out = 24)[hourly_step_num]
  
  # strat to interpolate
  # sd.u and sd.v
  xsi<-(lon-sd0.lon)/(sd1.lon-sd0.lon)
  target.sd.u_1<-sd.u0 + xsi*(sd.u1-sd.u0)
  target.sd.v_1<-sd.v0 + xsi*(sd.v1-sd.v0)
  target.sd.wmp_1<-sd.wmp0 + xsi*(sd.wmp1-sd.wmp0)
  target.sd.swh_1<-sd.swh0 + xsi*(sd.swh1-sd.swh0)
  
  xsi<-(lon-sd2.lon)/(sd3.lon-sd2.lon)
  target.sd.u_2<-sd.u2 + xsi*(sd.u3-sd.u2)
  target.sd.v_2<-sd.v2 + xsi*(sd.v3-sd.v2)
  target.sd.wmp_2<-sd.wmp2 + xsi*(sd.wmp3-sd.wmp2)
  target.sd.swh_2<-sd.swh2 + xsi*(sd.swh3-sd.swh2)
  
  xsi<-(lat-sd1.lat)/(sd3.lat-sd1.lat)
  target.sd.u<-target.sd.u_1+xsi*(target.sd.u_2 - target.sd.u_1)
  target.sd.v<-target.sd.v_1+xsi*(target.sd.v_2 - target.sd.v_1)
  target.sd.wmp<-target.sd.wmp_1+xsi*(target.sd.wmp_2 - target.sd.wmp_1)
  target.sd.swh<-target.sd.swh_1+xsi*(target.sd.swh_2 - target.sd.swh_1)
  
  output<-c(target.sd.u,target.sd.v,target.sd.wmp,target.sd.swh)
  return(output)
}

# 2.4 define function to calculate velocity of stokes drift at the desired depth
# !!!!!!!!!!!! depth is positive !!!!!!!!!!!!!!!!!!!!!
cal_sd_at_desired_depth<-function(us, vs, wmp, swh, depth){
  stokes_surface_speed<-sqrt(us^2+vs^2)
  if(stokes_surface_speed == 0){
    output<-c(0,0)
    return(output)
  }else{
    fm02<-1/wmp
    total_transport <- (2*pi/16)*fm02*(swh^2)
    k <- stokes_surface_speed/(2*total_transport)
    stokesspeed<-stokes_surface_speed*exp(-2*k*depth)
    usz <- stokesspeed*us/stokes_surface_speed
    vsz <- stokesspeed*vs/stokes_surface_speed
    output<-c(usz, vsz)
    return(output)
  }
}

# 3.1 define function of updating location (first 2D, this is the same as 2D)------------
update_dlon_dlat_3D<-function(lon, lat, dt, 
                              velocity, 
                              depth_level){
  dlon<-velocity[1] * dt / (1852*60*cos(lat*pi/180))
  dlat<-velocity[2] * dt / (1852*60)
  
  # because we used the hourly time step, we need to check whether it will exceed the boundary by RK4
  new_lon<-lon + dlon
  new_lat<-lat + dlat
  new_cell_num<-trans_lonlat_to_cellnumber(lon = new_lon,
                                           lat = new_lat)
  
  if (!new_cell_num %in% blacksea_points$Cell_num) {
    dlon <- 0
    dlat <- 0
  }else{
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==new_cell_num),'max_depth_level']
    if(depth_level>max_depth_level_new){
      dlon <- 0
      dlat <- 0
    }
  }
  
  return(c(dlon, dlat))
}

# 3.2 define function of updating d_depth (1D) ---------
update_d_depth<-function(lon, lat, dlon, dlat, 
                         depth,
                         dt, 
                         w_vertical){
  
  # define the current location of particle
  target.cell.num<-trans_lonlat_to_cellnumber(lon = lon+dlon,
                                              lat = lat+dlat)
  
  d_depth<-w_vertical * dt
  
  # because we used the hourly time step, we need to check whether it will exceed the boundary by RK4
  new_depth<-depth + d_depth*(-1)
  
  new_depth_level<-define_depth_level(depth = new_depth)
  
  if (new_depth < 0) {
    d_depth <- depth # let the new_depth == 0
  }else{ 
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==target.cell.num),'max_depth_level']
    if(new_depth_level>max_depth_level_new){ # let the particle keep upper the sea floor
      d_depth <- depth-all_depth_levels_w[max_depth_level_new]+0.1 # still retain the particle near the sea floo with 0.1m
    }
  }
  return(d_depth)
}

# 4.1 define function of Runge-Kutta 4 (3D) considering wave effect-------------
RK4_PTM_temporal_interpolation_sd_3D<-function(lon,lat, cell_number,
                                               dt,
                                               current_uo_data,
                                               next_uo_data,
                                               current_vo_data,
                                               next_vo_data,
                                               current_wo_data,
                                               next_wo_data,
                                               current_Temp_data,
                                               next_Temp_data,
                                               current_Sal_data,
                                               next_Sal_data,
                                               current_rho_data,
                                               next_rho_data,
                                               current_sd_data,
                                               next_sd_data,
                                               depth,
                                               depth_level,
                                               hourly_step_num){
  
  # calculate k1, k2, k3, k4
  # k1
  # for variable in T - center of the cell
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell_number),'max_depth_level']
  if(depth<=all_depth_levels_t[1] | depth>=all_depth_levels_t[max_depth_level_new] | depth == all_depth_levels_t[depth_level]){ # upper the top cell or below the bottom cell
    next_depth_level_UVTSR<-depth_level
    
  }else{
    
    if(depth < all_depth_levels_t[depth_level]){
      next_depth_level_UVTSR<-depth_level-1
    }else{
      next_depth_level_UVTSR<-depth_level+1
    }
  }
  
  # for variable in w - upper edge of the cell
  if(depth==all_depth_levels_w[1] | depth>=all_depth_levels_w[max_depth_level_new]){
    next_depth_level_W<-depth_level
  }else{
    next_depth_level_W<-depth_level+1
  }
  
  UVTSR_f0_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon,lat, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = depth_level,
                                                               hourly_step_num)
  W_f0_level<-interpolation_velocity_and_temporal_W_2D(lon,lat, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = depth_level,
                                                       hourly_step_num)
  
  UVTSR_f1_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon,lat, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = next_depth_level_UVTSR,
                                                               hourly_step_num)
  W_f1_level<-interpolation_velocity_and_temporal_W_2D(lon,lat, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = next_depth_level_W,
                                                       hourly_step_num)
  
  UVTSR_interpolation_vertical<-interpolation_velocity_UVTSR_vertical(f0_level_velocity = UVTSR_f0_level,
                                                                      f1_level_velocity = UVTSR_f1_level,
                                                                      cell_number,
                                                                      depth,
                                                                      depth_level)  # current depth level
  W_interpolation_vertical<-interpolation_velocity_W_vertical(f0_level_velocity = W_f0_level,
                                                              f1_level_velocity = W_f1_level,
                                                              cell_number,
                                                              depth,
                                                              depth_level)  # current depth level
  
  velocity1<-update_UVW_buoyancy_specific_3D_location(vertical_UVTSR_velocity = UVTSR_interpolation_vertical,
                                                      vertical_W_velocity = W_interpolation_vertical)
  
  interpolated_sd<-interpolation_sd_2D(lon, lat, 
                                       current_sd_data,
                                       next_sd_data,
                                       hourly_step_num)
  sd_velocity1<-cal_sd_at_desired_depth(us = interpolated_sd[1], 
                                        vs = interpolated_sd[2], 
                                        wmp = interpolated_sd[3], 
                                        swh = interpolated_sd[4], 
                                        depth)
  
  # update u and v velocity at the desired depth
  velocity1[1]<-velocity1[1]+sd_velocity1[1]
  velocity1[2]<-velocity1[2]+sd_velocity1[2]
  
  k1_UV<-update_dlon_dlat_3D(lon, lat, dt, 
                             velocity = velocity1,
                             depth_level = depth_level)
  k1_W<-update_d_depth(lon, lat, 
                       dlon = k1_UV[1], 
                       dlat = k1_UV[2], 
                       depth,
                       dt,
                       w_vertical = velocity1[3])
  
  # k2
  new_depth_k2 <- depth + k1_W/2*(-1)  # this is only for calculating k2
  new_depth_level_k2 <- define_depth_level(depth = new_depth_k2)
  
  # judge whether it will enter the land cell
  lon_k2<-lon+k1_UV[1]/2
  lat_k2<-lat+k1_UV[2]/2
  cell.number_k2<-trans_lonlat_to_cellnumber(lon = lon_k2,
                                             lat = lat_k2)
  
  if(!cell.number_k2 %in% blacksea_points$Cell_num){
    lon_k2<-lon
    lat_k2<-lat
    cell.number_k2<-cell_number
  }else{
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell.number_k2),'max_depth_level']
    if(new_depth_level_k2>max_depth_level_new){
      lon_k2<-lon
      lat_k2<-lat
      cell.number_k2<-cell_number
    }
  }
  
  # for variable in T - center of the cell
  if(new_depth_k2<=all_depth_levels_t[1] | new_depth_k2>=all_depth_levels_t[max_depth_level_new] | new_depth_k2 == all_depth_levels_t[new_depth_level_k2]){ # upper the top cell or below the bottom cell
    next_depth_level_UVTSR_k2<-new_depth_level_k2
    
  }else{
    
    if(new_depth_k2 < all_depth_levels_t[new_depth_level_k2]){
      next_depth_level_UVTSR_k2<-new_depth_level_k2-1
    }else{
      next_depth_level_UVTSR_k2<-new_depth_level_k2+1
    }
  }
  
  # for variable in w - upper edge of the cell
  if(new_depth_k2==all_depth_levels_w[1] | new_depth_k2>=all_depth_levels_w[max_depth_level_new]){
    next_depth_level_W_k2<-new_depth_level_k2
  }else{
    next_depth_level_W_k2<-new_depth_level_k2+1
  }
  
  UVTSR_f0_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k2,
                                                               lat = lat_k2, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = new_depth_level_k2,
                                                               hourly_step_num)
  W_f0_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k2,
                                                       lat = lat_k2, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = new_depth_level_k2,
                                                       hourly_step_num)
  
  UVTSR_f1_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k2,
                                                               lat = lat_k2, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = next_depth_level_UVTSR_k2,
                                                               hourly_step_num)
  W_f1_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k2,
                                                       lat = lat_k2, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = next_depth_level_W_k2,
                                                       hourly_step_num)
  
  UVTSR_interpolation_vertical<-interpolation_velocity_UVTSR_vertical(f0_level_velocity = UVTSR_f0_level,
                                                                      f1_level_velocity = UVTSR_f1_level,
                                                                      cell_number = cell.number_k2,
                                                                      depth = new_depth_k2,
                                                                      depth_level = new_depth_level_k2)  # current depth level
  
  W_interpolation_vertical<-interpolation_velocity_W_vertical(f0_level_velocity = W_f0_level,
                                                              f1_level_velocity = W_f1_level,
                                                              cell_number = cell.number_k2,
                                                              depth = new_depth_k2,
                                                              depth_level = new_depth_level_k2)  # current depth level
  
  velocity2<-update_UVW_buoyancy_specific_3D_location(vertical_UVTSR_velocity = UVTSR_interpolation_vertical,
                                                      vertical_W_velocity = W_interpolation_vertical)
  interpolated_sd<-interpolation_sd_2D(lon = lon_k2, 
                                       lat = lat_k2, 
                                       current_sd_data,
                                       next_sd_data,
                                       hourly_step_num)
  sd_velocity2<-cal_sd_at_desired_depth(us = interpolated_sd[1], 
                                        vs = interpolated_sd[2], 
                                        wmp = interpolated_sd[3], 
                                        swh = interpolated_sd[4], 
                                        depth = new_depth_k2)
  
  # update u and v velocity at the desired depth
  velocity2[1]<-velocity2[1]+sd_velocity2[1]
  velocity2[2]<-velocity2[2]+sd_velocity2[2]
  
  k2_UV<-update_dlon_dlat_3D(lon = lon_k2, 
                             lat = lat_k2, 
                             dt, 
                             velocity = velocity2,
                             depth_level = new_depth_level_k2)
  
  k2_W<-update_d_depth(lon = lon_k2,
                       lat = lat_k2, 
                       dlon = k2_UV[1], 
                       dlat = k2_UV[2], 
                       depth = new_depth_k2,
                       dt, 
                       w_vertical = velocity2[3])
  
  # k3
  new_depth_k3 <- depth + k2_W/2*(-1)  # this is only for calculating k3
  new_depth_level_k3 <- define_depth_level(depth = new_depth_k3)
  
  # judge whether it will enter the land cell
  lon_k3<-lon+k2_UV[1]/2
  lat_k3<-lat+k2_UV[2]/2
  cell.number_k3<-trans_lonlat_to_cellnumber(lon = lon_k3,
                                             lat = lat_k3)
  
  if(!cell.number_k3 %in% blacksea_points$Cell_num){
    lon_k3<-lon
    lat_k3<-lat
    cell.number_k3<-cell_number
  }else{
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell.number_k3),'max_depth_level']
    if(new_depth_level_k3>max_depth_level_new){
      lon_k3<-lon
      lat_k3<-lat
      cell.number_k3<-cell_number
    }
  }
  
  # for variable in T - center of the cell
  if(new_depth_k3<=all_depth_levels_t[1] | new_depth_k3>=all_depth_levels_t[max_depth_level_new] | new_depth_k3 == all_depth_levels_t[new_depth_level_k3]){ # upper the top cell or below the bottom cell
    next_depth_level_UVTSR_k3<-new_depth_level_k3
    
  }else{
    
    if(new_depth_k3 < all_depth_levels_t[new_depth_level_k3]){
      next_depth_level_UVTSR_k3<-new_depth_level_k3-1
    }else{
      next_depth_level_UVTSR_k3<-new_depth_level_k3+1
    }
  }
  
  # for variable in w - upper edge of the cell
  if(new_depth_k3==all_depth_levels_w[1] | new_depth_k3>=all_depth_levels_w[max_depth_level_new]){
    next_depth_level_W_k3<-new_depth_level_k3
  }else{
    next_depth_level_W_k3<-new_depth_level_k3+1
  }
  
  UVTSR_f0_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k3,
                                                               lat = lat_k3, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = new_depth_level_k3,
                                                               hourly_step_num)
  W_f0_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k3,
                                                       lat = lat_k3, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = new_depth_level_k3,
                                                       hourly_step_num)
  
  UVTSR_f1_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k3,
                                                               lat = lat_k3, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = next_depth_level_UVTSR_k3,
                                                               hourly_step_num)
  W_f1_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k3,
                                                       lat = lat_k3, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = next_depth_level_W_k3,
                                                       hourly_step_num)
  
  UVTSR_interpolation_vertical<-interpolation_velocity_UVTSR_vertical(f0_level_velocity = UVTSR_f0_level,
                                                                      f1_level_velocity = UVTSR_f1_level,
                                                                      cell_number = cell.number_k3,
                                                                      depth = new_depth_k3,
                                                                      depth_level = new_depth_level_k3)  # current depth level
  
  W_interpolation_vertical<-interpolation_velocity_W_vertical(f0_level_velocity = W_f0_level,
                                                              f1_level_velocity = W_f1_level,
                                                              cell_number = cell.number_k3,
                                                              depth = new_depth_k3,
                                                              depth_level = new_depth_level_k3)  # current depth level
  
  velocity3<-update_UVW_buoyancy_specific_3D_location(vertical_UVTSR_velocity = UVTSR_interpolation_vertical,
                                                      vertical_W_velocity = W_interpolation_vertical)
  
  interpolated_sd<-interpolation_sd_2D(lon = lon_k3, 
                                       lat = lat_k3, 
                                       current_sd_data,
                                       next_sd_data,
                                       hourly_step_num)
  sd_velocity3<-cal_sd_at_desired_depth(us = interpolated_sd[1], 
                                        vs = interpolated_sd[2], 
                                        wmp = interpolated_sd[3], 
                                        swh = interpolated_sd[4], 
                                        depth = new_depth_k3)
  
  # update u and v velocity at the desired depth
  velocity3[1]<-velocity3[1]+sd_velocity3[1]
  velocity3[2]<-velocity3[2]+sd_velocity3[2]
  
  k3_UV<-update_dlon_dlat_3D(lon = lon_k3, 
                             lat = lat_k3, 
                             dt, 
                             velocity = velocity3,
                             depth_level = new_depth_level_k3)
  
  k3_W<-update_d_depth(lon = lon_k3,
                       lat = lat_k3, 
                       dlon = k3_UV[1], 
                       dlat = k3_UV[2], 
                       depth = new_depth_k3,
                       dt, 
                       w_vertical = velocity3[3])
  
  # k4
  new_depth_k4 <- depth + k3_W*(-1)  # this is only for calculating k4
  new_depth_level_k4 <- define_depth_level(depth = new_depth_k4)
  
  # judge whether it will enter the land cell
  lon_k4<-lon+k3_UV[1]
  lat_k4<-lat+k3_UV[2]
  cell.number_k4<-trans_lonlat_to_cellnumber(lon = lon_k4,
                                             lat = lat_k4)
  
  if(!cell.number_k4 %in% blacksea_points$Cell_num){
    lon_k4<-lon
    lat_k4<-lat
    cell.number_k4<-cell_number
  }else{
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell.number_k4),'max_depth_level']
    if(new_depth_level_k4>max_depth_level_new){
      lon_k4<-lon
      lat_k4<-lat
      cell.number_k4<-cell_number
    }
  }
  
  # for variable in T - center of the cell
  if(new_depth_k4<=all_depth_levels_t[1] | new_depth_k4>=all_depth_levels_t[max_depth_level_new] | new_depth_k4 == all_depth_levels_t[new_depth_level_k4]){ # upper the top cell or below the bottom cell
    next_depth_level_UVTSR_k4<-new_depth_level_k4
    
  }else{
    
    if(new_depth_k4 < all_depth_levels_t[new_depth_level_k4]){
      next_depth_level_UVTSR_k4<-new_depth_level_k4-1
    }else{
      next_depth_level_UVTSR_k4<-new_depth_level_k4+1
    }
  }
  
  # for variable in w - upper edge of the cell
  if(new_depth_k4==all_depth_levels_w[1] | new_depth_k4>=all_depth_levels_w[max_depth_level_new]){
    next_depth_level_W_k4<-new_depth_level_k4
  }else{
    next_depth_level_W_k4<-new_depth_level_k4+1
  }
  
  UVTSR_f0_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k4,
                                                               lat = lat_k4, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = new_depth_level_k4,
                                                               hourly_step_num)
  W_f0_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k4,
                                                       lat = lat_k4, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = new_depth_level_k4,
                                                       hourly_step_num)
  
  UVTSR_f1_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k4,
                                                               lat = lat_k4, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = next_depth_level_UVTSR_k4,
                                                               hourly_step_num)
  W_f1_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k4,
                                                       lat = lat_k4, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = next_depth_level_W_k4,
                                                       hourly_step_num)
  
  UVTSR_interpolation_vertical<-interpolation_velocity_UVTSR_vertical(f0_level_velocity = UVTSR_f0_level,
                                                                      f1_level_velocity = UVTSR_f1_level,
                                                                      cell_number = cell.number_k4,
                                                                      depth = new_depth_k4,
                                                                      depth_level = new_depth_level_k4)  # current depth level
  
  W_interpolation_vertical<-interpolation_velocity_W_vertical(f0_level_velocity = W_f0_level,
                                                              f1_level_velocity = W_f1_level,
                                                              cell_number = cell.number_k4,
                                                              depth = new_depth_k4,
                                                              depth_level = new_depth_level_k4)  # current depth level
  
  velocity4<-update_UVW_buoyancy_specific_3D_location(vertical_UVTSR_velocity = UVTSR_interpolation_vertical,
                                                      vertical_W_velocity = W_interpolation_vertical)
  
  interpolated_sd<-interpolation_sd_2D(lon = lon_k4, 
                                       lat = lat_k4, 
                                       current_sd_data,
                                       next_sd_data,
                                       hourly_step_num)
  sd_velocity4<-cal_sd_at_desired_depth(us = interpolated_sd[1], 
                                        vs = interpolated_sd[2], 
                                        wmp = interpolated_sd[3], 
                                        swh = interpolated_sd[4], 
                                        depth = new_depth_k4)
  
  # update u and v velocity at the desired depth
  velocity4[1]<-velocity4[1]+sd_velocity4[1]
  velocity4[2]<-velocity4[2]+sd_velocity4[2]
  
  k4_UV<-update_dlon_dlat_3D(lon = lon_k4, 
                             lat = lat_k4, 
                             dt, 
                             velocity = velocity4,
                             depth_level = new_depth_level_k4)
  
  k4_W<-update_d_depth(lon = lon_k4,
                       lat = lat_k4, 
                       dlon = k4_UV[1], 
                       dlat = k4_UV[2], 
                       depth = new_depth_k4,
                       dt, 
                       w_vertical = velocity4[3])
  
  
  # calculate the final new location of latitude and longitude
  dlon <- (k1_UV[1] + 2*k2_UV[1] + 2*k3_UV[1] + k4_UV[1])/6
  dlat <- (k1_UV[2] + 2*k2_UV[2] + 2*k3_UV[2] + k4_UV[2])/6
  d_depth <- (k1_W + 2*k2_W + 2*k3_W + k4_W)/6
  
  new_lon <- lon + dlon
  new_lat <- lat + dlat
  new_depth <- depth + d_depth*(-1)
  
  new_cell_num<-trans_lonlat_to_cellnumber(lon = new_lon,
                                           lat = new_lat)
  new_depth_level <- define_depth_level(depth = new_depth)
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==new_cell_num),'max_depth_level']
  
  # when eggs stunk on the seaside or the bottom floor side - retain the last timestep
  if (!new_cell_num %in% blacksea_points$Cell_num) {
    dlon <- 0
    dlat <- 0
    new_lon <- lon
    new_lat <- lat
    new_cell_num <- cell_number
    
  }else{
    # if particle would enter the right or left stone on the seafloor
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==new_cell_num),'max_depth_level']
    if(new_depth_level>max_depth_level_new){
      dlon <- 0
      dlat <- 0
      new_lon <- lon
      new_lat <- lat
      new_cell_num <- cell_number
      
      # let the particle keep upper the sea floor
      new_depth <- all_depth_levels_w[max_depth_level_new+1]-0.1 # still retain the particle near the sea floo with 0.1m
      d_depth <- depth - new_depth 
      new_depth_level <- max_depth_level_new
    }
  }
  
  # when eggs float above the sea surface
  if (new_depth < 0) {
    new_depth<-0
    d_depth <- depth - new_depth
    new_depth_level <- 1
  }
  
  return(c(dlon, dlat, d_depth, new_lon,new_lat,new_cell_num,new_depth,new_depth_level))
  
}

# 4.2 define function of Runge-Kutta 4 (3D) without considering wave effect-------------
RK4_PTM_temporal_interpolation_3D<-function(lon,lat, cell_number,
                                            dt,
                                            current_uo_data,
                                            next_uo_data,
                                            current_vo_data,
                                            next_vo_data,
                                            current_wo_data,
                                            next_wo_data,
                                            current_Temp_data,
                                            next_Temp_data,
                                            current_Sal_data,
                                            next_Sal_data,
                                            current_rho_data,
                                            next_rho_data,
                                            depth,
                                            depth_level,
                                            hourly_step_num){
  # calculate k1, k2, k3, k4
  # k1
  # for variable in T - center of the cell
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell_number),'max_depth_level']
  if(depth<=all_depth_levels_t[1] | depth>=all_depth_levels_t[max_depth_level_new] | depth == all_depth_levels_t[depth_level]){ # upper the top cell or below the bottom cell
    next_depth_level_UVTSR<-depth_level
    
  }else{
    
    if(depth < all_depth_levels_t[depth_level]){
      next_depth_level_UVTSR<-depth_level-1
    }else{
      next_depth_level_UVTSR<-depth_level+1
    }
  }
  
  # for variable in w - upper edge of the cell
  if(depth==all_depth_levels_w[1] | depth>=all_depth_levels_w[max_depth_level_new]){
    next_depth_level_W<-depth_level
  }else{
    next_depth_level_W<-depth_level+1
  }
  
  UVTSR_f0_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon,lat, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = depth_level,
                                                               hourly_step_num)
  W_f0_level<-interpolation_velocity_and_temporal_W_2D(lon,lat, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = depth_level,
                                                       hourly_step_num)
  
  UVTSR_f1_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon,lat, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = next_depth_level_UVTSR,
                                                               hourly_step_num)
  W_f1_level<-interpolation_velocity_and_temporal_W_2D(lon,lat, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = next_depth_level_W,
                                                       hourly_step_num)
  
  UVTSR_interpolation_vertical<-interpolation_velocity_UVTSR_vertical(f0_level_velocity = UVTSR_f0_level,
                                                                      f1_level_velocity = UVTSR_f1_level,
                                                                      cell_number,
                                                                      depth,
                                                                      depth_level)  # current depth level
  
  W_interpolation_vertical<-interpolation_velocity_W_vertical(f0_level_velocity = W_f0_level,
                                                              f1_level_velocity = W_f1_level,
                                                              cell_number,
                                                              depth,
                                                              depth_level)  # current depth level
  
  velocity1<-update_UVW_buoyancy_specific_3D_location(vertical_UVTSR_velocity = UVTSR_interpolation_vertical,
                                                      vertical_W_velocity = W_interpolation_vertical)
  
  k1_UV<-update_dlon_dlat_3D(lon, lat, dt, 
                             velocity = velocity1,
                             depth_level = depth_level)
  k1_W<-update_d_depth(lon, lat, 
                       dlon = k1_UV[1], 
                       dlat = k1_UV[2], 
                       depth = depth,
                       dt,
                       w_vertical = velocity1[3])
  
  # k2
  new_depth_k2 <- depth + k1_W/2*(-1)  # this is only for calculating k2
  new_depth_level_k2 <- define_depth_level(depth = new_depth_k2)
  
  # judge whether it will enter the land cell
  lon_k2<-lon+k1_UV[1]/2
  lat_k2<-lat+k1_UV[2]/2
  cell.number_k2<-trans_lonlat_to_cellnumber(lon = lon_k2,
                                             lat = lat_k2)
  
  if(!cell.number_k2 %in% blacksea_points$Cell_num){
    lon_k2<-lon
    lat_k2<-lat
    cell.number_k2<-cell_number
  }else{
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell.number_k2),'max_depth_level']
    if(new_depth_level_k2>max_depth_level_new){
      lon_k2<-lon
      lat_k2<-lat
      cell.number_k2<-cell_number
    }
  }
  
  # for variable in T - center of the cell
  if(new_depth_k2<=all_depth_levels_t[1] | new_depth_k2>=all_depth_levels_t[max_depth_level_new] | new_depth_k2 == all_depth_levels_t[new_depth_level_k2]){ # upper the top cell or below the bottom cell
    next_depth_level_UVTSR_k2<-new_depth_level_k2
    
  }else{
    
    if(new_depth_k2 < all_depth_levels_t[new_depth_level_k2]){
      next_depth_level_UVTSR_k2<-new_depth_level_k2-1
    }else{
      next_depth_level_UVTSR_k2<-new_depth_level_k2+1
    }
  }
  
  # for variable in w - upper edge of the cell
  if(new_depth_k2==all_depth_levels_w[1] | new_depth_k2>=all_depth_levels_w[max_depth_level_new]){
    next_depth_level_W_k2<-new_depth_level_k2
  }else{
    next_depth_level_W_k2<-new_depth_level_k2+1
  }
  
  UVTSR_f0_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k2,
                                                               lat = lat_k2, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = new_depth_level_k2,
                                                               hourly_step_num)
  W_f0_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k2,
                                                       lat = lat_k2, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = new_depth_level_k2,
                                                       hourly_step_num)
  
  UVTSR_f1_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k2,
                                                               lat = lat_k2, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = next_depth_level_UVTSR_k2,
                                                               hourly_step_num)
  W_f1_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k2,
                                                       lat = lat_k2, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = next_depth_level_W_k2,
                                                       hourly_step_num)
  
  UVTSR_interpolation_vertical<-interpolation_velocity_UVTSR_vertical(f0_level_velocity = UVTSR_f0_level,
                                                                      f1_level_velocity = UVTSR_f1_level,
                                                                      cell_number = cell.number_k2,
                                                                      depth = new_depth_k2,
                                                                      depth_level = new_depth_level_k2)  # current depth level
  
  W_interpolation_vertical<-interpolation_velocity_W_vertical(f0_level_velocity = W_f0_level,
                                                              f1_level_velocity = W_f1_level,
                                                              cell_number = cell.number_k2,
                                                              depth = new_depth_k2,
                                                              depth_level = new_depth_level_k2)  # current depth level
  
  velocity2<-update_UVW_buoyancy_specific_3D_location(vertical_UVTSR_velocity = UVTSR_interpolation_vertical,
                                                      vertical_W_velocity = W_interpolation_vertical)
  
  k2_UV<-update_dlon_dlat_3D(lon = lon_k2, 
                             lat = lat_k2, 
                             dt, 
                             velocity = velocity2,
                             depth_level = new_depth_level_k2)
  
  k2_W<-update_d_depth(lon = lon_k2,
                       lat = lat_k2, 
                       dlon = k2_UV[1], 
                       dlat = k2_UV[2], 
                       depth = new_depth_k2,
                       dt, 
                       w_vertical = velocity2[3])
  
  # k3
  new_depth_k3 <- depth + k2_W/2*(-1)  # this is only for calculating k3
  new_depth_level_k3 <- define_depth_level(depth = new_depth_k3)
  
  # judge whether it will enter the land cell
  lon_k3<-lon+k2_UV[1]/2
  lat_k3<-lat+k2_UV[2]/2
  cell.number_k3<-trans_lonlat_to_cellnumber(lon = lon_k3,
                                             lat = lat_k3)
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell.number_k3),'max_depth_level']
  
  if(!cell.number_k3 %in% blacksea_points$Cell_num){
    lon_k3<-lon
    lat_k3<-lat
    cell.number_k3<-cell_number
  }else{
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell.number_k3),'max_depth_level']
    if(new_depth_level_k3>max_depth_level_new){
      lon_k3<-lon
      lat_k3<-lat
      cell.number_k3<-cell_number
    }
  }
  
  # for variable in T - center of the cell
  if(new_depth_k3<=all_depth_levels_t[1] | new_depth_k3>=all_depth_levels_t[max_depth_level_new] | new_depth_k3 == all_depth_levels_t[new_depth_level_k3]){ # upper the top cell or below the bottom cell
    next_depth_level_UVTSR_k3<-new_depth_level_k3
    
  }else{
    
    if(new_depth_k3 < all_depth_levels_t[new_depth_level_k3]){
      next_depth_level_UVTSR_k3<-new_depth_level_k3-1
    }else{
      next_depth_level_UVTSR_k3<-new_depth_level_k3+1
    }
  }
  
  # for variable in w - upper edge of the cell
  if(new_depth_k3==all_depth_levels_w[1] | new_depth_k3>=all_depth_levels_w[max_depth_level_new]){
    next_depth_level_W_k3<-new_depth_level_k3
  }else{
    next_depth_level_W_k3<-new_depth_level_k3+1
  }
  
  UVTSR_f0_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k3,
                                                               lat = lat_k3, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = new_depth_level_k3,
                                                               hourly_step_num)
  W_f0_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k3,
                                                       lat = lat_k3, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = new_depth_level_k3,
                                                       hourly_step_num)
  
  UVTSR_f1_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k3,
                                                               lat = lat_k3, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = next_depth_level_UVTSR_k3,
                                                               hourly_step_num)
  W_f1_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k3,
                                                       lat = lat_k3, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = next_depth_level_W_k3,
                                                       hourly_step_num)
  
  UVTSR_interpolation_vertical<-interpolation_velocity_UVTSR_vertical(f0_level_velocity = UVTSR_f0_level,
                                                                      f1_level_velocity = UVTSR_f1_level,
                                                                      cell_number = cell.number_k3,
                                                                      depth = new_depth_k3,
                                                                      depth_level = new_depth_level_k3)  # current depth level
  
  W_interpolation_vertical<-interpolation_velocity_W_vertical(f0_level_velocity = W_f0_level,
                                                              f1_level_velocity = W_f1_level,
                                                              cell_number = cell.number_k3,
                                                              depth = new_depth_k3,
                                                              depth_level = new_depth_level_k3)  # current depth level
  
  velocity3<-update_UVW_buoyancy_specific_3D_location(vertical_UVTSR_velocity = UVTSR_interpolation_vertical,
                                                      vertical_W_velocity = W_interpolation_vertical)
  
  k3_UV<-update_dlon_dlat_3D(lon = lon_k3, 
                             lat = lat_k3, 
                             dt, 
                             velocity = velocity3,
                             depth_level = new_depth_level_k3)
  
  k3_W<-update_d_depth(lon = lon_k3,
                       lat = lat_k3, 
                       dlon = k3_UV[1], 
                       dlat = k3_UV[2], 
                       depth = new_depth_k3,
                       dt, 
                       w_vertical = velocity3[3])
  
  # k4
  new_depth_k4 <- depth + k3_W*(-1)  # this is only for calculating k4
  new_depth_level_k4 <- define_depth_level(depth = new_depth_k4)
  
  # judge whether it will enter the land cell
  lon_k4<-lon+k3_UV[1]
  lat_k4<-lat+k3_UV[2]
  cell.number_k4<-trans_lonlat_to_cellnumber(lon = lon_k4,
                                             lat = lat_k4)
  
  if(!cell.number_k4 %in% blacksea_points$Cell_num){
    lon_k4<-lon
    lat_k4<-lat
    cell.number_k4<-cell_number
  }else{
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell.number_k4),'max_depth_level']
    if(new_depth_level_k4>max_depth_level_new){
      lon_k4<-lon
      lat_k4<-lat
      cell.number_k4<-cell_number
    }
  }
  
  # for variable in T - center of the cell
  if(new_depth_k4<=all_depth_levels_t[1] | new_depth_k4>=all_depth_levels_t[max_depth_level_new] | new_depth_k4 == all_depth_levels_t[new_depth_level_k4]){ # upper the top cell or below the bottom cell
    next_depth_level_UVTSR_k4<-new_depth_level_k4
    
  }else{
    
    if(new_depth_k4 < all_depth_levels_t[new_depth_level_k4]){
      next_depth_level_UVTSR_k4<-new_depth_level_k4-1
    }else{
      next_depth_level_UVTSR_k4<-new_depth_level_k4+1
    }
  }
  
  # for variable in w - upper edge of the cell
  if(new_depth_k4==all_depth_levels_w[1] | new_depth_k4>=all_depth_levels_w[max_depth_level_new]){
    next_depth_level_W_k4<-new_depth_level_k4
  }else{
    next_depth_level_W_k4<-new_depth_level_k4+1
  }
  
  UVTSR_f0_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k4,
                                                               lat = lat_k4, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = new_depth_level_k4,
                                                               hourly_step_num)
  W_f0_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k4,
                                                       lat = lat_k4, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = new_depth_level_k4,
                                                       hourly_step_num)
  
  UVTSR_f1_level<-interpolation_velocity_and_temporal_UVTSR_2D(lon = lon_k4,
                                                               lat = lat_k4, 
                                                               current_uo_data,
                                                               next_uo_data,
                                                               current_vo_data,
                                                               next_vo_data,
                                                               current_Temp_data,
                                                               next_Temp_data,
                                                               current_Sal_data,
                                                               next_Sal_data,
                                                               current_rho_data,
                                                               next_rho_data,
                                                               depth_level = next_depth_level_UVTSR_k4,
                                                               hourly_step_num)
  W_f1_level<-interpolation_velocity_and_temporal_W_2D(lon = lon_k4,
                                                       lat = lat_k4, 
                                                       current_wo_data,
                                                       next_wo_data,
                                                       depth_level = next_depth_level_W_k4,
                                                       hourly_step_num)
  
  UVTSR_interpolation_vertical<-interpolation_velocity_UVTSR_vertical(f0_level_velocity = UVTSR_f0_level,
                                                                      f1_level_velocity = UVTSR_f1_level,
                                                                      cell_number = cell.number_k4,
                                                                      depth = new_depth_k4,
                                                                      depth_level = new_depth_level_k4)  # current depth level
  
  W_interpolation_vertical<-interpolation_velocity_W_vertical(f0_level_velocity = W_f0_level,
                                                              f1_level_velocity = W_f1_level,
                                                              cell_number = cell.number_k4,
                                                              depth = new_depth_k4,
                                                              depth_level = new_depth_level_k4)  # current depth level
  
  velocity4<-update_UVW_buoyancy_specific_3D_location(vertical_UVTSR_velocity = UVTSR_interpolation_vertical,
                                                      vertical_W_velocity = W_interpolation_vertical)
  
  
  k4_UV<-update_dlon_dlat_3D(lon = lon_k4, 
                             lat = lat_k4, 
                             dt, 
                             velocity = velocity4,
                             depth_level = new_depth_level_k4)
  
  k4_W<-update_d_depth(lon = lon_k4,
                       lat = lat_k4, 
                       dlon = k4_UV[1], 
                       dlat = k4_UV[2], 
                       depth = new_depth_k4,
                       dt, 
                       w_vertical = velocity4[3])
  
  
  # calculate the final new location of latitude and longitude
  dlon <- (k1_UV[1] + 2*k2_UV[1] + 2*k3_UV[1] + k4_UV[1])/6
  dlat <- (k1_UV[2] + 2*k2_UV[2] + 2*k3_UV[2] + k4_UV[2])/6
  d_depth <- (k1_W + 2*k2_W + 2*k3_W + k4_W)/6
  
  new_lon <- lon + dlon
  new_lat <- lat + dlat
  new_depth <- depth + d_depth*(-1)
  
  new_cell_num<-trans_lonlat_to_cellnumber(lon = new_lon,
                                           lat = new_lat)
  new_depth_level <- define_depth_level(depth = new_depth)
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==new_cell_num),'max_depth_level']
  
  # when eggs stunk on the seaside or the bottom floor side - retain the last timestep
  if (!new_cell_num %in% blacksea_points$Cell_num) {
    dlon <- 0
    dlat <- 0
    new_lon <- lon
    new_lat <- lat
    new_cell_num <- cell_number
    
  }else{
    # if particle would enter the right or left stone on the seafloor
    max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==new_cell_num),'max_depth_level']
    if(new_depth_level>max_depth_level_new){
      dlon <- 0
      dlat <- 0
      new_lon <- lon
      new_lat <- lat
      new_cell_num <- cell_number
      
      # let the particle keep upper the sea floor
      new_depth <- all_depth_levels_w[max_depth_level_new+1]-0.1 # still retain the particle near the sea floo with 0.1m
      d_depth <- depth - new_depth
      new_depth_level <- max_depth_level_new
    }
  }
  
  # when eggs float above the sea surface
  if (new_depth < 0) {
    new_depth<-0
    d_depth <- depth - new_depth
    new_depth_level <- 1
  }
  
  return(c(dlon, dlat, d_depth, new_lon,new_lat,new_cell_num,new_depth,new_depth_level))
  
}

# 5. running PTM--------------
PTM_function_3D<-function(input.stoke.drift,
                          state_variables,
                          dt,
                          current_uo_data,
                          next_uo_data,
                          current_vo_data,
                          next_vo_data,
                          current_wo_data,
                          next_wo_data,
                          current_Temp_data,
                          next_Temp_data,
                          current_Sal_data,
                          next_Sal_data,
                          current_rho_data,
                          next_rho_data,
                          current_sd_data,
                          next_sd_data,
                          hourly_step_num){
  
  running_status <- state_variables[1,'running_status']
  life_stage <- state_variables[1,'life_stage']
  
  lon <- state_variables[1,'Lon']
  lat <- state_variables[1,'Lat']
  depth <- state_variables[1,'depth']
  depth_level <- state_variables[1,'depth_level']
  PTM_X <- state_variables[1,'PTM_X']
  PTM_Y <- state_variables[1,'PTM_Y']
  PTM_Z <- state_variables[1,'PTM_Z']
  Cell_num <- state_variables[1,'Cell_num']
  dis_to_coast <- state_variables[1,'dis_to_coast']
  
  if(running_status == 'alive' & life_stage %in% life_stages_in_PTM){
    # if individual alive
    
    if(input.stoke.drift == 1){ # consider wave effect
      output<-RK4_PTM_temporal_interpolation_sd_3D(lon,lat, 
                                                   cell_number = Cell_num,
                                                   dt,
                                                   current_uo_data,
                                                   next_uo_data,
                                                   current_vo_data,
                                                   next_vo_data,
                                                   current_wo_data,
                                                   next_wo_data,
                                                   current_Temp_data,
                                                   next_Temp_data,
                                                   current_Sal_data,
                                                   next_Sal_data,
                                                   current_rho_data,
                                                   next_rho_data,
                                                   current_sd_data,
                                                   next_sd_data,
                                                   depth,
                                                   depth_level,
                                                   hourly_step_num)
    }else{ # do not consider stoke drift
      output<-RK4_PTM_temporal_interpolation_3D(lon,lat, 
                                                cell_number = Cell_num,
                                                dt,
                                                current_uo_data,
                                                next_uo_data,
                                                current_vo_data,
                                                next_vo_data,
                                                current_wo_data,
                                                next_wo_data,
                                                current_Temp_data,
                                                next_Temp_data,
                                                current_Sal_data,
                                                next_Sal_data,
                                                current_rho_data,
                                                next_rho_data,
                                                depth,
                                                depth_level,
                                                hourly_step_num)
      
    }
    
    new_dis_to_coast<-blacksea_points[match(output[6],blacksea_points$Cell_num),'dis_to_coast']
    
    PTM_output<-data.frame(PTM_X = output[1],
                           PTM_Y = output[2],
                           PTM_Z = output[3],
                           Lon = output[4],
                           Lat = output[5],
                           Cell_num = output[6],
                           depth = output[7],
                           depth_level = output[8],
                           dis_to_coast = new_dis_to_coast)
    
  }else{ # if they are dead or do not need PTM
    PTM_output<-data.frame(PTM_X = PTM_X,
                           PTM_Y = PTM_Y,
                           PTM_Z = PTM_Z,
                           Lon = lon,
                           Lat = lat,
                           Cell_num = Cell_num,
                           depth = depth,
                           depth_level = depth_level,
                           dis_to_coast = dis_to_coast)
  }
  
  return(PTM_output)
}
