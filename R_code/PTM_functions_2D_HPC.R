# Particle tracking model (PTM) functions 2D --- 20250211

# 1. consider stokes drift-------------
interpolation_velocity_and_temporal_sd_2D<-function(lon, lat, 
                                                    current_uo_data,
                                                    next_uo_data,
                                                    current_vo_data,
                                                    next_vo_data,
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
  
  # determine location of cell for u0, v0, and sd0
  u0.lon<-blacksea_points_uo[which(blacksea_points_uo$Cell_num==target.cell.num),c('Longitude')]
  u0.lat<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Latitude')]
  v0.lon<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Longitude')]
  v0.lat<-blacksea_points_vo[which(blacksea_points_vo$Cell_num==target.cell.num),c('Latitude')]
  sd0.lon<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Longitude')]
  sd0.lat<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Latitude')]
  
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
  
  # extract velocity to u0, u1, u2, u3; v0, v1, v2, v3; sd.u0, sd.u1, sd.u2, sd.u3; and sd.v0, sd.v1, sd.v2, sd.v3
  if(target.cell.num%in%blacksea_points$Cell_num){
    original.u0<-current_uo_data[which(current_uo_data$Cell_number==target.cell.num),'uo_surface']
    original.v0<-current_vo_data[which(current_vo_data$Cell_number==target.cell.num),'vo_surface']
    original.sd.u0<-current_sd_data[which(current_sd_data$Cell_number==target.cell.num),'VSDX']
    original.sd.v0<-current_sd_data[which(current_sd_data$Cell_number==target.cell.num),'VSDY']
    
    next.u0<-next_uo_data[which(next_uo_data$Cell_number==target.cell.num),'uo_surface']
    next.v0<-next_vo_data[which(next_vo_data$Cell_number==target.cell.num),'vo_surface']
    next.sd.u0<-next_sd_data[which(next_sd_data$Cell_number==target.cell.num),'VSDX']
    next.sd.v0<-next_sd_data[which(next_sd_data$Cell_number==target.cell.num),'VSDY']
  }else{
    original.u0<-0
    original.v0<-0
    original.sd.u0<-0
    original.sd.v0<-0
    
    next.u0<-0
    next.v0<-0
    next.sd.u0<-0
    next.sd.v0<-0
  }
  
  if(u1.cell.num %in% blacksea_points$Cell_num){
    original.u1<-current_uo_data[which(current_uo_data$Cell_number==u1.cell.num),'uo_surface']
    next.u1<-next_uo_data[which(next_uo_data$Cell_number==u1.cell.num),'uo_surface']
  }else{
    original.u1<-0
    next.u1<-0
  }
  
  if(u2.cell.num %in% blacksea_points$Cell_num){
    original.u2<-current_uo_data[which(current_uo_data$Cell_number==u2.cell.num),'uo_surface']
    next.u2<-next_uo_data[which(next_uo_data$Cell_number==u2.cell.num),'uo_surface']
  }else{
    original.u2<-0
    next.u2<-0
  }
  
  if(u3.cell.num %in% blacksea_points$Cell_num){
    original.u3<-current_uo_data[which(current_uo_data$Cell_number==u3.cell.num),'uo_surface']
    next.u3<-next_uo_data[which(next_uo_data$Cell_number==u3.cell.num),'uo_surface']
  }else{
    original.u3<-0
    next.u3<-0
  }
  
  if(v1.cell.num %in% blacksea_points$Cell_num){
    original.v1<-current_vo_data[which(current_vo_data$Cell_number==v1.cell.num),'vo_surface']
    next.v1<-next_vo_data[which(next_vo_data$Cell_number==v1.cell.num),'vo_surface']
  }else{
    original.v1<-0
    next.v1<-0
  }
  
  if(v2.cell.num %in% blacksea_points$Cell_num){
    original.v2<-current_vo_data[which(current_vo_data$Cell_number==v2.cell.num),'vo_surface']
    next.v2<-next_vo_data[which(next_vo_data$Cell_number==v2.cell.num),'vo_surface']
  }else{
    original.v2<-0
    next.v2<-0
  }
  
  if(v3.cell.num %in% blacksea_points$Cell_num){
    original.v3<-current_vo_data[which(current_vo_data$Cell_number==v3.cell.num),'vo_surface']
    next.v3<-next_vo_data[which(next_vo_data$Cell_number==v3.cell.num),'vo_surface']
  }else{
    original.v3<-0
    next.v3<-0
  }
  
  if(sd1.cell.num %in% blacksea_points$Cell_num){
    original.sd.u1<-current_sd_data[which(current_sd_data$Cell_number==sd1.cell.num),'VSDX']
    next.sd.u1<-next_sd_data[which(next_sd_data$Cell_number==sd1.cell.num),'VSDX']
    original.sd.v1<-current_sd_data[which(current_sd_data$Cell_number==sd1.cell.num),'VSDY']
    next.sd.v1<-next_sd_data[which(next_sd_data$Cell_number==sd1.cell.num),'VSDY']
  }else{
    original.sd.u1<-0
    original.sd.v1<-0
    next.sd.u1<-0
    next.sd.v1<-0
  }
  
  if(sd2.cell.num %in% blacksea_points$Cell_num){
    original.sd.u2<-current_sd_data[which(current_sd_data$Cell_number==sd2.cell.num),'VSDX']
    next.sd.u2<-next_sd_data[which(next_sd_data$Cell_number==sd2.cell.num),'VSDX']
    original.sd.v2<-current_sd_data[which(current_sd_data$Cell_number==sd2.cell.num),'VSDY']
    next.sd.v2<-next_sd_data[which(next_sd_data$Cell_number==sd2.cell.num),'VSDY']
  }else{
    original.sd.u2<-0
    original.sd.v2<-0
    next.sd.u2<-0
    next.sd.v2<-0
  }
  
  if(sd3.cell.num %in% blacksea_points$Cell_num){
    original.sd.u3<-current_sd_data[which(current_sd_data$Cell_number==sd3.cell.num),'VSDX']
    next.sd.u3<-next_sd_data[which(next_sd_data$Cell_number==sd3.cell.num),'VSDX']
    original.sd.v3<-current_sd_data[which(current_sd_data$Cell_number==sd3.cell.num),'VSDY']
    next.sd.v3<-next_sd_data[which(next_sd_data$Cell_number==sd3.cell.num),'VSDY']
  }else{
    original.sd.u3<-0
    original.sd.v3<-0
    next.sd.u3<-0
    next.sd.v3<-0
  }
  
  # start temporal interpolation
  u0<-seq(from = original.u0, to = next.u0, length.out = 24)[hourly_step_num]
  v0<-seq(from = original.v0, to = next.v0, length.out = 24)[hourly_step_num]
  u1<-seq(from = original.u1, to = next.u1, length.out = 24)[hourly_step_num]
  v1<-seq(from = original.v1, to = next.v1, length.out = 24)[hourly_step_num]
  u2<-seq(from = original.u2, to = next.u2, length.out = 24)[hourly_step_num]
  v2<-seq(from = original.v2, to = next.v2, length.out = 24)[hourly_step_num]
  u3<-seq(from = original.u3, to = next.u3, length.out = 24)[hourly_step_num]
  v3<-seq(from = original.v3, to = next.v3, length.out = 24)[hourly_step_num]
  
  sd.u0<-seq(from = original.sd.u0, to = next.sd.u0, length.out = 24)[hourly_step_num]
  sd.v0<-seq(from = original.sd.v0, to = next.sd.v0, length.out = 24)[hourly_step_num]
  sd.u1<-seq(from = original.sd.u1, to = next.sd.u1, length.out = 24)[hourly_step_num]
  sd.v1<-seq(from = original.sd.v1, to = next.sd.v1, length.out = 24)[hourly_step_num]
  sd.u2<-seq(from = original.sd.u2, to = next.sd.u2, length.out = 24)[hourly_step_num]
  sd.v2<-seq(from = original.sd.v2, to = next.sd.v2, length.out = 24)[hourly_step_num]
  sd.u3<-seq(from = original.sd.u3, to = next.sd.u3, length.out = 24)[hourly_step_num]
  sd.v3<-seq(from = original.sd.v3, to = next.sd.v3, length.out = 24)[hourly_step_num]
  
  # strat to interpolate
  # u and v
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
  
  # sd.u and sd.v
  xsi<-(lon-sd0.lon)/(sd1.lon-sd0.lon)
  target.sd.u_1<-sd.u0 + xsi*(sd.u1-sd.u0)
  target.sd.v_1<-sd.v0 + xsi*(sd.v1-sd.v0)
  
  xsi<-(lon-sd2.lon)/(sd3.lon-sd2.lon)
  target.sd.u_2<-sd.u2 + xsi*(sd.u3-sd.u2)
  target.sd.v_2<-sd.v2 + xsi*(sd.v3-sd.v2)
  
  xsi<-(lat-sd1.lat)/(sd3.lat-sd1.lat)
  target.sd.u<-target.sd.u_1+xsi*(target.sd.u_2 - target.sd.u_1)
  target.sd.v<-target.sd.v_1+xsi*(target.sd.v_2 - target.sd.v_1)
  
  target.u<-target.u+target.sd.u
  target.v<-target.v+target.sd.v
  
  output<-c(target.u,target.v)
  return(output)
}

RK4_PTM_temporal_interpolation_sd_2D<-function(lon,lat, 
                                               dt,
                                               current_uo_data,
                                               next_uo_data,
                                               current_vo_data,
                                               next_vo_data,
                                               current_sd_data,
                                               next_sd_data,
                                               hourly_step_num){
  
  # calculate k1, k2, k3, k4
  # k1
  velocity1<-interpolation_velocity_and_temporal_sd_2D(lon,lat, 
                                                       current_uo_data,
                                                       next_uo_data,
                                                       current_vo_data,
                                                       next_vo_data,
                                                       current_sd_data,
                                                       next_sd_data,
                                                       hourly_step_num)
  k1<-update_dlon_dlat_2D(lon, lat, dt, velocity1)
  
  # k2
  lon_k2<-lon+k1[1]/2
  lat_k2<-lat+k1[2]/2
  cell.number_k2<-trans_lonlat_to_cellnumber(lon = lon_k2,
                                             lat = lat_k2)
  if(!cell.number_k2 %in% blacksea_points$Cell_num){
    lon_k2<-lon
    lat_k2<-lat
  }
  
  velocity2<-interpolation_velocity_and_temporal_sd_2D(lon = lon_k2, lat = lat_k2, 
                                                       current_uo_data,
                                                       next_uo_data,
                                                       current_vo_data,
                                                       next_vo_data,
                                                       current_sd_data,
                                                       next_sd_data,
                                                       hourly_step_num)
  k2<-update_dlon_dlat_2D(lon = lon_k2, lat = lat_k2, dt, velocity2)
  
  # k3
  lon_k3<-lon+k2[1]/2
  lat_k3<-lat+k2[2]/2
  cell.number_k3<-trans_lonlat_to_cellnumber(lon = lon_k3,
                                             lat = lat_k3)
  if(!cell.number_k3 %in% blacksea_points$Cell_num){
    lon_k3<-lon
    lat_k3<-lat
  }
  velocity3<-interpolation_velocity_and_temporal_sd_2D(lon = lon_k3, lat = lat_k3, 
                                                       current_uo_data,
                                                       next_uo_data,
                                                       current_vo_data,
                                                       next_vo_data,
                                                       current_sd_data,
                                                       next_sd_data,
                                                       hourly_step_num)
  k3<-update_dlon_dlat_2D(lon = lon_k3, lat = lat_k3, dt, velocity3)
  
  # k4
  lon_k4<-lon+k3[1]
  lat_k4<-lat+k3[2]
  cell.number_k4<-trans_lonlat_to_cellnumber(lon = lon_k4,
                                             lat = lat_k4)
  if(!cell.number_k4 %in% blacksea_points$Cell_num){
    lon_k4<-lon
    lat_k4<-lat
  }
  velocity4<-interpolation_velocity_and_temporal_sd_2D(lon = lon_k4, lat = lat_k4, 
                                                       current_uo_data,
                                                       next_uo_data,
                                                       current_vo_data,
                                                       next_vo_data,
                                                       current_sd_data,
                                                       next_sd_data,
                                                       hourly_step_num)
  k4<-update_dlon_dlat_2D(lon = lon_k4, lat = lat_k4, dt, velocity4)
  
  # calculate the final new location of latitude and longitude
  dlon <- (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
  dlat <- (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
  new_lon<-lon + dlon
  new_lat<-lat + dlat
  
  new_cell_num<-trans_lonlat_to_cellnumber(lon = new_lon,
                                           lat = new_lat)
  
  # when eggs stunk on the seaside
  if (!new_cell_num %in% blacksea_points$Cell_num) {
    dlon <- 0
    dlat <- 0
    new_lon <- lon
    new_lat <- lat
    new_cell_num <- trans_lonlat_to_cellnumber(lon = new_lon,
                                               lat = new_lat)
  }
  
  return(c(dlon,dlat,new_lon,new_lat,new_cell_num))
  
}


# 2. Do not consider Stokes drift-------------------
# ！！！！！！！！！！！
# ！！！！！！！！！
# ！！！！！
# This should be examined by different life stages
# #########################
# if we don't use PTM for late_larvae, we need to change the codes here!!!!!!!!!!!!!!!!!!!!!!!!!
interpolation_velocity_and_temporal_2D<-function(lon, lat, 
                                                 current_uo_data,
                                                 next_uo_data,
                                                 current_vo_data,
                                                 next_vo_data,
                                                 hourly_step_num,
                                                 life_stage){
  
  # bilineary interpolation and temporal interpolation
  
  # determine cell number
  target.cell.num<-trans_lonlat_to_cellnumber(lat = lat,
                                              lon = lon)
  
  # determine row number and column number
  target.row.con.num<-trans_lonlat_to_rowcol(lat = lat,
                                             lon = lon)
  
  # determine location of cell for u0 and v0
  u0.lon<-blacksea_points_uo[which(blacksea_points_uo$Cell_num==target.cell.num),c('Longitude')]
  u0.lat<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Latitude')]
  v0.lon<-blacksea_points[which(blacksea_points$Cell_num==target.cell.num),c('Longitude')]
  v0.lat<-blacksea_points_vo[which(blacksea_points_vo$Cell_num==target.cell.num),c('Latitude')]
  
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
  
  # extract velocity to u0, u1, u2, u3 and v0, v1, v2, v3
  if(life_stage %in% c('egg','yolksac_larvae')){
    if(target.cell.num%in%blacksea_points$Cell_num){
      original.u0<-current_uo_data[which(current_uo_data$Cell_number==target.cell.num),'uo_surface']
      original.v0<-current_vo_data[which(current_vo_data$Cell_number==target.cell.num),'vo_surface']
      
      next.u0<-next_uo_data[which(next_uo_data$Cell_number==target.cell.num),'uo_surface']
      next.v0<-next_vo_data[which(next_vo_data$Cell_number==target.cell.num),'vo_surface']
    }else{
      original.u0<-0
      original.v0<-0
      
      next.u0<-0
      next.v0<-0
    }
    
    if(u1.cell.num %in% blacksea_points$Cell_num){
      original.u1<-current_uo_data[which(current_uo_data$Cell_number==u1.cell.num),'uo_surface']
      next.u1<-next_uo_data[which(next_uo_data$Cell_number==u1.cell.num),'uo_surface']
    }else{
      original.u1<-0
      next.u1<-0
    }
    
    if(u2.cell.num %in% blacksea_points$Cell_num){
      original.u2<-current_uo_data[which(current_uo_data$Cell_number==u2.cell.num),'uo_surface']
      next.u2<-next_uo_data[which(next_uo_data$Cell_number==u2.cell.num),'uo_surface']
    }else{
      original.u2<-0
      next.u2<-0
    }
    
    if(u3.cell.num %in% blacksea_points$Cell_num){
      original.u3<-current_uo_data[which(current_uo_data$Cell_number==u3.cell.num),'uo_surface']
      next.u3<-next_uo_data[which(next_uo_data$Cell_number==u3.cell.num),'uo_surface']
    }else{
      original.u3<-0
      next.u3<-0
    }
    
    if(v1.cell.num %in% blacksea_points$Cell_num){
      original.v1<-current_vo_data[which(current_vo_data$Cell_number==v1.cell.num),'vo_surface']
      next.v1<-next_vo_data[which(next_vo_data$Cell_number==v1.cell.num),'vo_surface']
    }else{
      original.v1<-0
      next.v1<-0
    }
    
    if(v2.cell.num %in% blacksea_points$Cell_num){
      original.v2<-current_vo_data[which(current_vo_data$Cell_number==v2.cell.num),'vo_surface']
      next.v2<-next_vo_data[which(next_vo_data$Cell_number==v2.cell.num),'vo_surface']
    }else{
      original.v2<-0
      next.v2<-0
    }
    
    if(v3.cell.num %in% blacksea_points$Cell_num){
      original.v3<-current_vo_data[which(current_vo_data$Cell_number==v3.cell.num),'vo_surface']
      next.v3<-next_vo_data[which(next_vo_data$Cell_number==v3.cell.num),'vo_surface']
    }else{
      original.v3<-0
      next.v3<-0
    }
    
  }else if(life_stage %in% c('early_larvae','late_larvae')){
    
    
    
    
    
    # be careful if we want to include late larvae for PTM later
    # this if condition should be updated
    
    
    
    
    
    # early larvae use the water column from 0 to 30 meters
    
    if(target.cell.num%in%blacksea_points$Cell_num){
      original.u0<-current_uo_data[which(current_uo_data$Cell_number==target.cell.num),'uo_0_30m']
      original.v0<-current_vo_data[which(current_vo_data$Cell_number==target.cell.num),'vo_0_30m']
      
      next.u0<-next_uo_data[which(next_uo_data$Cell_number==target.cell.num),'uo_0_30m']
      next.v0<-next_vo_data[which(next_vo_data$Cell_number==target.cell.num),'vo_0_30m']
    }else{
      original.u0<-0
      original.v0<-0
      
      next.u0<-0
      next.v0<-0
    }
    
    if(u1.cell.num %in% blacksea_points$Cell_num){
      original.u1<-current_uo_data[which(current_uo_data$Cell_number==u1.cell.num),'uo_0_30m']
      next.u1<-next_uo_data[which(next_uo_data$Cell_number==u1.cell.num),'uo_0_30m']
    }else{
      original.u1<-0
      next.u1<-0
    }
    
    if(u2.cell.num %in% blacksea_points$Cell_num){
      original.u2<-current_uo_data[which(current_uo_data$Cell_number==u2.cell.num),'uo_0_30m']
      next.u2<-next_uo_data[which(next_uo_data$Cell_number==u2.cell.num),'uo_0_30m']
    }else{
      original.u2<-0
      next.u2<-0
    }
    
    if(u3.cell.num %in% blacksea_points$Cell_num){
      original.u3<-current_uo_data[which(current_uo_data$Cell_number==u3.cell.num),'uo_0_30m']
      next.u3<-next_uo_data[which(next_uo_data$Cell_number==u3.cell.num),'uo_0_30m']
    }else{
      original.u3<-0
      next.u3<-0
    }
    
    if(v1.cell.num %in% blacksea_points$Cell_num){
      original.v1<-current_vo_data[which(current_vo_data$Cell_number==v1.cell.num),'vo_0_30m']
      next.v1<-next_vo_data[which(next_vo_data$Cell_number==v1.cell.num),'vo_0_30m']
    }else{
      original.v1<-0
      next.v1<-0
    }
    
    if(v2.cell.num %in% blacksea_points$Cell_num){
      original.v2<-current_vo_data[which(current_vo_data$Cell_number==v2.cell.num),'vo_0_30m']
      next.v2<-next_vo_data[which(next_vo_data$Cell_number==v2.cell.num),'vo_0_30m']
    }else{
      original.v2<-0
      next.v2<-0
    }
    
    if(v3.cell.num %in% blacksea_points$Cell_num){
      original.v3<-current_vo_data[which(current_vo_data$Cell_number==v3.cell.num),'vo_0_30m']
      next.v3<-next_vo_data[which(next_vo_data$Cell_number==v3.cell.num),'vo_0_30m']
    }else{
      original.v3<-0
      next.v3<-0
    }
  }
  
  # start temporal interpolation
  u0<-seq(from = original.u0, to = next.u0, length.out = 24)[hourly_step_num]
  v0<-seq(from = original.v0, to = next.v0, length.out = 24)[hourly_step_num]
  u1<-seq(from = original.u1, to = next.u1, length.out = 24)[hourly_step_num]
  v1<-seq(from = original.v1, to = next.v1, length.out = 24)[hourly_step_num]
  u2<-seq(from = original.u2, to = next.u2, length.out = 24)[hourly_step_num]
  v2<-seq(from = original.v2, to = next.v2, length.out = 24)[hourly_step_num]
  u3<-seq(from = original.u3, to = next.u3, length.out = 24)[hourly_step_num]
  v3<-seq(from = original.v3, to = next.v3, length.out = 24)[hourly_step_num]
  
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
  
  output<-c(target.u,target.v)
  return(output)
}


RK4_PTM_temporal_interpolation_2D<-function(lon,lat, 
                                            dt,
                                            current_uo_data,
                                            next_uo_data,
                                            current_vo_data,
                                            next_vo_data,
                                            hourly_step_num,
                                            life_stage){
  
  # calculate k1, k2, k3, k4
  # k1
  velocity1<-interpolation_velocity_and_temporal_2D(lon,lat, 
                                                    current_uo_data,
                                                    next_uo_data,
                                                    current_vo_data,
                                                    next_vo_data,
                                                    hourly_step_num,
                                                    life_stage)
  k1<-update_dlon_dlat_2D(lon, lat, dt, velocity1)
  
  # k2
  lon_k2<-lon+k1[1]/2
  lat_k2<-lat+k1[2]/2
  cell.number_k2<-trans_lonlat_to_cellnumber(lon = lon_k2,
                                             lat = lat_k2)
  if(!cell.number_k2 %in% blacksea_points$Cell_num){
    lon_k2<-lon
    lat_k2<-lat
  }
  
  velocity2<-interpolation_velocity_and_temporal_2D(lon = lon_k2, lat = lat_k2, 
                                                    current_uo_data,
                                                    next_uo_data,
                                                    current_vo_data,
                                                    next_vo_data,
                                                    hourly_step_num,
                                                    life_stage)
  k2<-update_dlon_dlat_2D(lon = lon_k2, lat = lat_k2, dt, velocity2)
  
  # k3
  lon_k3<-lon+k2[1]/2
  lat_k3<-lat+k2[2]/2
  cell.number_k3<-trans_lonlat_to_cellnumber(lon = lon_k3,
                                             lat = lat_k3)
  if(!cell.number_k3 %in% blacksea_points$Cell_num){
    lon_k3<-lon
    lat_k3<-lat
  }
  velocity3<-interpolation_velocity_and_temporal_2D(lon = lon_k3, lat = lat_k3, 
                                                    current_uo_data,
                                                    next_uo_data,
                                                    current_vo_data,
                                                    next_vo_data,
                                                    hourly_step_num,
                                                    life_stage)
  k3<-update_dlon_dlat_2D(lon = lon_k3, lat = lat_k3, dt, velocity3)
  
  # k4
  lon_k4<-lon+k3[1]
  lat_k4<-lat+k3[2]
  cell.number_k4<-trans_lonlat_to_cellnumber(lon = lon_k4,
                                             lat = lat_k4)
  if(!cell.number_k4 %in% blacksea_points$Cell_num){
    lon_k4<-lon
    lat_k4<-lat
  }
  velocity4<-interpolation_velocity_and_temporal_2D(lon = lon_k4, lat = lat_k4, 
                                                    current_uo_data,
                                                    next_uo_data,
                                                    current_vo_data,
                                                    next_vo_data,
                                                    hourly_step_num,
                                                    life_stage)
  k4<-update_dlon_dlat_2D(lon = lon_k4, lat = lat_k4, dt, velocity4)
  # calculate the final new location of latitude and longitude
  dlon <- (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
  dlat <- (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
  new_lon<-lon + dlon
  new_lat<-lat + dlat
  
  new_cell_num<-trans_lonlat_to_cellnumber(lon = new_lon,
                                           lat = new_lat)
  
  # when eggs stunk on the seaside
  if (!new_cell_num %in% blacksea_points$Cell_num) {
    dlon <- 0
    dlat <- 0
    new_lon <- lon
    new_lat <- lat
    new_cell_num <- trans_lonlat_to_cellnumber(lon = new_lon,
                                               lat = new_lat)
  }
  
  return(c(dlon, dlat,new_lon,new_lat,new_cell_num))
  
}

# 3. runnig PTM
update_dlon_dlat_2D<-function(lon, lat, dt, velocity){
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
  }
  return(c(dlon, dlat))
}

PTM_function_2D<-function(input.stoke.drift,
                          state_variables,
                          dt,
                          current_uo_data,
                          next_uo_data,
                          current_vo_data,
                          next_vo_data,
                          current_sd_data,
                          next_sd_data,
                          hourly_step_num){
  
  running_status <- state_variables[1,'running_status']
  life_stage <- state_variables[1,'life_stage']
  
  lon <- state_variables[1,'Lon']
  lat <- state_variables[1,'Lat']
  PTM_X <- state_variables[1,'PTM_X']
  PTM_Y <- state_variables[1,'PTM_Y']
  PTM_Z <- state_variables[1,'PTM_Z']
  Cell_num <- state_variables[1,'Cell_num']
  dis_to_coast <- state_variables[1,'dis_to_coast']
  
  if(running_status == 'alive' & life_stage %in% life_stages_in_PTM){
    # if individual alive
    
    if(input.stoke.drift == 1 & life_stage %in% life_stages_in_Stokes_drift){ # consider wave effect
      output<-RK4_PTM_temporal_interpolation_sd_2D(lon,lat, 
                                                   dt,
                                                   current_uo_data,
                                                   next_uo_data,
                                                   current_vo_data,
                                                   next_vo_data,
                                                   current_sd_data,
                                                   next_sd_data,
                                                   hourly_step_num)
    }else{ # do not consider stoke drift, OR if early_larvae with PTM also DOES NOT consider wave effect
      output<-RK4_PTM_temporal_interpolation_2D(lon,lat, 
                                                dt,
                                                current_uo_data,
                                                next_uo_data,
                                                current_vo_data,
                                                next_vo_data,
                                                hourly_step_num,
                                                life_stage)
      
    }
    
    new_dis_to_coast<-blacksea_points[match(output[5],blacksea_points$Cell_num),'dis_to_coast']
    
    PTM_output<-data.frame(PTM_X = output[1],
                           PTM_Y = output[2],
                           PTM_Z = PTM_Z,
                           Lon = output[3],
                           Lat = output[4],
                           Cell_num = output[5],
                           dis_to_coast = new_dis_to_coast)
    
  }else{ # if they are dead or not belong to the target life stage
    PTM_output<-data.frame(PTM_X = PTM_X,
                           PTM_Y = PTM_Y,
                           PTM_Z = PTM_Z,
                           Lon = lon,
                           Lat = lat,
                           Cell_num = Cell_num,
                           dis_to_coast = dis_to_coast)
  }
  
  return(PTM_output)
}