trans_lonlat_to_rowcol<-function(lon,lat){
  lon_interval<-0.025
  lat_interval<-0.025
  
  # according to raster of lon and lat, calculate lon and lat of cell[1,1], cell[1,578], cell[258,1],cell[258.578]
  # 1st col = 29.875 - (100 - 1) * 0.025 = 27.4-0.0125
  # 578th col = 37.375 + (578 - 400) * 0.025 = 41.825-0.0125
  # 1st row =  46.1 + (50 - 1) * 0.025 = 47.325+0.0125
  # 258th row = 42.1 - (258 - 210) * 0.025 = 40.9+0.0125
  
  min.lon<-27.4-0.0125
  max.lon<-41.825-0.0125
  min.lat<-40.9+0.0125
  max.lat<-47.325+0.0125
  
  lon_col <- floor((lon - min.lon) / lon_interval) + 1
  lat_row <- floor((max.lat - lat) / lat_interval) + 1
  
  if(lon_col<1 | lon_col > 578 | lat_row<1 | lat_row>258){
    lon_col<-NA
    lat_row<-NA
  }
  
  aaa<-cbind(lat_row,lon_col)
  return(aaa)
}

trans_lonlat_to_cellnumber<-function(lon,lat){
  lon_interval<-0.025
  lat_interval<-0.025
  
  # according to raster of lon and lat, calculate lon and lat of cell[1,1], cell[1,578], cell[258,1],cell[258.578]
  # 1st col = 29.875 - (100 - 1) * 0.025 = 27.4-0.0125
  # 578th col = 37.375 + (578 - 400) * 0.025 = 41.825-0.0125
  # 1st row =  46.1 + (50 - 1) * 0.025 = 47.325+0.0125
  # 258th row = 42.1 - (258 - 210) * 0.025 = 40.9+0.0125
  
  min.lon<-27.4-0.0125
  max.lon<-41.825-0.0125
  min.lat<-40.9+0.0125
  max.lat<-47.325+0.0125
  total.ncol<-578
  
  lon_col <- floor((lon - min.lon) / lon_interval) + 1
  lat_row <- floor((max.lat - lat) / lat_interval) + 1
  
  cell_number<-(lat_row-1)*total.ncol + lon_col
  
  if(lon<min.lon|lon>max.lon|lat<min.lat|lat>max.lat){
    cell_number<-0
  }
  
  return(cell_number)
}

trans_cellnum_to_rowcol<-function(cell.number){
  total.ncol<-578
  
  row_numbers <- floor((cell.number - 1) / total.ncol) + 1
  col_numbers <- (cell.number - 1) %% total.ncol + 1
  aaa<-cbind(row_numbers,col_numbers)
  return(aaa)
}

trans_rowcol_to_cellnum<-function(nrow_number,col_number){
  total.ncol<-578
  
  cell_number<-(nrow_number-1)*total.ncol + col_number
  return(cell_number)
}

trans_rowcol_to_lonlat <- function(lat_row, lon_col) {
  lon_interval <- 0.025
  lat_interval <- 0.025
  
  # 已知网格的边界经纬度
  min_lon <- 27.4 - 0.0125
  max_lon <- 41.825 - 0.0125
  min_lat <- 40.9 + 0.0125
  max_lat <- 47.325 + 0.0125
  
  # 网格中心点的经纬度计算
  lon <- min_lon + (lon_col - 0.5) * lon_interval
  lat <- max_lat - (lat_row - 0.5) * lat_interval
  
  # 检查行列号是否超出范围
  if (lon_col < 1 | lon_col > 578 | lat_row < 1 | lat_row > 258) {
    lon <- NA
    lat <- NA
  }
  
  result <- cbind(lon, lat)
  return(result)
}

trans_cellnum_to_lonlat <- function(cell.number){
  total.ncol<-578
  
  row_numbers <- floor((cell.number - 1) / total.ncol) + 1
  col_numbers <- (cell.number - 1) %% total.ncol + 1
  
  lon_interval <- 0.025
  lat_interval <- 0.025
  
  # 已知网格的边界经纬度
  min_lon <- 27.4 - 0.0125
  max_lon <- 41.825 - 0.0125
  min_lat <- 40.9 + 0.0125
  max_lat <- 47.325 + 0.0125
  
  # 网格中心点的经纬度计算
  lon <- min_lon + (col_numbers - 0.5) * lon_interval
  lat <- max_lat - (row_numbers - 0.5) * lat_interval
  
  # 检查行列号是否超出范围
  if (col_numbers < 1 | col_numbers > 578 | row_numbers < 1 | row_numbers > 258) {
    lon <- NA
    lat <- NA
  }
  
  result <- cbind(lon, lat)
  return(result)
}


rotate_clockwise <- function(x, y, theta) {
  a1<-rnorm(1,0,1)
  if(a1>0){
    x_new <- x * cos(theta) + y * sin(theta)
    y_new <- -x * sin(theta) + y * cos(theta)
  }else{
    x_new <- x * cos(theta) - y * sin(theta)
    y_new <- x * sin(theta) + y * cos(theta)
  }
  
  return(c(x_new, y_new))
}
