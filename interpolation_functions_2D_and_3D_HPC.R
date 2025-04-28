# create interpolation functions 20250120

Temp_2D_spatial_interpolation<-function(var_data, life_stage, cell_number, lon, lat){
  # this is not used for PTM, this is for Kinesis and DEB for temperature and food
  # PTM has its specific interpolation function
  
  # bilineary interpolation
  
  # determine cell number
  target.cell.num<-cell_number
  
  # determine row number and column number
  target.row.con.num<-trans_lonlat_to_rowcol(lat = lat,
                                             lon = lon)
  target.row.con.num<-as.data.frame(target.row.con.num)
  
  # determine location of cell
  idx <- which(blacksea_points$Cell_num == target.cell.num)
  var0.lon<-blacksea_points$Longitude[idx]
  var0.lat<-blacksea_points$Latitude[idx]
  
  # 预计算 row 和 col，避免重复计算
  row_num <- target.row.con.num$lat_row
  col_num <- target.row.con.num$lon_col
  
  # locating 4 cells for var0, var1, var2, var3
  if(lat>= var0.lat & lon<= var0.lon){ # on the upper left panel
    
    var_shift <- list(c(0, -1), c(-1, 0), c(-1, -1))  # 左边、上边、左上角
    
  }else if(lat>= var0.lat & lon > var0.lon){ # on the upper right panel
    
    var_shift <- list(c(0, 1), c(-1, 0), c(-1, 1))  # 右边、上边、右上角
    
  }else if(lat < var0.lat & lon <= var0.lon){ # on the bottom left panel
    
    var_shift <- list(c(0, -1), c(1, 0), c(1, -1))  # 左边、下边、左下角
    
  }else{ # on the bottom right panel
    
    var_shift <- list(c(0, 1), c(1, 0), c(1, 1))  # 右边、下边、右下角
    
  }
  
  get_var_cell_info <- function(shift) {
    row_new <- row_num + shift[1]
    col_new <- col_num + shift[2]
    cell_new <- trans_rowcol_to_cellnum(row_new, col_new)
    lonlat <- trans_rowcol_to_lonlat(row_new, col_new)
    return(list(cell = cell_new, lon = lonlat[1], lat = lonlat[2]))
  }
  
  # 获取 var1, var2, var3
  var_info <- lapply(var_shift, get_var_cell_info)
  
  var1 <- var_info[[1]]
  var2 <- var_info[[2]]
  var3 <- var_info[[3]]
  
  # extract var0, var1, var2, var3
  # target.cell.num cannot be the cell in the land!!!
  nrow_num<-which(var_data$Cell_number==target.cell.num)
  if(life_stage %in% c('early_larvae', 'late_larvae')){
    var0.val<-var_data[nrow_num,'Temp_0_30m']
  }else if(life_stage %in% c('juvenile', 'adult')){
    var0.val<-var_data[nrow_num,'Temp_0_100m']
  }else if(life_stage %in% c('egg', 'yolksac_larvae')){
    var0.val<-var_data[nrow_num,'Temp_0_10m']
  }else{
    var0.val<-NA
  }
  var0_Temp_spawn<-var_data[nrow_num,'Temp_0_10m']
  
  if(var1$cell %in% blacksea_points$Cell_num){
    nrow_num<-which(var_data$Cell_number==var1$cell)
    
    if(life_stage %in% c('early_larvae', 'late_larvae')){
      var1.val<-var_data[nrow_num,'Temp_0_30m']
    }else if(life_stage %in% c('juvenile', 'adult')){
      var1.val<-var_data[nrow_num,'Temp_0_100m']
    }else if(life_stage %in% c('egg', 'yolksac_larvae')){
      var1.val<-var_data[nrow_num,'Temp_0_10m']
    }else{
      var1.val<-NA
    }
    var1_Temp_spawn<-var_data[nrow_num,'Temp_0_10m']
    
  }else{
    var1.val<-var0.val
    var1_Temp_spawn<-var0_Temp_spawn
  }
  
  if(var2$cell %in% blacksea_points$Cell_num){
    nrow_num<-which(var_data$Cell_number==var2$cell)
    
    if(life_stage %in% c('early_larvae', 'late_larvae')){
      var2.val<-var_data[nrow_num,'Temp_0_30m']
    }else if(life_stage %in% c('juvenile', 'adult')){
      var2.val<-var_data[nrow_num,'Temp_0_100m']
    }else if(life_stage %in% c('egg', 'yolksac_larvae')){
      var2.val<-var_data[nrow_num,'Temp_0_10m']
    }else{
      var2.val<-NA
    }
    var2_Temp_spawn<-var_data[nrow_num,'Temp_0_10m']
    
  }else{
    var2.val<-var0.val
    var2_Temp_spawn<-var0_Temp_spawn
  }
  
  if(var3$cell %in% blacksea_points$Cell_num){
    nrow_num<-which(var_data$Cell_number==var3$cell)
    
    if(life_stage %in% c('early_larvae', 'late_larvae')){
      var3.val<-var_data[nrow_num,'Temp_0_30m']
    }else if(life_stage %in% c('juvenile', 'adult')){
      var3.val<-var_data[nrow_num,'Temp_0_100m']
    }else if(life_stage %in% c('egg', 'yolksac_larvae')){
      var3.val<-var_data[nrow_num,'Temp_0_10m']
    }else{
      var3.val<-NA
    }
    var3_Temp_spawn<-var_data[nrow_num,'Temp_0_10m']
    
  }else{
    var3.val<-var0.val
    var3_Temp_spawn<-var0_Temp_spawn
  }
  
  # strat to interpolate
  xsi<-(lon-var0.lon)/(var1$lon-var0.lon)
  target.var_1<-var0.val + xsi*(var1.val-var0.val)
  target.var.spawn_1<-var0_Temp_spawn + xsi*(var1_Temp_spawn-var0_Temp_spawn)
  
  xsi<-(lon-var2$lon)/(var3$lon-var2$lon)
  target.var_2<-var2.val + xsi*(var3.val-var2.val)
  target.var.spawn_2<-var2_Temp_spawn + xsi*(var3_Temp_spawn-var2_Temp_spawn)
  
  xsi<-(lat-var1$lat)/(var3$lat-var1$lat)
  target.var<-target.var_1+xsi*(target.var_2 - target.var_1)
  target.var.spawn<-target.var.spawn_1+xsi*(target.var.spawn_2 - target.var.spawn_1)
  
  return(c(target.var,
           target.var.spawn))
}

Food_2D_spatial_interpolation<-function(var_data, life_stage, cell_number, lon, lat){
  # this is not used for PTM, this is for Kinesis and DEB for temperature and food
  # PTM has its specific interpolation function
  
  # bilineary interpolation
  # determine cell number
  target.cell.num<-cell_number
  
  # determine row number and column number
  target.row.con.num<-trans_lonlat_to_rowcol(lat = lat,
                                             lon = lon)
  target.row.con.num<-as.data.frame(target.row.con.num)
  
  # determine location of cell
  idx <- which(blacksea_points$Cell_num == target.cell.num)
  MIC0.lon<-blacksea_points$Longitude[idx]
  MIC0.lat<-blacksea_points$Latitude[idx]
  
  # 预计算 row 和 col，避免重复计算
  row_num <- target.row.con.num$lat_row
  col_num <- target.row.con.num$lon_col
  
  # locating 4 cells for var0, var1, var2, var3
  if(lat>= MIC0.lat & lon<= MIC0.lon){ # on the upper left panel
    
    var_shift <- list(c(0, -1), c(-1, 0), c(-1, -1))  # 左边、上边、左上角
    
  }else if(lat>= MIC0.lat & lon > MIC0.lon){ # on the upper right panel
    
    var_shift <- list(c(0, 1), c(-1, 0), c(-1, 1))  # 右边、上边、右上角
    
  }else if(lat < MIC0.lat & lon <= MIC0.lon){ # on the bottom left panel
    
    var_shift <- list(c(0, -1), c(1, 0), c(1, -1))  # 左边、下边、左下角
    
  }else{ # on the bottom right panel
    
    var_shift <- list(c(0, 1), c(1, 0), c(1, 1))  # 右边、下边、右下角
    
  }
  
  get_var_cell_info <- function(shift) {
    row_new <- row_num + shift[1]
    col_new <- col_num + shift[2]
    cell_new <- trans_rowcol_to_cellnum(row_new, col_new)
    lonlat <- trans_rowcol_to_lonlat(row_new, col_new)
    return(list(cell = cell_new, lon = lonlat[1], lat = lonlat[2]))
  }
  
  # 获取 var1, var2, var3
  var_info <- lapply(var_shift, get_var_cell_info)
  
  MIC1 <- var_info[[1]]
  MIC2 <- var_info[[2]]
  MIC3 <- var_info[[3]]
  
  # extract MIC0, MIC1, MIC2, MIC3 and MES0, MES1, MES2, MES3
  # target.cell.num cannot be the cell in the land!!!
  nrow_num<-which(var_data$Cell_number==target.cell.num)
  if(life_stage %in% c('early_larvae', 'late_larvae')){
    MIC0.val<-var_data[nrow_num,'Food_MIC_0_30m']
    MES0.val<-var_data[nrow_num,'Food_MES_0_30m']
  }else if(life_stage %in% c('juvenile', 'adult')){
    MIC0.val<-var_data[nrow_num,'Food_MIC_0_100m']
    MES0.val<-var_data[nrow_num,'Food_MES_0_100m']
  }
  
  if(MIC1$cell %in% blacksea_points$Cell_num){
    nrow_num<-which(var_data$Cell_number==MIC1$cell)
    
    if(life_stage %in% c('early_larvae', 'late_larvae')){
      MIC1.val<-var_data[nrow_num,'Food_MIC_0_30m']
      MES1.val<-var_data[nrow_num,'Food_MES_0_30m']
    }else if(life_stage %in% c('juvenile', 'adult')){
      MIC1.val<-var_data[nrow_num,'Food_MIC_0_100m']
      MES1.val<-var_data[nrow_num,'Food_MES_0_100m']
    }
    
  }else{
    MIC1.val<-MIC0.val
    MES1.val<-MES0.val
  }
  
  if(MIC2$cell %in% blacksea_points$Cell_num){
    nrow_num<-which(var_data$Cell_number==MIC2$cell)
    
    if(life_stage %in% c('early_larvae', 'late_larvae')){
      MIC2.val<-var_data[nrow_num,'Food_MIC_0_30m']
      MES2.val<-var_data[nrow_num,'Food_MES_0_30m']
    }else if(life_stage %in% c('juvenile', 'adult')){
      MIC2.val<-var_data[nrow_num,'Food_MIC_0_100m']
      MES2.val<-var_data[nrow_num,'Food_MES_0_100m']
    }
    
  }else{
    MIC2.val<-MIC0.val
    MES2.val<-MES0.val
  }
  
  if(MIC3$cell %in% blacksea_points$Cell_num){
    nrow_num<-which(var_data$Cell_number==MIC3$cell)
    
    if(life_stage %in% c('early_larvae', 'late_larvae')){
      MIC3.val<-var_data[nrow_num,'Food_MIC_0_30m']
      MES3.val<-var_data[nrow_num,'Food_MES_0_30m']
    }else if(life_stage %in% c('juvenile', 'adult')){
      MIC3.val<-var_data[nrow_num,'Food_MIC_0_100m']
      MES3.val<-var_data[nrow_num,'Food_MES_0_100m']
    }
    
  }else{
    MIC3.val<-MIC0.val
    MES3.val<-MES0.val
  }
  
  # strat to interpolate
  xsi<-(lon-MIC0.lon)/(MIC1$lon-MIC0.lon)
  target.MIC_1<-MIC0.val + xsi*(MIC1.val-MIC0.val)
  target.MES_1<-MES0.val + xsi*(MES1.val-MES0.val)
  
  xsi<-(lon-MIC2$lon)/(MIC3$lon-MIC2$lon)
  target.MIC_2<-MIC2.val + xsi*(MIC3.val-MIC2.val)
  target.MES_2<-MES2.val + xsi*(MES3.val-MES2.val)
  
  xsi<-(lat-MIC1$lat)/(MIC3$lat-MIC1$lat)
  target.MIC<-target.MIC_1+xsi*(target.MIC_2 - target.MIC_1)
  target.MES<-target.MES_1+xsi*(target.MES_2 - target.MES_1)
  
  output<-data.frame(MIC=target.MIC,
                     MES=target.MES)
  
  return(output)
}

var_3D_one_level_horizontal_spatial_interpolation <- function(lon, lat, var_data, depth_level) {
  
  # 预计算 cell_number 及 row_col_number
  target.cell.num <- trans_lonlat_to_cellnumber(lat = lat, lon = lon)
  target.row.con.num <- trans_lonlat_to_rowcol(lat = lat, lon = lon)
  target.row.con.num <- as.data.frame(target.row.con.num)
  
  # 提前存储查询到的 cell 数据，避免重复 `which()`
  idx <- which(blacksea_points$Cell_num == target.cell.num)
  if (length(idx) == 0) return(NA) # 直接返回 NA 避免后续错误
  
  var0.lon <- blacksea_points$Longitude[idx]
  var0.lat <- blacksea_points$Latitude[idx]
  
  # 预计算 row 和 col，避免重复计算
  row_num <- target.row.con.num$lat_row
  col_num <- target.row.con.num$lon_col
  
  # locating 4 cells for var0, var1, var2, var3
  if(lat>= var0.lat & lon<= var0.lon){ # on the upper left panel
    
    var_shift <- list(c(0, -1), c(-1, 0), c(-1, -1))  # 左边、上边、左上角
    
  }else if(lat>= var0.lat & lon > var0.lon){ # on the upper right panel
    
    var_shift <- list(c(0, 1), c(-1, 0), c(-1, 1))  # 右边、上边、右上角
    
  }else if(lat < var0.lat & lon <= var0.lon){ # on the bottom left panel
    
    var_shift <- list(c(0, -1), c(1, 0), c(1, -1))  # 左边、下边、左下角
    
  }else{ # on the bottom right panel
    
    var_shift <- list(c(0, 1), c(1, 0), c(1, 1))  # 右边、下边、右下角
    
  }
  
  get_var_cell_info <- function(shift) {
    row_new <- row_num + shift[1]
    col_new <- col_num + shift[2]
    cell_new <- trans_rowcol_to_cellnum(row_new, col_new)
    lonlat <- trans_rowcol_to_lonlat(row_new, col_new)
    return(list(cell = cell_new, lon = lonlat[1], lat = lonlat[2]))
  }
  
  # 获取 var1, var2, var3
  var_info <- lapply(var_shift, get_var_cell_info)
  
  var1 <- var_info[[1]]
  var2 <- var_info[[2]]
  var3 <- var_info[[3]]
  
  # 取出变量数据，避免 `var_data[cell_number,depth_level]` 重复索引
  var_data_depth <- var_data[,depth_level]
  
  # 定义一个安全索引函数，防止 NA 访问
  safe_get_var <- function(cell_num, default) {
    if (cell_num %in% blacksea_points$Cell_num) {
      return(ifelse(is.na(var_data_depth[cell_num]), 0, var_data_depth[cell_num]))
    } else {
      return(default)
    }
  }
  
  var0.val <- safe_get_var(target.cell.num, 0)
  var1.val <- safe_get_var(var1$cell, var0.val)
  var2.val <- safe_get_var(var2$cell, var0.val)
  var3.val <- safe_get_var(var3$cell, var0.val)
  
  # 计算插值
  xsi1 <- (lon - var0.lon) / (var1$lon - var0.lon)
  target.var_1 <- var0.val + xsi1 * (var1.val - var0.val)
  
  xsi2 <- (lon - var2$lon) / (var3$lon - var2$lon)
  target.var_2 <- var2.val + xsi2 * (var3.val - var2.val)
  
  xsi3 <- (lat - var1$lat) / (var3$lat - var1$lat)
  target.var <- target.var_1 + xsi3 * (target.var_2 - target.var_1)
  
  return(target.var)
}

var_3D_vertical_spatial_interpolation<-function(f0_level_var,
                                                f1_level_var,
                                                cell_number,
                                                depth,
                                                depth_level){ 
  
  max_depth_level_new<-max_depth_level[which(max_depth_level$Cell_num==cell_number),'max_depth_level']
  
  # for variable in T - center of the cell
  if(depth<=all_depth_levels_t[1] | depth>=all_depth_levels_t[max_depth_level_new] | depth == all_depth_levels_t[depth_level]){ # upper the top cell or below the bottom cell
    target.var.vertical<-f0_level_var
    
  }else{
    
    if(depth < all_depth_levels_t[depth_level]){
      next_depth_level<-depth_level-1
    }else{
      next_depth_level<-depth_level+1
    }
    f0_level_depth_t<-all_depth_levels_t[depth_level] 
    f1_level_depth_t<-all_depth_levels_t[next_depth_level]
    
    zeta_t<-(depth-f0_level_depth_t)/(f1_level_depth_t - f0_level_depth_t)
    target.var.vertical<-f0_level_var+zeta_t*(f1_level_var - f0_level_var)
    
  }
  
  return(target.var.vertical)
}

# temporal interpolation for hourly time step
temporal_interpolation_2D_3D<-function(current_value,
                                       next_value,
                                       hourly_step_num){
  output.value<-seq(from = current_value, to = next_value, length.out = 24)[hourly_step_num]
}
