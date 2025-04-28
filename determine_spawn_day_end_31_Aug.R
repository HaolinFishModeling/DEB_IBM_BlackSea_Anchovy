# A function to determine the earliest and last day of possible spawning days in each year 20240715
determine.spawn.day<-function(current.year){
  
  path1 <- '/gpfs/home/acad/ulg-mast/yuhaolin/rtest'
  path.temp<-paste0(path1,'/outputs_BAMHBI_real_values/Temp/')
  blacksea_reference_points<-readRDS(paste0(path1,'/blacksea_points.rds'))
  colnames(blacksea_reference_points)[9]<-'Cell_number'
  folder.sel <- list.files(path.temp)
  year.files <- grep(as.character(current.year), folder.sel, value = TRUE)
  
  all.temp<-data.frame()
  for (i in 1:length(year.files)) {
    aaa <- readRDS(paste0(path.temp, year.files[i]))
    bbb <- aaa %>%
      left_join(blacksea_reference_points[,c('EEZs_adj', 'Cell_number')], by='Cell_number') %>%
      group_by(EEZs_adj, Date) %>%
      summarise(min_Temp_0_10m=min(Temp_0_10m),
                mean_Temp_0_10m = mean(Temp_0_10m, na.rm = TRUE), 
                max_Temp_0_10m=max(Temp_0_10m), .groups = 'drop')
    all.temp <- bind_rows(all.temp, bbb)
  }
  
  all.temp$Date<-as.Date(all.temp$Date)
  # Step 1: find the earliest day when Temp_spawn > 16 and continue to increase in the following 3 days
  earliest_day<-
    all.temp %>%
      group_by(EEZs_adj) %>%
      arrange(Date) %>%
      filter(max_Temp_0_10m > 16) %>%
      mutate(rise_3_days = lead(mean_Temp_0_10m, 1) > mean_Temp_0_10m &
               lead(mean_Temp_0_10m, 2) > lead(mean_Temp_0_10m, 1) &
               lead(mean_Temp_0_10m, 3) > lead(mean_Temp_0_10m, 2)) %>%
      filter(rise_3_days) %>%
      summarise(first_rising_date = first(Date)) %>%
      ungroup() %>%
      summarise(earliest_rising_date = min(first_rising_date, na.rm = TRUE))%>%
      pull(earliest_rising_date)
  
  # Step 2: find the last day when Temp_spawn < 16 and previous three days' mean temperature continuously decreased 
  latest_day<-as.Date(paste0(current.year,'-08-31'))
  
  day.spawn<-data.frame(earliest_day=earliest_day,
                        latest_day=latest_day)
  return(day.spawn)
}