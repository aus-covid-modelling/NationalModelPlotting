# # # # # # # # # # # # # # # # # # # # 
# Functions to produce occupancy figures
# 
# David J Price
# 19 July 2021
# 

# Convert array from MATLAB into data.frame
convert.mat <- function(array.dat, name){
  df <- plyr::adply(array.dat, c(1,2,3)) %>% 
    select(-X2, -X3)
  
  colnames(df) <- c("var","sim","day",age.groups)
  
  long.df <- df %>% 
    pivot_longer(cols = `0_5`:`80+`, names_to = "age.groups", values_to = "value") %>% 
    # separate(col = age.groups, into = c("lower","upper"), sep = "_") %>% 
    mutate(group = name)
  
  return(long.df)
}



# Convert array from MATLAB into data.frame, aggregating across age groups
convert.mat.no.age <- function(array.dat, name){
  df <- plyr::adply(array.dat, c(1,2,3)) %>% 
    select(-X2, -X3)
  
  colnames(df) <- c("var","sim","day", age.groups)
  
  df <- df %>% mutate(
    value = `0_5`+`5_10`+`10_15`+
      `15_20`+`20_25`+`25_30`+`30_35`+`35_40`+
      `40_45`+`45_50`+`50_55`+`55_60`+
      `60_65`+`65_70`+`70_75`+`75_80`+`80+`) %>% 
    dplyr::select(var, sim, day, value) %>% 
    mutate(group = name)
  
  return(df)
}


# Automating y-axis limits
my_lims <- function(x){
  c(0, max(x,20))
}
