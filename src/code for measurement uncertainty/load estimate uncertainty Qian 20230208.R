## Estimating uncertainty in constituent annual load estimation based on interpolation and beale ratio methods
## Baihua Fu  ## Edited by Qian Wang 
## 9/12/2021  ## 08/02/2023
## This edited version changes the way of calculating load uncertainty assuming that there is a linear relationship between the uncertainty of L_k-1 and L_k

# Set work space----
### set word directory
setwd("C:\\Users\\u1066632\\OneDrive - Australian National University\\WORKANU\\Projects\\GraphEvalUnc\\src\\code for measurement uncertainty\\")
### load libraries
library(dplyr)
library(ggplot2)
library(zoo)
library(gridExtra)
library(lubridate)

# Generate functions ---------------------------------------------------------------

### Calculate hydrology year function
as.year <- function(x) as.integer(as.yearmon(x) - 6/12)

### Linear interpolation function 
# Test the functions below
conclinear <- function(flowdf, concdf, minconc){
  flowzoo <- read.zoo(flowdf, header = TRUE)
  conczoo <- read.zoo(concdf, header = TRUE, aggregate = mean) #take average for same time step
  
  # merge flow and concentration data
  flowconczoo <- merge(flowzoo, conczoo)
  names(flowconczoo) <- c("flow", "conc")
  
  # assign tiedowns: begin with half of min value of concentration, and end of the min value.
  flowconczoo[1,"conc"] <- minconc/2
  flowconczoo$conc[nrow(flowconczoo)] <- minconc
  
  # interpolate concentration
  flowconczoo$conc <- na.approx(flowconczoo$conc, rule = 1)
  return(flowconczoo)
  
}


loadlinear <- function(flowconczoo, seconds){
  
  # hourly load, mg/L * m3/s = g/s
  flowconczoo$load <- flowconczoo$conc*flowconczoo$flow
  
 # aggregate to annual mean instant flow, concentration and load
  yearly <- aggregate(flowconczoo, as.year(time(flowconczoo)), mean)
  yearly <- as.data.frame(yearly)
  yearly <- yearly[-1,]  # removed the first first row (10 data points, i.e. 10 hours) due to unresolved issue on timezone
  
  # calculate annual total load via: multiply mean instant load by number of seconds in that year, then convert unit to t/yr
  yearly$load <- round(yearly$load*seconds * 0.000001, 3) #g/s * s/yr = g/yr. 10^-6g/yr = t/yr
  return(yearly)
}


loadunc <- function(flowconczoo, concunc, flowunc, concdf, minconc, covar=1){
  conczoo <- read.zoo(concdf, header = TRUE, aggregate = mean)
  # Add tiedowns to the observed concentration data so as to enable interpolation
  concdfnew <- data.frame(Datetime = index(conczoo), conc = as.data.frame(conczoo)$conc)
  concdfnew <- rbind(concdfnew, data.frame(Datetime=flowunc$Datetime[1], conc=minconc/2))
  concdfnew <- rbind(concdfnew, data.frame(Datetime=flowunc$Datetime[length(flowunc$Datetime)], conc=minconc))
  concdfnew <- concdfnew[order(concdfnew$Datetime),]
  df_temp <- as.data.frame(flowconczoo)
  flowcondf <- data.frame(Datetime=flowunc$Datetime, flow=df_temp$flow, conc=df_temp$conc)
  flowcondf$loadraw <- flowcondf$flow * flowcondf$conc
  time_conc_sample = concdfnew$Datetime
  sample_conc_length <- length(time_conc_sample)
  # assign initial values used for calculation
  flowcondf$slope = 1
  concunc$concuncsq <- concunc$concuncbound**2
  names(concdfnew) <- c('Datetime', 'conc')
  flowcondf$flowuncbound <- flowunc$flowuncbound # Assign values of flow uncertainty bounds
  flowcondf$concuncbound <- 0
  flowcondf$conc_delta_i <- 0
  flowcondf$conc_delta_i2 <- 0
  time_end <- time_conc_sample[sample_conc_length]
  flow_time_end_bool <- flowcondf$Datetime == time_end
  flowcondf[flowcondf$Datetime == time_conc_sample[1],"concuncbound"] <- max(concunc$concuncbound)
  flowcondf[flow_time_end_bool, "concuncbound"] <- max(concunc$concuncbound)
  flowcondf[flow_time_end_bool, "conc_delta_i"] <-
    concdfnew[concdfnew$Datetime == time_end, "conc"] * flowcondf[flow_time_end_bool, "concuncbound"]
  
  flowcondf[flow_time_end_bool, "conc_delta_i2"] <-
    concdfnew[concdfnew$Datetime == time_end, "conc"] * flowcondf[flow_time_end_bool, "concuncbound"]
  
  
  for (i in 1:(sample_conc_length - 1)){
    if (i> 1){
      
      flowcondf[flowcondf$Datetime == time_conc_sample[i], "concuncbound"] <- mean(concunc[concunc$Datetime==time_conc_sample[i], "concuncbound"])
    }
  }
  
  for (i in 1:(sample_conc_length - 1)){
    # If a data is obtained by interpolation, calculate the concentration uncertainty using the eq. in chapter 4.
    diff_time_sample <- as.integer(difftime(time_conc_sample[i + 1], time_conc_sample[i], units = 'hours'))
    
    if (diff_time_sample > 1){
      # browser()
      t1_bool <- flowcondf$Datetime >= time_conc_sample[i]
      t2_bool <- flowcondf$Datetime < time_conc_sample[i+1]
      # # Vectorize the calculation to speed up
      flowcondf[t1_bool & t2_bool, "slope"] <- c(seq(diff_time_sample, 1) / diff_time_sample)
      flowcondf[t1_bool & t2_bool, "conc_delta_i"] <-
        concdfnew[concdfnew$Datetime == time_conc_sample[i], "conc"] * flowcondf[flowcondf$Datetime == time_conc_sample[i], "concuncbound"]
      
      flowcondf[t1_bool & t2_bool, "conc_delta_i2"] <-
        concdfnew[concdfnew$Datetime == time_conc_sample[i+1], "conc"] * flowcondf[flowcondf$Datetime == time_conc_sample[i+1], "concuncbound"]
    }
    # END if
  }
  if (covar==1){
    flowcondf[, "concuncbound"] <- ((flowcondf[, "slope"]*flowcondf[, "conc_delta_i"])** 2 + ((1 - flowcondf[, "slope"])*flowcondf[, "conc_delta_i2"])** 2)
    flowcondf[, "concuncbound"] <- round(sqrt(flowcondf[, "concuncbound"]) / flowcondf[, "conc"], 3)
    # Note that flowcondf$concuncbound is not squared.
    # The code in the line below calculate the square of load uncertainty at each time step with co-variance.
    flowcondf$loadbase = flowcondf$loadraw
    flowcondf$loaduncbound = (flowcondf$concuncbound + flowcondf$flowuncbound) * flowcondf$loadbase
  } else if (covar==0){
    flowcondf[, "concuncbound"] <- ((flowcondf[, "slope"]*flowcondf[, "conc_delta_i"])** 2 + ((1 - flowcondf[, "slope"])*flowcondf[, "conc_delta_i2"])** 2) / (flowcondf[, "conc"]**2)
    # The code in the line below calculate the square of load uncertainty at each time step without covariance.
    # Note that flowcondf$concuncbound is already the squared values.
    flowcondf[, "concuncbound"] <- round(flowcondf[, "concuncbound"], 3)
    flowcondf$loadbase = flowcondf$loadraw**2
    flowcondf$loaduncbound = (flowcondf$concuncbound + flowcondf$flowuncbound ** 2) * flowcondf$loadbase
  }
  else{
    flowcondf[, "concuncbound"] <- ((flowcondf[, "slope"]*flowcondf[, "conc_delta_i"])** 2 + ((1 - flowcondf[, "slope"])*flowcondf[, "conc_delta_i2"])** 2) / (flowcondf[, "conc"]**2)
    # The code in the line below calculate the square of load uncertainty at each time step without covariance.
    # Note that flowcondf$concuncbound is already the squared values.
    flowcondf[, "concuncbound"] <- round(flowcondf[, "concuncbound"], 3)
    flowcondf$loadbase = flowcondf$loadraw
    flowcondf$loaduncbound = sqrt(flowcondf$concuncbound + flowcondf$flowuncbound ** 2) * flowcondf$loadbase
    
  }
  # End if
  return(flowcondf)
}


loaduncannual <- function(flowcondf){
  loadunczoo <- read.zoo(flowcondf, header = TRUE)
  
  # aggregate to annual mean instant flow, concentration and load 
  yearlyunc <- aggregate(loadunczoo, as.year(time(loadunczoo)), sum)
  yearlyunc <- as.data.frame(yearlyunc)
  yearlyunc <- yearlyunc[-1,] 
  yearlyunc$loaduncbound <- round(yearlyunc$loaduncbound / yearlyunc$loadbase, 3)

  return(yearlyunc)
}

# Use a function to call the linear interpolation related functions
calloadunc <- function(flowdf, concdf, minconc, seconds, flowunc, concunc, yearlyloadbase, covar, boundtype){
  flowconczoo <- conclinear(flowdf, concdf, minconc)
  # browser()
  flowcondf <- loadunc(flowconczoo, concunc, flowunc, concdf, minconc, covar)
  yearlyloadunc <- loaduncannual(flowcondf)
  if (boundtype == "upper"){
    yearlyloadunc$load <- yearlyloadbase$load * (1 + yearlyloadunc$loaduncbound)
  } else {
    yearlyloadunc$load <- yearlyloadbase$load * (1 - yearlyloadunc$loaduncbound)
    }
  
  return(yearlyloadunc)
}


### Random sampling functions for sample frequency analysis
# random sample n, repeat N times. concentration data
randomall <- function(concdf, n, N){
  names(concdf) <- c("Datetime", "conc")
  randomout <- bind_rows(replicate(N, concdf %>% sample_n(n), simplify=F), .id="Obs")
  return(randomout)
}

# random sample n, repeat N times for a given hydrology year. concentration data
randomyear <- function(concdf, nprop, N){
  names(concdf) <- c("Datetime", "conc")
  concdf$hyear <- as.year(concdf$Datetime)
  set.seed(200)
  randomout <- bind_rows(replicate(N, concdf %>%
                                     group_by(hyear) %>%
                                     mutate(yearsamples = length(conc)) %>%
                                     slice_sample(prop=nprop, weight_by=yearsamples),
                                   simplify = FALSE),
                         .id="Obs")
  return(randomout)
}


### Estimating combined uncertainty with choice of random sample methods, data, and load estimation methods.
## Qian to change the settings of arguments and inside loop
samplefreqfun <- function(flowdf, concdf, seconds,minconc, loadfun, biasfun, randomfun, 
                          randomseeds, randomseedprop, randomtimes, flowunc, concunc, covar, boundtype, subsample_tf){
  # get random samples
  if (subsample_tf){
    if (randomfun == "randomall"){
      concsamplefreq <- randomall(concdf, randomseeds, randomtimes)
    } else if (randomfun == "randomyear"){
      concsamplefreq <- randomyear(concdf, randomseedprop, randomtimes)
    }
  } else {
    concsamplefreq <- concdf
    concsamplefreq$Obs <- 1
  }
    
  hyear <- as.year(flowdf$Datetime)
  
  conclist <- data.frame(matrix(0, randomtimes, length(unique(hyear))), row.names = seq(1, randomtimes))
  names(conclist) <- unique(hyear)
  
  if (loadfun == "linear"){
    for (i in (unique(concsamplefreq$Obs))){
      concsub <- subset(concsamplefreq, concsamplefreq$Obs == i)
      concsub <- as.data.frame(concsub[,c("DateTime", "conc")])
      # conclist[[i]] <- loadlinear(flowdf, concsub, seconds, minconc)
      # Calculate the new base load with re-sampled concentration
      # browser()
      flowconcbase_temp_day <- conclinear(flowdf, concsub, minconc)
      flowconcbase_temp_day$load <- flowconcbase_temp_day$conc * flowconcbase_temp_day$flow
      flowconcbase_temp_df <- data.frame(Date = as.Date(index(flowconcbase_temp_day)), coredata(flowconcbase_temp_day))
      dayflowconcbase <- aggregate(load ~ Date, data = flowconcbase_temp_df , FUN = mean)
      dayflowconcbase$load <- round(dayflowconcbase$load*86400 * 0.000001, 8)
      flowconcbase_temp <- conclinear(flowdf, concsub, minconc)
      yearlyloadbase_temp <- loadlinear(flowconcbase_temp, seconds)
      load_temp <- calloadunc(flowdf, concsub, minconc, seconds, flowunc, concunc, yearlyloadbase_temp, covar, boundtype) # Use Qian's new settings
      conclist[i, ] <- yearlyloadbase_temp$load
      # conclist[[i]]$hyear <- unique(hyear) ## TODO
    }
  } else if (loadfun == "beale"){
    for (i in (unique(concsamplefreq$Obs))){
      concsub <- subset(concsamplefreq, concsamplefreq$Obs == i)
      concsub <- as.data.frame(concsub[,c("Datetime", "conc")])
      conclist[[i]] <- loadbeale(flowdf, concsub, seconds, biasfun)
      conclist[[i]]$hyear <- unique(hyear) ## TODO
      dayflowconcbase <- NaN
    }
  }
  
  # concsampleall <- bind_rows(conclist, .id = "samples")
  # browser()
  # concsampleminmax <- concsampleall %>%
  #   group_by(hyear, flow) %>%
  #   dplyr::summarise(across(.cols = c("load"),list (min = min, max = max), na.rm=TRUE))
  # 
  return(list(var1 = conclist, var2 = dayflowconcbase))
  # return(list(var1 = yearlyloadbase_temp))
  
}


# Import data -------------------------------------------------------------
### Import and format flow data
# import flow data
flow <- read.csv("126001A_flows.csv", skip = 1, header = F)

# tidy flow data to date, flow, qc format
colnames(flow) <- rep(c("Datetime", "flow", "QC","X"), 10)
flow <- rbind(flow[,c(1:3)], flow[,c(5:7)], flow[,c(9:11)], flow[,c(13:15)], flow[,c(17:19)], 
              flow[,c(21:23)], flow[,c(25:27)], flow[,c(29:31)], flow[,c(33:35)])

flow <- na.omit(flow) 

flow$Datetime <- strptime(flow$Datetime, format="%d/%m/%Y %H:%M", tz = "Australia/Brisbane")

# Qian's edit
flowunc <- flow
flowunc$flowuncbound <- 0
  
# create flow upper bounds and lower bounds data
flowqc <- flow
flow <- flow[,-3]

# create flow upper bounds and lower bounds data
flowuncup <- mutate(flowunc, flowuncbound = case_when(
  QC == 9 ~ 0.05,
  QC == 10 ~ 0.1,
  QC == 20 ~ 0.2,
  QC == 30 ~ 0.5,
  QC == 60 ~ 0.75,
  TRUE ~ flowuncbound))

flowunclow <- mutate(flowunc, flowuncbound = case_when(
  QC == 9 ~ 0.05,
  QC == 10 ~ 0.1,
  QC == 20 ~ 0.2,
  QC == 30 ~ 0.5,
  QC == 60 ~ 0.75,
  TRUE ~ flowuncbound))
# Form the flowunc* into two columns (Datatime, flowuncbound)
flowuncup <- flowuncup[,c(1, 4)]
flowunclow <- flowunclow[,c(1, 4)]
flowuncbest <- flowunc[,c(1, 4)]
## For demonstration in graphical evaluation paper --> set flow unc to 0.0
flowuncup$flowuncbound <- 0
flowunclow$flowuncbound <- 0

## setting other input data
yearlyseconds <- c(31532400,31532400,31618800,31532400,31532400,31532400,31618800,31532400,31532400)

MinNH4 = 0.002
MinNOx = 0.001
MinDIN = MinNH4 + MinNOx


# ### scenarios data

# Estimate and plot combined uncertainty -----------------------------------------------------------
# Read input from new data provided by Danlu
### import and format concentration data
conc <- read.csv("126001A_din_concentrations_conditions.csv", header=TRUE, stringsAsFactors=FALSE) #TSS mg/L, NOx mg/L

# round the concentration data to nearest hour
conc$DateTime <- strptime(conc$DateTime, format="%d/%m/%Y %H:%M", tz = "Australia/Brisbane")
conc$DateTime <- round(conc$DateTime, units = "hours")
# Define the target date to filter by (convert to Date or POSIXct)
target_date <- as.POSIXct("2018-07-01", format="%Y-%m-%d")

# Filter the data frame for rows where 'Datetime' is before the target date
conc <- conc[conc$DateTime < target_date, ]

# convert NOx, NH4 and DIn data from character to numeric
conc["conc"] <- lapply(conc["conc"], function(x) as.numeric(x))
conc$concuncbound <- 0
# create separate TSS and Nox concentration data
# List of columns you want to keep
columns_to_keep <- c("DateTime", "conc")
dinraw<- conc[, columns_to_keep]
noxuncup <- conc[, c("DateTime", "concuncbound")]
# ## generating combined uncertainty for lower bound for din
# Generate random samples between lower and upper bounds
set.seed(123)
samples <- matrix(runif(nrow(dinraw)*100, min = 0.458, max = 1.629), nrow = nrow(dinraw), ncol = 100)
# Calculate base realisation without noise added.
din_copy <- dinraw
index <- 1:101
columns <- 2009:2018

# Create a sequence of dates from 2009-07-01 to 2018-06-30
date_sequence <- seq.Date(from = as.Date("2009-07-01"), to = as.Date("2018-06-30"), by = "day")
# Create an empty data frame with this date sequence as row names
dayload <- matrix(nrow = length(date_sequence), ncol = length(index))

dayload_zoo <- zoo(dayload, order.by = date_sequence)

# Create annual load dataframe
data <- matrix(0, nrow = length(index), ncol = length(columns))
dinlabunc <- as.data.frame(data)
colnames(dinlabunc) <- columns
# Set covar as TRUE
return_temp <- samplefreqfun(flowdf = flow, concdf = dinraw, seconds = yearlyseconds, minconc = MinDIN,
                         loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 1.0,
                         randomseeds = 10, randomtimes = 1, flowunc = flowuncup,
                         concunc = noxuncup, covar = TRUE, boundtype = "upper", subsample_tf = FALSE)
dinlabunc[1, ] <- return_temp$var1
dayload[, 1] <- (return_temp$var2)$load[-1]
# Add random noise to concentration and calculate annual loads.
for (ii in 1:ncol(samples)) {
  din_copy$conc <- dinraw$conc * samples[,ii]
  return_temp <- samplefreqfun(flowdf = flow, concdf = din_copy, seconds = yearlyseconds, minconc = MinDIN,
                           loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 1.0,
                           randomseeds = 10, randomtimes = 1, flowunc = flowuncup,
                           concunc = noxuncup, covar = TRUE, boundtype = "upper", subsample_tf = FALSE)
  dinlabunc[ii + 1, ] <- return_temp$var1
  dayload[, ii+1] <- (return_temp$var2)$load[-1]
}
# Save annual load unc data
write.csv(dinlabunc, 'output/resample_freq100_new_assumption/dinlabunc_inc_individual_sample.csv')
# Create the zoo object with date_sequence as the index
dayload_zoo <- zoo(dayload, order.by = date_sequence)
write.zoo(dayload_zoo, file = 'output/resample_freq100_new_assumption/dailyloadunc.csv', sep = ",", index.name = "Date")
