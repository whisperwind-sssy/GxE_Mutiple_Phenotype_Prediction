

capitalize <- function(string){
  index <- is.na(string)
  out <- rep(NA, length(string))
  out[!index] <- unlist(lapply(strsplit(as.character(string[!index]),""),function(x){
                    x[1] <- toupper(x[1])
                    paste(x, collapse="")
                  }))
  out
}

# Function to rename columns

# x=INFO[,CC[[i]][[1]]]; names_list=CC[[i]][2:3]
# pattern.replace=FALSE; check.repeated=TRUE; verbose=TRUE

#-------------------------------------------------------------------
#  replace:           Whether replace in 'x' only the substring that
#                     matches the pattern in names_list[k].
#                     Otherwise, if the pattern is found, the whole
#                     string in 'x' will be replaced
#  multiple.existing: When FALSE, will produce an error if the new
#                     name for the pattern exists more than once in 'x'
#  multiple.patterns: When FALSE, will produce an error if the pattern
#                     was found in 'x' more than once
#-------------------------------------------------------------------
Rename <- function(x, names_list,
                   multiple.existing=FALSE,
                   multiple.patterns=FALSE,
                   replace=c("pattern","all"),
                   verbose=TRUE) {

  replace <- match.arg(replace)
  out <- x[]
  unmatched <- c()  # To store names that were not found

  if(!is.list(names_list) | is.null(names(names_list))){
     stop("'names_list' must be list type with no null names(names_list)")
  }

  for(i in 1:length(names_list)){
    patterns <- names_list[[i]]
    new_names <- names(names_list)[i]
    index <- grep(patterns, out)   # Match the old_name
    if(length(index) > 0L){
      if(length(index) > 1L){
        if(length(grep("\\*",new_names)) == 1){
          new_names2 <- unlist(lapply(seq(length(index)),function(j)gsub("\\*",j,new_names)))
        }else{
          if(!multiple.patterns){
            stop("Multiple items were found for pattern(s) '",paste(patterns,collapse=","),
                 "'.\n       No replacement is performed")
          }
          new_names2 <- rep(new_names, length(index))
        }

        if(replace=="all"){
          # new_names will have multiple entries only when the whole word
          # will be replaced
          new_names <- new_names2
        }
      }

      if(replace=="pattern"){
        if(verbose){
          message("  Pattern '",patterns,"' in ",length(index)," item(s) of 'x'",
                  " was replaced by '",new_names,"'")
        }
        out[index] <- gsub(patterns, new_names, out[index])

      }else{
        for(j in 1:length(index)){
          if(!multiple.existing){
            if(sum(x %in% new_names[j]) > 1){
               stop("The new name '",new_names[j],"' was found multiple times in 'x'")
            }
          }
          if(out[index[j]]!=new_names[j]){
            if(verbose){
              message("  Entry '",out[index[j]],"' was replaced by '",new_names[j],"'")
            }
            out[index[j]] <- new_names[j]
          }
        }
      }

    }else{
      unmatched <- c(unmatched, paste(patterns, collapse=","))
    }
  }

  if((length(unmatched) > 0) & verbose){
     message("  Pattern(s) ",paste(unmatched,collapse=",")," could not be found")
  }
  return(out)
}

# ------------------------------------------------------------------------------

Rename2 <- function(x, names_list, check.repeated=FALSE,
                       pattern.replace=FALSE, verbose=TRUE) {

  out <- x[]
  unmatched <- c()  # To store names that were not found

  if(sum(unlist(lapply(names_list,length))>0) < 2){
     stop("'names_list[[",i,"]]' is not a list with 2 elements")
  }

  if(pattern.replace & length(names_list[[1]])>1){
     warning("Only one element is considered in 'names_list[[",i,"]]' when 'pattern.replace=TRUE'", immediate.=T)
  }
  index <- grep(names_list[[2]], out)   # match the old_name
  if(length(index) > 0){
    if(length(names_list[[1]]) != length(index)){
      if(check.repeated){
        stop("Multiple items '",paste(unique(out[index]),collapse=","),
           "' could not be matched to '", paste(names_list[[1]],collapse=","),"'")

      }else{
        if(length(names_list[[1]])==1){
          names_list[[1]] <- rep(names_list[[1]][1], length(index))
        }
      }
    }
    if(pattern.replace){
      if(verbose){
        cat("   Pattern '",names_list[[2]],"' was replaced by '",names_list[[1]][1],"'.\n", sep="")
      }
      out[index] <- gsub(names_list[[2]], names_list[[1]][1], out[index])  # gsub(pattern, replacement, x)

    }else{
      for(j in 1:length(index)){
        if(check.repeated){
          if(sum(out %in% names_list[[1]][j]) > 1){
             stop("The new name '",names_list[[1]][j],"' exists multiple times in 'x'")
          }
        }
        if(verbose){
          cat("   Field '",out[index[j]],"' was replaced by '",names_list[[1]][j],"'.\n", sep="")
        }
        out[index[j]] <- names_list[[1]][j]
      }
    }

  }else{
    unmatched <- c(unmatched, paste(names_list[[1]], collapse=","))
  }

  if((length(unmatched) > 0) & verbose){
     cat("   Names",paste(unmatched,collapse=","),"could not be found.\n")
  }
  return(out)
}

#=============================================================================

# Fix date format
fix_date <- function(date_as_character) {
  tmp <- strsplit(date_as_character, '/')
  nas <- sapply(tmp, length) != 3
  tmp2 <- tmp[!nas]

  d1 <- sapply(tmp2, function(x) x[1])
  d2 <- sapply(tmp2, function(x) x[2])
  d3 <- sapply(tmp2, function(x) x[3])
  d3[nchar(d3) == 2] <- paste0('20', d3[nchar(d3) == 2])
  swap <- which(as.numeric(d1) > 12)
  if(length(swap) > 0){
    if(any(as.numeric(d2) > 12))stop("Some values are greater than 12 simultaneously in days and months")
    d1_swap <- d1[swap]
    d1[swap] <- d2[swap]
    d2[swap] <- d1_swap
  }
  res_dates <- as.Date(paste0(d1,'/', d2, '/', d3), format = '%m/%d/%Y')
  res <- as.Date(rep(NA, length(date_as_character)))
  res[!nas] <- res_dates
  return(res)
}

#=============================================================================

# Calculate daily rainfall from a dataset with columns 'valid' and 'rainfall'

calculate_daily <- function(x) {
  require(tidyverse)
  options(dplyr.summarise.inform = FALSE) # added
  if (!all(c('rainfall', 'valid') %in% colnames(x))) stop('columns named (rainfall) and (valid) are needed.')
  precip <- x %>%
    filter(!is.na(rainfall)) %>%
    select(valid, rainfall) %>%
    mutate(date = with_tz(as.POSIXct(valid, tz='UCT'),
                          ""),
           date = date-dst(date)*3600,
           year = year(date),
           month = month(date),
           day = day(date),
           hour = hour(date),
           minute = minute(date)) %>%
    group_by(year, month, day, hour)

  modalminute <- precip %>%
    arrange(desc(minute), .by_group=TRUE) %>%
    mutate(maxminute = minute[which.max(rainfall)]) %>%
    {which.max(tabulate(.$maxminute))}

  prec <- precip %>%
    filter(minute <= modalminute) %>%
    summarize(hourlyprecip = max(rainfall, na.rm=TRUE)) %>%
    group_by(year, month, day) %>%
    summarize('rainfall' = sum(hourlyprecip, na.rm=TRUE))
  prec <- as.data.frame(prec)
  prec$date <- date(ISOdate(prec$year, prec$month, prec$day))
  prec[,5:4]
}

#=============================================================================

# Get weather data from ASOS/AWOS network through https://mesonet.agron.iastate.edu

getWeatherASOS <- function(time_period = NA, network = "IA_ASOS", sid = NA) {
  require(jsonlite)
  time_period <- as.Date(time_period)
  if (is.na(sid) | is.na(time_period)[1]) {
    uri <- paste("https://mesonet.agron.iastate.edu/geojson/network/", network, ".geojson", sep = "")
    data <- url(uri)
    jdict <- fromJSON(data)
    return(jdict$features$properties)
  } else {
    if(time_period[1] > time_period[2]) stop('error: time period 1 must be before 2')
    service <- "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"
    service <- paste(service, "data=all&tz=Etc/UTC&format=comma&latlon=yes&", sep = "")
    service <- paste(service, "year1=", year(time_period)[1], "&month1=", month(time_period)[1], "&day1=", mday(time_period)[1], "&", sep = "")
    service <- paste(service, "year2=", year(time_period)[2], "&month2=", month(time_period)[2], "&day2=", mday(time_period)[2], "&", sep = "")

    uri2 <- paste(service, "station=", sid, sep = '')
    data_table <- read.table(uri2, header = T, sep = ',')

    return(data_table)
  }
}

#=============================================================================

# Get weather station information from ASOS/AWOS

getStationASOS <- function(network){
  stations <- getWeatherASOS(network = network)
  if (is.null(stations[1])) stop('network not found')
  stations$lat <- NA
  stations$lon <- NA
  for (i in 1:nrow(stations)){
    tmp <- getWeatherASOS(time_period = c('2018-10-10', '2018-10-10'), network, stations$sid[i])
    stations$lat[i] <- tmp$lat[1]
    stations$lon[i] <- tmp$lon[1]
  }
  stations[,c('elevation', 'sname', 'county', 'state', 'country', 'sid', 'lat', 'lon')]
}

#=============================================================================

# Get weather data from NOAA

getWeatherNOAA <- function(time_period = c('2018-01-01', '2018-12-31'), sid = 'ZI000067983', dataid = 'GHCND') {
  sid2 <- paste0(dataid, ':', sid)
  weather <- list()
  tmp_offset <- 0
  while(T){
    tmp <- ncdc(datasetid = dataid, stationid = sid2, startdate = time_period[1], enddate = time_period[2], limit = 1000, offset = tmp_offset)$data
    tmp_offset <- tmp_offset + 1000
    weather[[length(weather) + 1]] <- tmp
    if (nrow(tmp) < 1000) break
  }
  weather <- as.data.frame(do.call(rbind, weather))
  dtypes <- grep('PRCP|TMAX|TMIN|TOBS', unique(weather$datatype), value = T)
  dweather <- data.frame(date = unique(weather$date))
  for (ty in dtypes){
    tmpw <- weather[weather$datatype == ty, c('date', 'value')]
    dweather[,ty] <- tmpw$value[match(dweather$date, tmpw$date)] / 10
  }
  return(dweather)
}

#=============================================================================

# Function to calculate GDD from weather data

# type=c("C","F")[1]; Tbase=10; Tmax=30; date_names=c("anthesis","silking")
get_GDD_old <- function(dat,type=c("C","F"),Tbase=10,Tmax=30,date_names=c("anthesis","silking"))
{
  type <- match.arg(type)

  end_date <- NA
  if(any(!is.na(dat$harvesting))) end_date <- max(as.Date(dat$harvesting),na.rm=T)
  if(is.na(end_date)) end_date <- max(as.Date(dat$silking),as.Date(dat$anthesis),na.rm=T) + 100
  if(is.na(end_date)) end_date <- min(as.Date(dat$sowing)) + 500

  pwr <- data.frame(nasapower::get_power(community = 'AG',
            lonlat = c(dat$lon[1], dat$lat[1]),
            pars = c("T2M","T2M_MIN","T2M_MAX"),
            dates = c(min(as.Date(dat$sowing))-5, end_date),
            temporal_api = 'daily'))
  pwr$Tmin <- ifelse(pwr$T2M_MIN < Tbase, Tbase, pwr$T2M_MIN)
  pwr$Tmax <- ifelse(pwr$T2M_MAX > Tmax, Tmax, pwr$T2M_MAX)
  pwr$Tavg <- (pwr$Tmin + pwr$Tmax)/2
  pwr$GDD <- ifelse(pwr$Tavg<Tbase, Tbase, pwr$Tavg-Tbase)

  out <- matrix(NA,ncol=length(date_names),nrow=nrow(dat))
  colnames(out) <- paste0(date_names,"_GDD",type)
  for(i in 1:nrow(dat)){
      isowing <- which(pwr$YYYYMMDD == dat$sowing[i])
      for(k in 1:length(date_names)){
        if(!is.na(dat[i,date_names[k]])){
          idate <- which(pwr$YYYYMMDD == dat[i,date_names[k]])
          out[i,k] <- sum(pwr[isowing:idate,"GDD"])
        }
      }
  }
  out
}

# Get the GDD from a met file
# Tbase=0; Tmax=34; Tmax2=44; method=c("Simple","Jones")[2]
get_GDD <- function(met, data, sowing=NULL,
            date_names=c("date_silking","date_harvest"),
            Tbase=10, Tmax=34, Tmax2=44,
            method=c("Simple","Jones"))
{
  method <- match.arg(method)

  # Temperture factor: Jones, C.A., Kiniry, J.R., Dyke, P.T., 1986. (pp51)
  # "CERES-Maize: a simulation model of maize growth and development"
  tmfact_fn <- function(i){
    0.931 + 0.114*i - 0.0703*(i^2) + 0.0053*(i^3)
  }

  if(any(!date_names %in% colnames(data))){
    stop("Some 'date_names' were not found in column names of 'data'")
  }

  # Convert days to date yyyy-mm-dd
  met <- data.frame(date=as.Date(paste0(met$year,"-01-01"))+met$day-1, met)
  stopifnot(all(lubridate::yday(met$date) == met$day))

  stopifnot(all(unique(as.character(data[,date_names])) %in% as.character(met$date)))

  if(method=="Simple"){
    # If the Max and/or Min Temp < Tbase, it's set equal to Tbase
    met$mint <- ifelse(met$mint < Tbase, Tbase, met$mint)

    # If the daily Max Temp < Tbase, it's set equal to Tbase
    # If the daily Max Temp > Tmax, it's set equal to Tmax
    met$maxt <- ifelse(met$maxt < Tbase, Tbase,
                ifelse(met$maxt > Tmax, Tmax, met$maxt))

    met$avgt <- (met$mint + met$maxt)/2

    met$TT <- ifelse(met$avgt < Tbase, Tbase, met$avgt-Tbase)
  }
  if(method=="Jones"){
    tempfactor <- tmfact_fn(1:8)
    met$avgt <- (met$mint + met$maxt)/2
    met$TT <- NA

    #pp <- data.frame(x=c(0,18,26,34,44),y=c(0,10,18,26,0),cut=NA,m=NA)
    #for(k in 1:(nrow(pp)-1)){
    #  pp[k+1,'m'] <- (pp$y[k+1]-pp$y[k])/(pp$x[k+1]-pp$x[k])
    #  pp[k+1,'cut'] <- paste0(ifelse(k==1,"[","("),pp$x[k],",",pp$x[k+1],"]")
    #}
    pp <- data.frame(x=c(0,18,34,44),y=c(0,10,26,0),cut=NA,m=NA)
    for(k in 1:(nrow(pp)-1)){
      pp[k,'m'] <- (pp$y[k+1]-pp$y[k])/(pp$x[k+1]-pp$x[k])
      pp[k,'cut'] <- paste0("[",pp$x[k],",",pp$x[k+1],")")
    }
    breaks0 <- pp$x
    pp <- pp[-nrow(pp),]

    for(i in 1:nrow(met)){
       # If Tbase<=tmin[i]<=Tmax & Tbase<=tmax[i]<=Tmax
       flag1 <- (Tbase<=met$mint[i])&(met$mint[i]<=Tmax)
       flag2 <- (Tbase<=met$maxt[i])&(met$maxt[i]<=Tmax)
       if(flag1&flag2){
         met$TT[i] <- met$avgt[i]-Tbase
       }else{
         if((met$mint[i] < Tbase) | (met$maxt[i] > Tmax))
         {
           tmp <- met$mint[i] + tempfactor*(met$maxt[i] - met$mint[i])

           ## Page 66-67 in Jones 1986
           tt1 <- ifelse(Tbase<=tmp & tmp<=Tmax,
                  tmp-Tbase,
                  ifelse(Tmax < tmp & tmp<=Tmax2,
                         (Tmax-Tbase)*(1-(tmp-Tmax)/10),
                         ifelse(tmp<Tbase | tmp>Tmax2,
                         0,
                         NA # This value should not appear
                       )))
           stopifnot(all(!is.na(tt1)))
           br <- cut(tt1, breaks=breaks0, include.lowest=FALSE,right=FALSE)
           stopifnot(all(!is.na(br)))
           stopifnot(all(levels(br) == pp$cut))

           tt2 <- unlist(lapply(1:length(tt1),function(k){
              p0 <- pp[which(pp$cut==br[k]),]
              p0$m*(tt1[k]-p0$x) + p0$y
           }))
           met$TT[i] <- mean(tt2)
         }else{
           stop("Thermal time could not be calculated at row ",i," of met file")
         }
      }
    }
  }

  out <- matrix(NA,ncol=length(date_names),nrow=nrow(data))
  colnames(out) <- paste0(date_names,"_GDD")

  if(is.null(sowing)){
    for(i in 1:nrow(data)){
      for(k in 1:length(date_names)){
        idate <- which(as.character(met$date) == as.character(data[i,date_names[k]]))
        out[i,k] <- met$TT[idate]
      }
    }
  }else{
    stopifnot(all(as.character(data[,sowing]) %in% as.character(met$date)))
    for(i in 1:nrow(data)){
      isowing <- which(as.character(met$date) == as.character(data[i,sowing]))
      for(k in 1:length(date_names)){
        idate <- which(as.character(met$date) == as.character(data[i,date_names[k]]))
        out[i,k] <- sum(met$TT[isowing:idate])
      }
    }
  }

  out
}


#=============================================================================

calcGDD <- function(wdata, phenotype, envCol = 'env', intCol = c('date_plant', 'date_silking'), basetemp = 10) {
  # works with celsius degrees only
  if (!(envCol %in% colnames(wdata) & envCol %in% colnames(phenotype))) stop('envCol must be in colnames of wdata and phenotype')

  phenotype <- phenotype[order(phenotype$env),]
  phenotype$GDD <- NA
  phenotype[,intCol[1]] <- as.Date(phenotype[,intCol[1]])
  phenotype[,intCol[2]] <- as.Date(phenotype[,intCol[2]])
  # Remove missing dates and environments
  phenotype <- phenotype[!is.na(phenotype[,intCol[1]]),]
  phenotype <- phenotype[!is.na(phenotype[,intCol[2]]),]
  phenotype <- phenotype[phenotype[,envCol] %in% unique(wdata[,envCol]),]

  for (env in unique(phenotype$env)) {
    envi <- wdata[wdata$env == env, c('date', 'temp')]
    if (nrow(envi) > 0){
      # missing mean temperature is filled in with the environment's mean temp
      envi$temp[is.na(envi$temp)] <- mean(envi$temp, na.rm = T)
      envi$temp[envi$temp > 30] <- 30
      envi$temp[envi$temp < basetemp] <- basetemp
      envi$gdd <- envi$temp - basetemp
    }
    for (i in which(phenotype$env == env)){
      envrows <- match(as.Date(phenotype[i, intCol[1]]:phenotype[i, intCol[2]], '1970-1-1'), as.Date(envi$date))
      phenotype$GDD[i] <- sum(envi$gdd[envrows], na.rm = T)
    }
  }
  return(phenotype)
}

#=============================================================================

search_apsimx <- function(tmp, keyword, return_lvls = FALSE) {
  lvls <- vector()
  for (i in 1:5) {
    sapp <- sapply(tmp$Children, function(x) length(grep(keyword, unlist(capture.output(str(x))))) > 0)
    if (length(sapp) == 0) break
    if (length(sapp) == 1) if (!sapp) break
    lvls[i] <- which(sapp)
    tmp <- tmp$Children[[lvls[i]]]
  }
  if (return_lvls) {
    return(lvls)
  } else {
    return(tmp)
  }
}

#=============================================================================
# New version of edit APSIM
# file=simfile; src.dir=simdir; node = 'Cultivar'; overwrite = T; verbose = F; parm = 'Name'; value = 'Custom'
# wrt.dir=NULL; soil.child="Metadata"; manager.child = NULL;  edit.tag = "-edited"; parm.path = NULL
edit_apsimx_new <- function (file, src.dir = ".", wrt.dir = NULL,
                             node = c("Clock", 'Report', 'Cultivar',
                                      "Weather", "Soil", "SurfaceOrganicMatter",
                                      "MicroClimate", "Crop", "Manager", "Other"),
                             soil.child = c("Metadata", "Water", "SoilWater",
                                            "Organic", "Physical", "Analysis",
                                            "Chemical", "InitialWater", "Sample"),
                             manager.child = NULL, parm = NULL, value = NULL, overwrite = FALSE,
                             edit.tag = "-edited", parm.path = NULL, root, verbose = TRUE)
{
  #.check_apsim_name(file)
  if (missing(wrt.dir))
    wrt.dir <- src.dir
  file.names <- dir(path = src.dir, pattern = ".apsimx$",
                    ignore.case = TRUE)
  if (length(file.names) == 0) {
    stop("There are no .apsimx files in the specified directory to edit.")
  }
  node <- match.arg(node)
  soil.child <- match.arg(soil.child)
  edited.child <- "none"
  file <- match.arg(file, file.names)
  if (apsimx_filetype(file = file, src.dir = src.dir) != "json")
    stop("This function only edits JSON files")
  apsimx_json <- jsonlite::read_json(paste0(src.dir, "/",
                                            file))
  wcore <- grep("Core.Simulation", apsimx_json$Children)
  if (length(wcore) > 1) {
    if (missing(root)) {
      cat("Simulation structure: \n")
      str_list(apsimx_json)
      stop("more than one simulation found and no root node label has been specified \n select one of the children names above")
    }
    else {
      if (length(root) == 1) {
        wcore1 <- grep(as.character(root), apsimx_json$Children)
        if (length(wcore1) == 0 || length(wcore1) > 1)
          stop("no root node label found or root is not unique")
        parent.node <- apsimx_json$Children[[wcore1]]$Children
      }
      else {
        root.node.0.names <- sapply(apsimx_json$Children,
                                    function(x) x$Name)
        wcore1 <- grep(as.character(root[1]), root.node.0.names)
        root.node.0 <- apsimx_json$Children[[wcore1]]
        root.node.0.child.names <- sapply(root.node.0$Children,
                                          function(x) x$Name)
        wcore2 <- grep(as.character(root[2]), root.node.0.child.names)
        parent.node <- apsimx_json$Children[[wcore1]]$Children[[wcore2]]$Children
      }
    }
  }else {
    parent.node <- apsimx_json$Children[[wcore]]$Children
  }

  if (node == "Clock") {
    parm.choices <- c("Start", "End")
    parm <- match.arg(parm, choices = parm.choices, several.ok = TRUE)
    wlc <- function(x) grepl("Clock", x$Name)
    wlcl <- sapply(parent.node, FUN = wlc)
    start <- grep("Start", names(parent.node[wlcl][[1]]),
                  ignore.case = TRUE, value = TRUE)
    end <- grep("End", names(parent.node[wlcl][[1]]),
                ignore.case = TRUE, value = TRUE)
    if (length(parm) == 1) {
      if (parm == "Start") {
        parent.node[wlcl][[1]][start] <- value
      }
      if (parm == "End") {
        parent.node[wlcl][[1]][end] <- value
      }
    }
    if (length(parm) == 2) {
      if (parm[1] == "Start") {
        parent.node[wlcl][[1]][start] <- value[1]
      }
      if (parm[2] == "End") {
        parent.node[wlcl][[1]][end] <- value[2]
      }
    }
    apsimx_json$Children[[1]]$Children <- parent.node
  }
  if (node == "Weather") {
    wlw <- function(x) grepl("Weather", x$Name)
    wlwl <- sapply(parent.node, FUN = wlw)
    parent.node[wlwl][[1]]$FileName <- value
  }
  wcz <- grepl("Models.Core.Zone", parent.node)
  core.zone.node <- parent.node[wcz][[1]]$Children
  if (node == "Soil") {
    wsn <- grepl("Models.Soils.Soil", core.zone.node)
    soil.node <- core.zone.node[wsn]
    soil.node0 <- soil.node[[1]]$Children
    if (soil.child == "Metadata") {
      edited.child <- soil.child
      metadata.parms <- c("RecordNumber", "ASCOrder",
                          "ASCSubOrder", "SoilType", "LocalName",
                          "Site", "NearestTown", "Region",
                          "State", "Country", "NaturalVegetation",
                          "ApsoilNumber", "Latitude", "Longitude",
                          "LocationAccuracy", "DataSource",
                          "Comments")
      if (!all(parm %in% metadata.parms))
        stop("parm name(s) might be wrong")
      for (i in seq_along(parm)) {
        soil.node[[1]][[parm[i]]] <- value[i]
      }
    }
    if (soil.child == "Water") {
      edited.child <- soil.child
      wwn <- grep("^Water", sapply(soil.node[[1]]$Children, function(x) x$Name))
      soil.water.node <- soil.node[[1]]$Children[[wwn]]
      if (soil.water.node$Name != "Water") {
        stop("Wrong node (Soil Water)")
      }
      crop.parms <- c("XF", "KL", "LL")
      if (parm %in% crop.parms) {
        for (i in 1:length(soil.water.node$Children[[1]][[parm]])) {
          soil.water.node$Children[[1]][[parm]][[i]] <- value[i]
        }
      }
      else {
        for (i in 1:length(soil.water.node[[parm]])) {
          soil.water.node[[parm]][[i]] <- value[i]
        }
      }
      soil.node[[1]]$Children[[wwn]] <- soil.water.node
    }
    if (soil.child == "Physical") {
      wpn <- grep("^Physical", sapply(soil.node[[1]]$Children, function(x) x$Name))
      if (!is.list(value)) value <- as.list(value)
      depth_length <- length(soil.node[[1]]$Children[[wpn]]$Depth)
      if (length(value) != depth_length) value <- value[1:depth_length]
      soil.node[[1]]$Children[[wpn]][[parm]] <- value
    }
    if (soil.child == "SoilWater") {
      edited.child <- soil.child
      wswn <- grep("^SoilWater", sapply(soil.node[[1]]$Children,
                                        function(x) x$Name))
      soil.soilwater.node <- soil.node[[1]]$Children[[wswn]]
      soilwat.parms <- c("SummerDate", "SummerU",
                         "SummerCona", "WinterDate", "WinterU",
                         "WinterCona", "DiffusConst", "DiffusSlope",
                         "Salb", "CN2Bare", "CNRed",
                         "CNCov", "Slope", "DischargeWidth",
                         "CatchmentArea")
      if (parm %in% soilwat.parms) {
        for (i in seq_along(parm)) {
          soil.soilwater.node[[parm[i]]] <- value[i]
        }
      }
      else {
        if (!parm %in% c("SWCON", "KLAT"))
          stop("parameter is likely incorrect")
        for (i in 1:length(soil.soilwater.node[[parm]])) {
          soil.soilwater.node[[parm]][[i]] <- value[i]
        }
      }
      soil.node[[1]]$Children[[wswn]] <- soil.soilwater.node
    }
    if (soil.child == "Nitrogen") {
      wnn <- grepl("Nitrogen", soil.node0)
      soil.nitrogen.node <- soil.node0[wnn][[1]]
      for (i in 1:length(soil.nitrogen.node[[parm]])) {
        soil.nitrogen.node[[parm]][[i]] <- value[i]
      }
      soil.node[[1]]$Children[wnn][[1]] <- soil.nitrogen.node
    }
    if (soil.child == "Organic") {
      edited.child <- "Organic"
      wsomn <- grepl("Organic", soil.node0)
      soil.om.node <- soil.node0[wsomn][[1]]
      som.parms1 <- c("RootCN", "EnrACoeff",
                      "EnrBCoeff")
      if (parm %in% som.parms1) {
        soil.om.node[[parm]] <- value
      }
      else {
        for (i in 1:length(soil.om.node[[parm]])) {
          soil.om.node[[parm]][[i]] <- value[i]
        }
      }
      soil.node[[1]]$Children[wsomn][[1]] <- soil.om.node
    }
    if (soil.child == "Analysis" || soil.child == "Chemical") {
      edited.child <- soil.child
      wan <- grepl(soil.child, soil.node0)
      soil.analysis.node <- soil.node0[wan][[1]]
      #if (parm != "PH")
      #  stop("only PH can be edited, use 'edit_apsimx_replace_soil_profile instead")
      #if (parm == "PH") {
      for (i in 1:length(soil.analysis.node[[parm]])) {
        soil.analysis.node[[parm]][[i]] <- value[i]
        #  }
      }
      soil.node[[1]]$Children[wan][[1]] <- soil.analysis.node
    }
    if (soil.child == "InitialWater") {
      edited.child <- "InitialWater"
      wiwn <- grepl("InitialWater", soil.node0)
      soil.initialwater.node <- soil.node0[wiwn][[1]]
      siw.parms <- c("PercentMethod", "FractionFull",
                     "DepthWetSoil")
      parm <- match.arg(parm, choices = siw.parms)
      soil.initialwater.node[[parm]] <- value
      soil.node[[1]]$Children[wiwn][[1]] <- soil.initialwater.node
    }
    if (soil.child == "Sample") {
      edited.child <- "Sample"
      wsn <- grepl("Sample", soil.node0)
      soil.sample.node <- soil.node0[wsn][[1]]
      for (i in 1:length(soil.sample.node[[parm]])) {
        soil.sample.node[[parm]][[i]] <- value[i]
      }
      soil.node[[1]]$Children[wsn][[1]] <- soil.sample.node
    }
    core.zone.node[wsn] <- soil.node
  }
  if (node == "SurfaceOrganicMatter") {
    wsomn <- grepl("Models.Surface.SurfaceOrganicMatter",
                   core.zone.node)
    som.node <- core.zone.node[wsomn][[1]]
    if (som.node$Name != "SurfaceOrganicMatter") {
      stop("Wrong node")
    }
    som.node[[parm]] <- value
    core.zone.node[wsomn][[1]] <- som.node
  }
  if (node == "MicroClimate") {
    wmcn <- grepl("Models.MicroClimate", core.zone.node)
    microclimate.node <- core.zone.node[wmcn][[1]]
    if (microclimate.node$Name != "MicroClimate") {
      stop("Wrong node")
    }
    microclimate.node[[parm]] <- value
    core.zone.node[wmcn][[1]] <- microclimate.node
  }
  if (node == "Crop") {
    wmmn <- grepl("Models.Manager", core.zone.node)
    manager.node <- core.zone.node[wmmn]
    wcn <- grepl("CultivarName", manager.node)
    crop.node <- manager.node[wcn][[1]]$Parameters
    for (i in 1:length(crop.node)) {
      if (crop.node[[i]]$Key == parm) {
        crop.node[[i]]$Value <- value
      }
    }
    core.zone.node[wmmn][wcn][[1]]$Parameters <- crop.node
  }
  if (node == "Other") {
    upp <- strsplit(parm.path, ".", fixed = TRUE)[[1]]
    upp.lngth <- length(upp)
    if (upp.lngth < 5)
      stop("Parameter path too short?")
    if (upp.lngth > 10)
      stop("Cannot handle this yet")
    if (apsimx_json$Name != upp[2])
      stop("Simulation root name does not match")
    wl3 <- which(upp[3] == sapply(apsimx_json$Children, function(x) x$Name))
    if (length(wl3) == 0)
      stop("Could not find parameter at level 3")
    n3 <- apsimx_json$Children[[wl3]]
    wl4 <- which(upp[4] == sapply(n3$Children, function(x) x$Name))
    if (length(wl4) == 0)
      stop("Could not find parameter at level 4")
    if (upp.lngth == 5) {
      if (upp[5] %in% names(n3$Children[[wl4]])) {
        apsimx_json$Children[[wl3]]$Children[[wl4]][[upp[5]]] <- value
      }
      else {
        wl5 <- which(upp[5] == sapply(n3$Children[[wl4]]$Children,
                                      function(x) x$Name))
        if (length(wl5) == 0)
          stop("Could not find parameter at level 5")
        apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]][[upp[5]]] <- value
      }
    }
    if (upp.lngth == 6) {
      n4 <- apsimx_json$Children[[wl3]]$Children[[wl4]]
      wl5 <- which(upp[5] == sapply(n4$Children, function(x) x$Name))
      if (length(wl5) == 0)
        stop("Could not find parameter at level 5")
      if (upp[6] %in% names(n4$Children[[wl5]])) {
        apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]][[upp[6]]] <- value
      }
      else {
        if ("Parameters" %in% names(n4$Children[[wl5]])) {
          wp <- grep(upp[6], n4$Children[[wl5]]$Parameters)
          if (length(wp) == 0)
            stop("Could not find parameter at level 6 (Parameter)")
          apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Parameters[[wp]]$Value <- value
        }
        else {
          wl6 <- which(upp[6] == sapply(n4$Children[[wl5]]$Children,
                                        function(x) x$Name))
          if (length(wl6) == 0)
            stop("Could not find parameter at level 6")
          apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]][[upp[6]]] <- value
        }
      }
    }
    if (upp.lngth == 7) {
      n4 <- apsimx_json$Children[[wl3]]$Children[[wl4]]
      wl5 <- grep(upp[5], n4$Children)
      n5 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]
      wl6 <- grep(upp[6], n4$Children)
      apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]][[upp[7]]] <- value
    }
    if (upp.lngth == 8) {
      n4 <- apsimx_json$Children[[wl3]]$Children[[wl4]]
      wl5 <- grep(upp[5], n4$Children)
      n5 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]
      wl6 <- grep(upp[6], n5$Children)
      n6 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]
      wl7 <- grep(upp[7], n6$Children)
      apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]$Children[[wl7]][[upp[8]]] <- value
    }
    if (upp.lngth == 9) {
      n4 <- apsimx_json$Children[[wl3]]$Children[[wl4]]
      wl5 <- grep(upp[5], n4$Children)
      n5 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]
      wl6 <- grep(upp[6], n5$Children)
      n6 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]
      wl7 <- grep(upp[7], n6$Children)
      n7 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]$Children[[wl7]]
      wl8 <- grep(upp[8], n7$Children)
      apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]$Children[[wl7]]$Children[[wl8]][[upp[9]]] <- value
    }
    if (upp.lngth == 10) {
      n4 <- apsimx_json$Children[[wl3]]$Children[[wl4]]
      wl5 <- grep(upp[5], n4$Children)
      n5 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]
      wl6 <- grep(upp[6], n5$Children)
      n6 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]
      wl7 <- grep(upp[7], n6$Children)
      n7 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]$Children[[wl7]]
      wl8 <- grep(upp[8], n7$Children)
      n8 <- apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]$Children[[wl7]]$Children[[wl8]]
      wl9 <- grep(upp[9], n8$Children)
      apsimx_json$Children[[wl3]]$Children[[wl4]]$Children[[wl5]]$Children[[wl6]]$Children[[wl7]]$Children[[wl8]]$Children[[wl9]][[upp[10]]] <- value
    }
  }
  if (node == 'Cultivar') {
    lvls <- search_apsimx(apsimx_json, 'Models.PMF.Plant, Models', T)

    if(length(grep('Models.PMF.Cultivar, Models', capture.output(str(apsimx_json)), ignore.case = T)) == 0) {
      x <- list(list())
      x[[1]]$`$type` <- 'Models.PMF.Cultivar, Models'
      x[[1]]$Command <- list('[Phenology].Juvenile.Target.FixedValue = 110',
                             '[Phenology].GrainFilling.Target.FixedValue = 550')
      x[[1]]$Name <- value
      x[[1]]$Children <- list()
      x[[1]]$IncludeInDocumentation <- T
      x[[1]]$Enabled <- T
      x[[1]]$ReadOnly <- F
      apsimx_json$Children[[lvls[1]]]$Children[[lvls[2]]]$Children[[lvls[3]]]$Children <- x
      parm <- 'New edited cultivar created'
      value <- 'Example values added'
    } else {
      # Edit existing cultivar
      if (parm == 'Command') if (!is.list(value)) value <- as.list(value)
      apsimx_json$Children[[lvls[1]]]$Children[[lvls[2]]]$Children[[lvls[3]]]$Children[[1]][[parm]] <- value
    }
    node <- "Other"
    parm.path <- ""
    print.path <- F
  }
  if (node == "Manager") {
    wmmn <- grepl("Models.Manager", core.zone.node)
    manager.node <- core.zone.node[wmmn]
    manager.node.names <- sapply(manager.node, FUN = function(x) x$Name)
    if (missing(manager.child))
      stop("need to specify manager.child")
    edited.child <- manager.child
    wmc <- grep(manager.child, manager.node.names)
    if (length(wmc) == 0) {
      manager.node.names <- sapply(manager.node[[1]]$Children,
                                   FUN = function(x) x$Name)
      wmc2 <- grep(manager.child, manager.node.names)
      manager.child.node <- manager.node[[1]]$Children[[wmc2]]$Parameters
    }
    else {
      manager.child.node <- manager.node[[wmc]]$Parameters
    }
    for (i in 1:length(manager.child.node)) {
      if (manager.child.node[[i]]$Key == parm) {
        manager.child.node[[i]]$Value <- value
      }
    }
    if (length(wmc) == 0) {
      manager.node[[1]]$Children[[wmc2]]$Parameters <- manager.child.node
    }
    else {
      manager.node[[wmc]]$Parameters <- manager.child.node
    }
    core.zone.node[wmmn] <- manager.node
  }

  if (node == "Report") {
    if (!parm %in% c("VariableNames", "EventNames"))
      stop ('When node = "Report", parm must be either "VariableNames" or "EventNames".')
    if (class(value) != 'list') {
      warning('value should be a list. It will be coerced.')
      value <- as.list(value)
    }
    tmp <- search_apsimx(apsimx_json, parm)
    lvls <- search_apsimx(apsimx_json, parm, T)

    if (length(lvls) == 2)
      apsimx_json$Children[[lvls[1]]]$Children[[lvls[2]]][[parm]] <- value
    if (length(lvls) == 3)
      apsimx_json$Children[[lvls[1]]]$Children[[lvls[2]]]$Children[[lvls[3]]][[parm]] <- value
    if (length(lvls) == 4)
      apsimx_json$Children[[lvls[1]]]$Children[[lvls[2]]]$Children[[lvls[3]]]$Children[[lvls[4]]][[parm]] <- value
    if (length(lvls) == 5)
      apsimx_json$Children[[lvls[1]]]$Children[[lvls[2]]]$Children[[lvls[3]]]$Children[[lvls[4]]]$Children[[lvls[5]]][[parm]] <- value

    node <- "Other"
    parm.path <- ""
    print.path <- F
  }

  if (node != "Other") {
    parent.node[wcz][[1]]$Children <- core.zone.node
    if (length(wcore) > 1) {
      if (length(root) == 1) {
        apsimx_json$Children[[wcore1]]$Children <- parent.node
      }
      else {
        apsimx_json$Children[[wcore1]]$Children[[wcore2]]$Children <- parent.node
      }
    }
    else {
      apsimx_json$Children[[wcore]]$Children <- parent.node
    }
  }
  if (overwrite == FALSE) {
    wr.path <- paste0(wrt.dir, "/", tools::file_path_sans_ext(file),
                      edit.tag, ".apsimx")
  }
  else {
    wr.path <- paste0(wrt.dir, "/", file)
  }
  jsonlite::write_json(apsimx_json, path = wr.path, pretty = TRUE,
                       digits = NA, auto_unbox = TRUE, null = "null")
  if (verbose) {
    cat("Edited (node): ", node, "\n")
    cat("Edited (child): ", edited.child, "\n")
    cat("Edited parameters: ", parm, "\n")
    cat("New values: ", unlist(value), "\n")
    cat("Created: ", wr.path, "\n")
  }
}

#=============================================================================

# New inspect_apsimx

inspect_apsimx_new <- function (file = "", src.dir = ".",
                                node = c("Clock","Weather", "Soil", "SurfaceOrganicMatter",
                                         "Cultivar", "MicroClimate", "Crop", "Manager", "Other", "Report"),
                                soil.child = c("Metadata", "Water", "InitialWater",
                                               "Chemical", "Physical", "Analysis",
                                               "SoilWater", "InitialN", "CERESSoilTemperature",
                                               "Sample", "Nutrient", "Organic"), parm = NULL,
                                digits = 3, print.path = FALSE, root)
{
  apsimx:::.check_apsim_name(file)
  file.names <- dir(path = src.dir, pattern = ".apsimx$",
                    ignore.case = TRUE)
  if (length(file.names) == 0) {
    stop("There are no .apsimx files in the specified directory to inspect.")
  }
  node <- match.arg(node)
  soil.child <- match.arg(soil.child)
  if (soil.child %in% c("Nutrient"))
    stop("Not implemented yet")
  file <- match.arg(file, file.names)
  apsimx_json <- jsonlite::read_json(paste0(src.dir, "/",
                                            file))
  parm.path.0 <- paste0(".", apsimx_json$Name)
  fcsn <- grep("Models.Core.Simulation", apsimx_json$Children,
               fixed = TRUE)
  if (length(fcsn) > 1) {
    if (missing(root)) {
      cat("Simulation structure: \n")
      str_list(apsimx_json)
      stop("more than one simulation found and no root node label has been specified \n select one of the children names above")
    }
    else {
      if (length(root) == 1) {
        nms <- vapply(apsimx_json$Children, FUN = function(x) x$Name,
                      FUN.VALUE = "character")
        fcsn <- grep(as.character(root), nms)
        parm.path.1 <- paste0(parm.path.0, ".",
                              apsimx_json$Children[[fcsn]]$Name)
        parent.node <- apsimx_json$Children[[fcsn]]$Children
        if (length(fcsn) == 0 || length(fcsn) > 1)
          stop("no root node label found or root is not unique")
      }
      else {
        nms1 <- vapply(apsimx_json$Children, FUN = function(x) x$Name,
                       FUN.VALUE = "character")
        fcsn1 <- grep(as.character(root[1]), nms1)
        root.node.0 <- apsimx_json$Children[[fcsn1]]
        root.node.0.child.names <- vapply(root.node.0$Children,
                                          function(x) x$Name, FUN.VALUE = "character")
        fcsn2 <- grep(as.character(root[2]), root.node.0.child.names)
        parent.node <- apsimx_json$Children[[fcsn1]]$Children[[fcsn2]]$Children
        parm.path.1 <- paste0(parm.path.0, ".",
                              apsimx_json$Children[[fcsn1]]$Children[[fcsn2]])
      }
    }
  }
  else {
    parent.node <- apsimx_json$Children[[fcsn]]$Children
    parm.path.1 <- paste0(parm.path.0, ".", apsimx_json$Children[[fcsn]]$Name)
  }

  if (node == "Clock") {
    wlc <- function(x) grepl("Clock", x$Name, ignore.case = TRUE)
    wlcl <- sapply(parent.node, FUN = wlc)
    clock.node <- as.list(parent.node[wlcl])[[1]]
    start.name <- grep("start", names(clock.node),
                       ignore.case = TRUE, value = TRUE)
    end.name <- grep("end", names(clock.node), ignore.case = TRUE,
                     value = TRUE)
    cat("Start:", clock.node[[start.name]], "\n")
    cat("End:", clock.node[[end.name]], "\n")
    parm.path <- paste0(parm.path.1, ".", parent.node[wlcl][[1]]$Name)
  }
  if (node == "Weather") {
    wlw <- function(x) grepl("Weather", x$Name)
    wlwl <- sapply(parent.node, FUN = wlw)
    weather.node <- parent.node[wlwl]
    gf1 <- function(x) grep(".met$", x, value = TRUE)
    cat("Met file:", as.character(sapply(weather.node,
                                         gf1)), "\n")
    parm.path <- paste0(parm.path.1, ".", parent.node[wlwl][[1]]$Name)
  }
  wcz <- grepl("Models.Core.Zone", parent.node)
  core.zone.node <- parent.node[wcz][[1]]$Children
  parm.path.2 <- paste0(parm.path.1, ".", parent.node[wcz][[1]]$Name)
  if (node == "Soil") {
    wsn <- grepl("Models.Soils.Soil", core.zone.node)
    soil.node <- core.zone.node[wsn]
    parm.path.2.1 <- paste0(parm.path.2, ".", soil.node[[1]]$Name)
    cat("Soil Type: ", soil.node[[1]]$SoilType, "\n")
    cat("Latitude: ", soil.node[[1]]$Latitude, "\n")
    cat("Longitude: ", soil.node[[1]]$Longitude, "\n")
    if (length(soil.node) != 1)
      stop("soil.node not equal to one")
    soil.children.names <- sapply(soil.node[[1]]$Children,
                                  function(x) x$Name)
    cat("Soil children:", soil.children.names, "\n")
    if (soil.child == "Metadata") {
      parm.path <- parm.path.2.1
      metadata <- NULL
      for (i in names(soil.node[[1]])) {
        if (i %in% c("Name", "Children",
                     "IncludeInDocumentation", "Enabled",
                     "ReadOnly"))
          next
        val <- as.character(ifelse(is.null(soil.node[[1]][[i]]),
                                   NA, soil.node[[1]][[i]]))
        if (!is.na(val) && nchar(val) > options()$width -
            30)
          val <- paste(strtrim(val, options()$width -
                                 30), "...")
        metadata <- rbind(metadata, data.frame(parm = i,
                                               value = val))
      }
      if (missing(parm)) {
        print(knitr::kable(metadata, longtable = FALSE))
      }
      else {
        if (!(parm %in% metadata[["parm"]]))
          stop("parm does not match a parameter in metadata")
        print(knitr::kable(metadata[metadata$parm ==
                                      parm, ]))
      }
    }
    else {
      wsc <- grep(soil.child, soil.children.names)
      if (length(wsc) == 0)
        stop("soil.child likely not present")
      selected.soil.node.child <- soil.node[[1]]$Children[wsc]
    }
    first.level.soil <- c("Water", "Physical",
                          "Chemical", "Analysis", "InitialWater",
                          "InitialN", "SoilWater", "Analysis",
                          "CERESSoilTemperature", "Organic")
    if (soil.child %in% first.level.soil) {
      parm.path <- paste0(parm.path.2.1, ".", selected.soil.node.child[[1]]$Name)
      enms <- c("IncludeInDocumentation", "Enabled",
                "ReadOnly", "Children", "Name")
      cnms <- setdiff(names(selected.soil.node.child[[1]]),
                      enms)
      soil.d1 <- NULL
      soil.d2 <- NULL
      col.nms <- NULL
      for (ii in cnms) {
        tmp <- selected.soil.node.child[[1]][ii][[1]]
        if (length(tmp) == 0)
          next
        if (length(tmp) == 1) {
          soil.d1 <- rbind(soil.d1, data.frame(parm = ii,
                                               value = as.character(tmp)))
        }
        if (length(tmp) > 1) {
          col.nms <- c(col.nms, ii)
          vals <- as.vector(unlist(tmp))
          soil.d2 <- cbind(soil.d2, vals)
        }
      }
      if (missing(parm)) {
        if (!is.null(soil.d1))
          print(knitr::kable(soil.d1, digits = digits))
        if (!is.null(soil.d2)) {
          soil.d2 <- as.data.frame(soil.d2)
          names(soil.d2) <- col.nms
          print(knitr::kable(soil.d2, digits = digits))
        }
      }
      else {
        if (!is.null(soil.d1))
          print(knitr::kable(soil.d1[soil.d1$parm ==
                                       parm, ], digits = digits))
        if (!is.null(soil.d2)) {
          soil.d2 <- as.data.frame(soil.d2)
          names(soil.d2) <- col.nms
          print(knitr::kable(soil.d2[soil.d2$parm ==
                                       parm, ], digits = digits))
        }
      }
    }
  }
  if (node == "SurfaceOrganicMatter") {
    wsomn <- grepl("Models.Surface.SurfaceOrganicMatter",
                   core.zone.node)
    som.node <- core.zone.node[wsomn][[1]]
    parm.path <- paste0(parm.path.2, ".", som.node$Name)
    som.d <- data.frame(parm = names(som.node)[2:8], value = as.vector(unlist(som.node)[2:8]))
    print(knitr::kable(som.d, digits = digits))
  }
  if (node == "MicroClimate") {
    wmcn <- grepl("Models.MicroClimate", core.zone.node)
    microclimate.node <- core.zone.node[wmcn][[1]]
    parm.path <- paste0(parm.path.2, ".", microclimate.node$Name)
    microclimate.d <- data.frame(parm = names(microclimate.node)[2:9],
                                 value = as.vector(unlist(microclimate.node)[2:9]))
    print(knitr::kable(microclimate.d, digits = digits))
  }
  if (node == "Crop") {
    wmmn <- grepl("Models.Manager", core.zone.node)
    manager.node <- core.zone.node[wmmn]
    wcn <- grepl("CultivarName", manager.node)
    crop.node <- manager.node[wcn][[1]]$Parameters
    parm.path <- paste0(parm.path.2, ".", manager.node[wcn][[1]]$Name)
    mat <- matrix(NA, nrow = length(crop.node), ncol = 2,
                  dimnames = list(NULL, c("parm", "value")))
    j <- 1
    for (i in 1:length(crop.node)) {
      mat[j, 1] <- crop.node[[i]]$Key
      mat[j, 2] <- crop.node[[i]]$Value
      j <- j + 1
    }
    print(knitr::kable(as.data.frame(mat), digits = digits))
  }
  if (node == 'Cultivar') {
    if (length(grep('Models.PMF.Cultivar, Models', capture.output(str(apsimx_json)), ignore.case = T)) > 0) {
      # There is at least one edited cultivar
      tmp <- search_apsimx(apsimx_json, keyword = 'Models.PMF.Cultivar, Models')
      cat('Name = ', tmp$Name)
      cat("\n")
      print(knitr::kable(data.frame(Command = unlist(tmp$Command))))
      cat("\n")
    } else {
      cat("No edited cultivars. Use function edit_apsimx() to add one.\n")
    }
    parm.path <- ""
    print.path <- F
  }
  if (node == "Manager") {
    wmmn <- grepl("Models.Manager", core.zone.node)
    manager.node <- core.zone.node[wmmn]
    parm.path <- parm.path.2
    manager.node.names <- sapply(manager.node, FUN = function(x) x$Name)
    cat("Management Scripts: ", manager.node.names,
        "\n\n")
    if (!is.null(parm)) {
      parm1 <- parm[[1]]
      position <- parm[[2]]
      find.manager <- grep(parm1, manager.node.names, ignore.case = TRUE)
      selected.manager.node <- manager.node.names[find.manager]
      parm.path <- paste0(parm.path.2, ".", selected.manager.node)
      if (is.na(position)) {
        ms.params <- manager.node[[find.manager]]$Parameters
        if (length(ms.params) == 0)
          warning("parameter not found")
        mat <- matrix(NA, ncol = 2, nrow = length(ms.params),
                      dimnames = list(NULL, c("parm", "value")))
        if (length(ms.params) > 0) {
          for (j in 1:length(ms.params)) {
            mat[j, 1] <- ms.params[[j]]$Key
            mat[j, 2] <- ms.params[[j]]$Value
          }
        }
        cat("Name: ", selected.manager.node, "\n")
        print(knitr::kable(as.data.frame(mat), digits = digits))
        cat("\n")
      }
      if (!is.na(position)) {
        ms.params <- manager.node[[find.manager]]$Parameters
        if (length(ms.params) == 0)
          warning("no parameters found")
        mat <- matrix(NA, ncol = 2, nrow = length(position),
                      dimnames = list(NULL, c("parm", "value")))
        k <- 1
        for (j in 1:length(ms.params)) {
          if (j == position) {
            mat[k, 1] <- ms.params[[j]]$Key
            mat[k, 2] <- ms.params[[j]]$Value
            k <- k + 1
          }
        }
        cat("Name: ", selected.manager.node, "\n")
        parm2 <- ms.params[[position]]$Key
        cat("Key:", ms.params[[position]]$Key,
            "\n")
        print(knitr::kable(as.data.frame(mat), digits = digits))
        cat("\n")
      }
    }
  }
  if (node == "Other") {
    tmp <- core.zone.node
    parm.path.2.1 <- parm.path.2
    if (is.null(parm))
      stop("'parm' should be provided when node = 'Other'")
    if (length(parm) == 1L)
      stop("'parm' should be a list of length 2 or more")
    for (i in 1:(length(parm) - 1)) {
      nms <- sapply(tmp, function(x) x$Name)
      wcp <- grep(parm[[i]], nms)
      if (length(wcp) == 0) {
        cat("Names: ", nms, "\n")
        cat("parm[[i]]", parm[[i]], "\n")
        stop("Parameter not found")
      }
      tmp <- tmp[[wcp]]
      if (!is.null(tmp$Children))
        tmp <- tmp$Children
      parm.path.2.1 <- paste0(parm.path.2.1, ".",
                              nms[wcp])
    }
    if (!is.null(tmp$Parameters)) {
      wp <- grep(parm[[length(parm)]], tmp$Parameters)
      tmp2 <- tmp$Parameters[[wp]]
      parm.path <- paste0(parm.path.2.1, ".", tmp2$Key)
      print(knitr::kable(as.data.frame(tmp2)))
    }
    else {
      parm.path <- parm.path.2.1
      unpack_node(tmp)
    }
  }
  if (node == "Report") {
    tmp <- search_apsimx(apsimx_json, keyword = 'VariableNames')
    print(knitr::kable(data.frame(VariableNames = unlist(tmp$VariableNames))))
    cat("\n")
    print(knitr::kable(data.frame(EventNames = unlist(tmp$EventNames))))
    parm.path <- ""
    print.path <- F
  }
  if (print.path && node != "Other") {
    if (!missing(parm)) {
      if (length(parm) == 1) {
        parm.path <- paste0(parm.path, ".", parm)
      }
      else {
        if (!is.na(position)) {
          parm.path <- paste0(parm.path, ".", parm2)
        }
      }
    }
    cat("Parm path:", parm.path, "\n")
  }
  else {
    if (print.path)
      cat("Parm path:", parm.path, "\n")
  }
  invisible(parm.path)
}

#=============================================================================

# Plot Yield predictions with conditional probabilities
# Function to plot

plot_conditional2 <- function(x, y, quantiles = c(.2, .5, .8), gp=NULL,
    point_size=0.9, vjust=-0.2){
  index <- !is.na(x) & !is.na(y)
  x <- x[index]
  y <- y[index]
  if(is.null(gp)){
    gp <- factor(rep(1,length(y)))
  }else{
    gp <- gp[index]
  }
  xq <- quantile(x, probs = quantiles)
  yq <- quantile(y, probs = quantiles)
  rho <- sprintf('%.3f', as.numeric(cor(x, y)))

  dat <- data.frame(x,y,gp)
  dat2 <- data.frame(qq=names(yq), xq, yq, ymin=-Inf, xmin=-Inf)

  ggplot(dat,aes(x,y)) + theme_bw() +
    geom_point(aes(fill=gp),shape=21,color="gray20",size=point_size) +
    geom_vline(xintercept=xq, color=4, linetype=3, size=0.5) +
    geom_label(data=dat2,aes(x=xq,y=ymin, label=qq), color=4, vjust=vjust, label.size = NA) +
    geom_hline(yintercept=yq, color=4, linetype=3, size=0.5) +
    geom_label(data=dat2,aes(x=xmin,y=yq, label=qq), color=4, hjust=vjust, label.size = NA)
}

#=============================================================================

plot_conditional <- function(x, y, quantiles = c(.2, .5, .8), xpos = .05, ypos = .05, cex.text=0.7, cex.point=0.6, ...) {
  index <- !is.na(x) & !is.na(y)
  x <- x[index]
  y <- y[index]
  xq <- quantile(x, probs = quantiles)
  yq <- quantile(y, probs = quantiles)
  rho <- sprintf('%.3f', as.numeric(cor(x, y)))
  plot(x, y, col = 'grey', cex=cex.point, ...)
  abline(0,1, col = 2)
  text(min(x), max(y) - 0.05*diff(range(y)), labels=bquote(rho * ' = ' * .(rho)), col = "green4", pos = 4)

  for(i in 1:length(xq)) {
    abline(h = yq[i], col = 4, lty = 3)
    abline(v = xq[i], col = 4, lty = 3)
    #
    chr_x_space <- diff(par("usr")[1:2]) * xpos
    chr_y_space <- diff(par("usr")[3:4]) * ypos
    tmp_x_pos <- par("usr")[1] + nchar(names(yq)[i]) * chr_x_space
    tmp_y_pos <- par("usr")[3] + chr_y_space
    #
    #rect(xleft = par("usr")[1] + chr_x_space * .3, xright = tmp_x_pos,
    #     ybottom = yq[i] - chr_y_space, yq[i] + chr_y_space, col = 'white', border = NA)
    text(tmp_x_pos, yq[i], labels = names(yq)[i], col = 4, pos = 2)
    #
    #rect(xleft = xq[i] - chr_x_space, xright = xq[i] + chr_x_space,
    #     ybottom = tmp_y_pos - chr_y_space, tmp_y_pos + chr_y_space, col = 'white', border = NA)
    text(xq[i], tmp_y_pos, labels = names(xq)[i], col = 4)
  }
  x_cut <- cut(x, c(min(x)-1, xq, max(x)))
  y_cut <- cut(y, c(min(y)-1, yq, max(y)))
  for (i in 1:length(levels(x_cut))) {
    for (j in 1:length(levels(y_cut))) {
      tmp <- x_cut == levels(x_cut)[i] & y_cut == levels(y_cut)[j]
      tmp_x_lim <- as.numeric(gsub('\\[|\\]|[()]','', unlist(strsplit(levels(x_cut)[i], ','))))
      tmp_y_lim <- as.numeric(gsub('\\[|\\]|[()]','', unlist(strsplit(levels(y_cut)[j], ','))))
      text(mean(tmp_x_lim), mean(tmp_y_lim), round(sum(tmp) / table(x_cut)[i], 3),cex=cex.text)
    }
  }
}

#=============================================================================

make_heatmap <- function(ec, W, layer=NULL, cluster_rows=TRUE, cluster_cols=TRUE, ...)
{
  phase <- unlist(lapply(strsplit(colnames(W),"_"),function(x)x[1])) # Phases
  Phases <- unique(phase)

  ec0 <- unlist(lapply(strsplit(colnames(W),"_"),function(x)x[2]))

  ec <- grep(paste(ec,collapse="|"),ec0,value=T)
  drop <- grep("FlowNO3|PAWmm",ec)
  if(length(drop)>0) ec <- ec[-drop]

  if(!is.null(layer)){
    tmp <- (1:10)[!(1:10) %in% layer]
    drop <- grep(paste(paste0("(",tmp,")"),collapse="|"),ec)
    if(length(drop)>0) ec <- ec[-drop]
  }

  index <- which(ec0 %in% ec)

  namesPhases <- paste0("P",1:length(Phases)); names(namesPhases) <- Phases
  ec2 <- paste0(namesPhases[phase][index],"-",ec0[index])

  W0 <- W[,index]
  colnames(W0) <- ec2

  annot <- data.frame(Phase=factor(namesPhases[phase][index],levels=namesPhases))
  rownames(annot) <- colnames(W0)

  pheatmap(cor(W0),cluster_rows = cluster_rows, cluster_cols = cluster_cols,
            show_rownames = T, show_colnames = T,
            annotation_col=annot, annotation_row=annot, ...)
}


#=============================================================================

get_frontier <- function(dat,x="ASI_GDDC",y="yield",
              probs=c(0.6,0.7,0.8,0.9), x.breaks_n=10,
              x.breaks_type=c("percentile","evenly_spaced"))
{
  x.breaks_type <- match.arg(x.breaks_type)

  dat <- dat[!is.na(dat[,x]) & !is.na(dat[,y]),]

  br <- seq(min(dat[,x]),max(dat[,x]),length=x.breaks_n+1)
  if(x.breaks_type=="percentile"){
     br <- quantile(dat[,x], probs=seq(0,1,length=x.breaks_n+1))
  }
  br2 <- cbind(br[1:x.breaks_n],br[2:(x.breaks_n+1)])

  out <- c()
  for(j in 1:x.breaks_n){
      a1 <- ifelse(j==1,br2[j,1]-1,br2[j,1])
      a2 <- br2[j,2]
      index <- which(dat[,x] > a1 & dat[,x] <= a2)
      #cat(range(dat[index,x]),"\n")
      q <- quantile(dat[index,y], probs)
      out <- rbind(out,data.frame(xmean=mean(br2[j,]),p=names(q),q))
  }
  pp <- ggplot(dat, aes_string(x,y)) +
          geom_point(color="gray70",size=0.4,shape=21) + theme_bw() +
          geom_line(data=out,aes(xmean,q,color=p,group=p),size=0.5) +
          geom_point(data=out,aes(xmean,q,color=p))
  pp
}

#=============================================================================

get_outliers <- function(y,kIQR=2){
  qq <- quantile(y,probs=c(0.25,0.75))
  IQR <- as.vector(qq[2] - qq[1])
  which((y < qq[1] - kIQR*IQR) | (y > qq[2] + kIQR*IQR))
}
#=============================================================================

get_variance_SDR <- function(sdr, rgTT=c(-100,100), a=-5,
                        qb=0.5,transform=TRUE, br)
{

  TT <- lapply(strsplit(colnames(sdr),"\\(|\\)-"),function(x)x[2])
  TT <- 100*as.numeric(unlist(TT))

  index <- which(TT >= rgTT[1] & TT <= rgTT[2])
  X <- sdr[, index, drop=FALSE]

  x0 <- apply(X,1,mean)
  if(transform){
    b=-1*as.vector(quantile(x0, probs=qb))
    eta <- a*(x0-b)
    x1 <- exp(eta)/(1+exp(eta))
  }else{
     x1 <- x0
     x1 <- x1 - min(x1)
     x1 <- x1/max(x1)
  }

  K <- length(br)
  gp <- rep(NA,K)

  rg <- range(unlist(br))
  for(k in 1:K){
    a1 <- br[[k]][1]
    a2 <- br[[k]][2]
    if(a1==rg[1]){
      index1 <- x1 >= a1
    }else index1 <- x1 > a1

    index2 <- x1 <= a2
    gp[which(index1 & index2)] <- names(br)[k]
  }
  if(any(is.na(gp))) stop("Some entries were not classified with the provided breaks")

  gp <- factor(gp, levels=names(br))
  names(gp) <- rownames(SDR)

  tmp <- data.frame(y=x1, Groups=gp)
  fm <- anova(lm(y~Groups, data=tmp))
  F0 <- fm["F value"][1,1]

  SS <- fm["Sum Sq"]
  DF <- fm["Df"]
  Pr <- fm["Pr(>F)"][1,1]

  list(Group=gp, DF=DF, SS=SS, F=F0, Pr=Pr)

}

get_variance_SDR2 <- function(sdr, rgTT=c(-100,100), K=2,
                        a=5, qb=0.5, transform=TRUE,
                        fun.stat=c("mean","prod"))
{
  fun.stat <- match.arg(fun.stat)
  TT <- lapply(strsplit(colnames(sdr),"\\(|\\)-"),function(x)x[2])
  TT <- 100*as.numeric(unlist(TT))

  index <- which(TT >= rgTT[1] & TT <= rgTT[2])
  X <- sdr[, index, drop=FALSE]

  if(fun.stat=="mean"){
    x0 <- apply(X,1,mean)
  }else{
    x0 <- apply(X,1,prod)
  }

  if(transform){
    b <- as.vector(quantile(x0, probs=qb))
    eta <- a*(x0-b)
    x1 <- exp(eta)/(1+exp(eta))
  }else{
     x1 <- x0
  }

  br <- quantile(x1,probs=seq(0,1,by=1/K))

  br2 <- cbind(br[1:K],br[-1])

  gp <- rep(NA,K)
  rg <- range(x1)

  for(k in 1:K){
    a1 <- br2[k,1]
    a2 <- br2[k,2]
    if(k==1){
      index1 <- x1 >= a1
    }else index1 <- x1 > a1

    index2 <- x1 <= a2
    gp[which(index1 & index2)] <- k
  }
  if(any(is.na(gp))) stop("Some entries were not classified with the provided breaks")

  names(gp) <- rownames(sdr)

  tmp <- data.frame(y=x1, Groups=factor(gp))
  fm <- anova(lm(y~Groups, data=tmp))
  F0 <- fm["F value"][1,1]

  SS <- fm["Sum Sq"]
  DF <- fm["Df"]
  Pr <- fm["Pr(>F)"][1,1]

  list(Group=gp, DF=DF, SS=SS, F=F0, Pr=Pr)

}

#=============================================================================
get_cluster_EC <- function(W, rgTT=c(-100,100),
                    method=c("hclust","kmeans"),k=4, decreasing=TRUE)
{
  method <- match.arg(method)
  TT <- lapply(strsplit(colnames(W),"\\(|\\)-"),function(x)x[2])
  TT <- 100*as.numeric(unlist(TT))

  index <- which(TT >= rgTT[1] & TT <= rgTT[2])
  W <- W[, index, drop=FALSE]

  # Grouping
  if(method=="hclust"){
    cl <- hclust(dist(W), method="ward.D")
    Group <- cutree(cl, k=k)
  }

  if(method=="kmeans"){
    cl <- kmeans(W, centers=k, nstart=50, iter.max=1000)
    Group <- cl$cluster
  }
  # Make boxplot by group and rearrenge acording to mean values
  tmp <- data.frame(Group,X=I(W))
  tmp <- split(tmp, tmp$Group)
  bp <- do.call(rbind,lapply(tmp, function(x){
    a1 <- apply(x$X,1,mean)
    data.frame(Group=x$Group[1],Var=names(a1),value=a1)
  }))
  rownames(bp) <- NULL
  tmp <- aggregate(value~Group, data=bp, median)
  tmp <- tmp[order(tmp$value, decreasing=decreasing),]

  # Change order of levels
  ll <- as.numeric(factor(as.character(Group),levels=tmp$Group))
  names(ll) <- names(Group)
  Group <- ll

  SVD <- svd(W)
  d <- SVD$d
  rk <- qr(W)$rank
  v <- d^2/sum(d^2)
  SDVe <- (-1/log(rk))*sum((v*log(v))[1:rk])
  PCvar <- round(100*d/sum(d),1)
  PC <- SVD$u
  for(j in 1:ncol(PC)) PC[,j] <- PC[,j]*d[j]
  rownames(PC) <- rownames(W)

  tmp <- ifelse(ncol(PC)>=2,2,ncol(PC))
  INFO <- data.frame(rownames(W),do.call(rbind,strsplit(rownames(W),"-")))
  dat <- data.frame(INFO,Group,SVD$u[,1:tmp, drop=F])
  colnames(dat) <- c("year_loc","year","location","Group",paste0("PC",1:tmp))
  dat$year_loc2 <- paste0(dat$location,"-",substr(dat$year,3,4))

  if(ncol(PC) > 1){
    tmp <- range(TT[index])
    title0 <- paste0("TT in [",tmp[1]-50,", ",tmp[2]+50,"]")
    pp <- ggplot(dat, aes(PC1,PC2,color=factor(Group))) +
      geom_label(aes(label=year_loc2),size=1.6) + theme_bw() +
      labs(x=paste0("PC1 (",PCvar[1],"%)"), color="Group",
           y=paste0("PC2 (",PCvar[2],"%)"), title=title0) +
      theme(legend.position="none", plot.title=element_text(hjust=0.5)) +
      guides(color = guide_legend(override.aes=list(size = 3))) +
      annotate("text",x=min(dat$PC1),y=max(dat$PC2),
        label=paste0("Entropy=",round(SDVe,4)),color="blue",hjust=-0.01)
  }else pp <- NULL

  list(pp=pp, Group=Group, PC=PC, W=W, Entropy=SDVe)
}

#=============================================================================
#dat=OUT; x="yHat"; y="y"; ylab="Observed"; xlab="Predicted"; showcor=TRUE;
#color="cluster"; title=vline=hline=label=NULL
my_scatterplot <- function(dat, x, y, xlab="x", ylab="y", title=NULL, showcor=TRUE,
            vline=NULL, hline=NULL, label=NULL, color=NULL, facet=NULL, scales=NULL)
{
  # xlab="x"; ylab="y"; title=NULL; showcor=TRUE; vline=NULL; hline=NULL; label=NULL; color=NULL
  drop <- which(apply(dat[,c(x,y)],1,function(x)any(is.na(x))))
  if(length(drop) > 0) dat <- dat[-drop,]

  rg <- range(dat[,x],dat[,y])
  rgx <- range(dat[,x])
  rgy <- range(dat[,y])

  if(is.null(color)){
    dat$color <- factor(1)
  }else dat$color <- dat[,color]

  if(is.null(facet)){
    dat$facet <- factor(" ")
  }else{
    tmp <- dat[,facet]; tmp <- ifelse(is.na(tmp),"none",tmp)
    dat$facet <- tmp
  }

  if(!is.null(label)){
    dat$label <- factor(as.character(dat[,label]))
  }
  #corre <- cor(dat[,x],dat[,y])
  #corre <- sprintf('%.2f', corre)
  if(is.null(scales)) scales <- "fixed"
  dat0 <- do.call(rbind,lapply(split(dat,dat$facet), function(tt){
    rr = lm(formula=paste(y,"~",x), data=tt)
    data.frame(facet=tt$facet[1], R2=summary(rr)[["r.squared"]],
       cor=cor(tt[,x],tt[,y]),
       minx=min(tt[,x]),maxx=max(tt[,x]),
       miny=min(tt[,y]),maxy=max(tt[,y]))
  }))
  dat0$label1 <- paste('R^2 ==',sprintf('%.3f',dat0$R2))
  dat0$label2 <- paste('rho ==',sprintf('%.3f',dat0$cor))

  if(all(!c("free","free_x")%in%scales)){
    dat0[,"minx"]=min(dat[,x]); dat0[,"maxx"]=max(dat[,x])
  }
  if(all(!c("free","free_y")%in%scales)){
    dat0[,"miny"]=min(dat[,y]); dat0[,"maxy"]=max(dat[,y])
  }

  pp <- ggplot(dat,aes_string(x,y)) +
     labs(x=xlab,y=ylab,fill=NULL,title=title) + theme_bw() + # lims(x=rg, y=rg) +
     geom_smooth(method='lm', formula= y~x, size=0.5,se=FALSE) +
     #stat_regline_equation(label.y=rgy[2],aes(label= ..rr.label..),size=3.8) +
     geom_text(data=dat0,aes(x=minx,y=maxy,
                   label=label1), parse=TRUE, hjust=0, vjust=0.7) +
     theme(legend.position=c(0.99,0.01),
           legend.justification=c(1,0),
           legend.key.height=unit(0.75,"line"),
           legend.margin=margin(t=-0.17,b=0.17,l=0.17,r=0.17,unit='line'),
           legend.box.background = element_rect(colour = "gray32",size=0.8),
           legend.background = element_rect(fill="gray92"),
           plot.title=element_text(hjust=0.5),
           strip.text.x = element_text(size = 10,margin = margin(t=2,b=2))
           )
  if(!is.null(facet)){
     pp <- pp + facet_wrap(~facet, scales=scales)
  }
  #else{

  #}

  if(showcor){
    pp <- pp + geom_text(data=dat0,aes(x=minx,y=maxy-0.09*(maxy-miny),
                  label=label2), parse=TRUE, hjust=0, vjust=0.7)
  }

  if(!is.null(vline)){
     pp <- pp + geom_vline(xintercept=vline, linetype="dashed",color="red",size=0.2)
  }
  if(!is.null(hline)){
     pp <- pp + geom_hline(yintercept=hline, linetype="dashed",color="red",size=0.2)
  }

  if(is.null(label)){
    pp <- pp + geom_point(aes(fill=color),shape=21,size=1.9)
    if(length(unique(dat$color)) == 1){
        pp <- pp + theme(legend.position="none")
    }
  }else{
    pp <- pp + geom_label(aes(label=label,fill=color),size=2.1) +
              theme(legend.position=c(0.99,0.01))
  }
  pp
}

#=============================================================================
# dat=dat0; x="Var2"; y="value"; ind="Var1"; group="Group"; xlab=NULL; highlight="MIH";
# rgTT=c(-300,300); ylab="y"; linewd=0.05; grand_mean=group_mean=FALSE; x.by=100; color.line="lightblue"
# legend.pos=c("topleft"); color=GroupColor; text.size=2; text.xpos=100; title=NULL; factorwd=8

my_xyplot_TT <- function(dat, x="Var2", y="value", ind="Var1",
                    group="Group", xlab=NULL, ylab="y", highlight=NULL, text.size=2,
                    text.vjust=0.5,text.color="blue",text.xpos=100,
                    linewd=0.2, factorwd=8, group_mean=TRUE, grand_mean=TRUE, color=NULL,
                    legend.pos=c("topleft","topright","bottomleft","bottomright","none"),
                    title=NULL, rgTT=NULL, x.by=200, color.line="lightblue")
{
  legend.pos <- match.arg(legend.pos)

  if(is.null(xlab)) xlab=expression(paste("TT around flowering (DD)"))
  dat$tt <- 100*as.numeric(unlist(lapply(strsplit(as.character(dat[,x]),"\\(|\\)-"),function(x)x[2])))

  dat <- do.call(rbind,lapply(split(dat,paste(dat$tt,dat[,ind])),function(z){
    rr <- z[1,]
    rr[,y] <- sum(z[,y])
    rr
  }))
  rownames(dat) <- NULL

  if(!is.null(rgTT)){
    stopifnot(length(rgTT)==2)
    stopifnot(rgTT[2] > rgTT[1])
    index <- which(dat$tt >= rgTT[1] & dat$tt <= rgTT[2])
    if(length(index)>0) dat <- dat[index,]
  }

  mm <- do.call(rbind,lapply(split(dat,paste(dat$tt)),function(z){
    data.frame(z[1,c(x,"tt")], t(apply(z[,y, drop=F],2,mean)))
  }))

  if(!is.null(group)){
    mm_gp <- do.call(rbind,lapply(split(dat,paste(dat$tt,dat[,group])),function(z){
      data.frame(z[1,c(x,"tt",group)], t(apply(z[,y, drop=F],2,mean)))
    }))
  }else mm_gp <- NULL

  if(is.null(color)){
    color <- gg_color_hue(nlevels(dat[,group]))
    names(color) <- levels(dat[,group])
  }else{
     stopifnot(all(levels(dat[,group]) %in% names(color)))
  }
  rg <- range(dat$tt)
  rg[1] <- rg[1]-50
  rg[2] <- rg[2]+50
  br <- seq(rg[1],rg[2], by=100)
  labels0 <- rep("", length(br))
  tmp <- seq(rg[1],rg[2], by=x.by)
  labels0[match(tmp,br)] <- tmp

  #legend.pos=c("topleft","topright","bottomleft","bottomright","none")[4]
  tmp <- switch(legend.pos,"topleft"=list(c(0,1),c(0.01,0.99)),
                     "topright"=list(c(1,1),c(0.99,0.99)),
                     "bottomleft"=list(c(0,0),c(0.01,0.01)),
                     "bottomright"=list(c(1,0),c(0.99,0.01)),
                     "none"=list(c(0,0),"none"))
  lgdjust <- tmp[[1]]
  lgdpos <- tmp[[2]]

  dat2 <- dat_label <- NULL
  if(!is.null(highlight)){
    tmp <- grep(highlight, dat[,ind])
    if(length(tmp)==0) stop("Entry '",highlight,"' was not found in column '",ind,"'")
    dat2 <- dat[tmp,]
    dat2$label <- unlist(lapply(strsplit(as.character(dat2[,ind]),"-"),function(x)x[1]))

    dat_label <- do.call(rbind,lapply(split(dat2, as.character(dat2[,ind])), function(x){
      tmp <- abs(x$tt - text.xpos)
      data.frame(x[1,c(ind,group,"label")],ypos=mean(x[abs(tmp-min(tmp)) < 1E-5, y]))
    }))
    dat_label <- dat_label[order(dat_label$ypos), ]
    dat_label$vjust <- text.vjust*(-1)^(1+(1:nrow(dat_label)))

  }

  pp <- ggplot(dat, aes_string("tt",y)) + theme_bw() +
        geom_vline(xintercept=0, color="gray75", linetype="dashed", size=0.4)

  if(is.null(group)){
    pp <- pp + geom_line(aes_string(group=ind),color=color.line,size=linewd)
  }else{
    pp <- pp + geom_line(aes_string(group=ind, color=group),size=linewd) +
               theme(legend.position=lgdpos,
                     legend.justification=lgdjust,
                     legend.key.height=unit(0.7,"line"),
                     legend.margin=margin(t=-0.15,b=0.15,l=0.15,unit='line'),
                     legend.background = element_rect(fill="gray92")) +
               guides(color = guide_legend(override.aes = list(size = 0.6))) +
               scale_color_manual(values=color)
    if(group_mean){
      pp <- pp + geom_line(data=mm_gp,aes_string(group=group, color=group),size=linewd*factorwd)
    }
  }
  if(!is.null(highlight)){
    pp <- pp + geom_line(data=dat2, aes_string(group=ind, color=group),size=linewd*factorwd) +
               geom_text(data=dat_label, aes(y=ypos, label=label, vjust=vjust),x=text.xpos,
                         color=text.color,size=text.size)
  }
  if(grand_mean){
    pp <- pp + geom_line(data=mm,color="gray40",size=linewd*factorwd)
  }

  pp <- pp + scale_x_continuous(breaks=br, labels=labels0,expand=c(0,0),
                  limits=c(rg[1],rg[2])) +
             labs(x=xlab, y=ylab, title=title) +
             theme(plot.title=element_text(hjust=0.5),
                  panel.grid.major.x=element_blank(),
                  panel.grid.minor.y=element_blank())


  pp

}

#=============================================================================

my_levelplot_TT <- function(COV,xlab="Phase j",ylab="Phase i",
                        title=NULL, show.legend=FALSE, show.value=FALSE,
                        legend.name="",text.col="gray30", digits=2)
{
  dat <- melt(t(COV))
  dat<- dat[!is.na(dat$value),]
  dat$Var1 <- factor(as.character(dat$Var1),levels=colnames(COV))
  dat$Var2 <- factor(as.character(dat$Var2),levels=rownames(COV))

  dat$x <- as.numeric(dat$Var1)
  dat$y <- max(as.numeric(dat$Var2)) - as.numeric(dat$Var2) +1
  ee <- 0.01
  pp <- ggplot(dat,aes(x,y)) + theme_bw() + geom_tile(aes(fill=value)) +
    scale_fill_gradientn(colours = rev(viridis(10))) +
    scale_y_continuous(breaks=1:nrow(COV),labels=rev(rownames(COV)), expand=c(ee,ee)) +
    scale_x_continuous(sec.axis=sec_axis(~.+0,xlab,breaks=1:ncol(COV),labels=colnames(COV)),
                       expand=c(ee,ee)) +
    labs(x=NULL, y=ylab, title=title, fill=legend.name) +
    theme(axis.ticks.x.bottom = element_blank(),
          axis.text.x.bottom = element_blank(),
          panel.grid = element_blank(),
          plot.title=element_text(hjust=0.5))
   if(show.value){
      pp <- pp + geom_text(aes(label=sprintf(paste0("%.",digits,"f"),value)),
                  color=text.col,size=3)
   }
   if(!show.legend) pp <- pp + theme(legend.position="none")

   pp
}

#=============================================================================

get_mean <- function(dat, x="Var2", y="value", group="Var1", stat=c("mean","sum"))
{
  stat <- match.arg(stat)
  dat$tt <- 100*as.numeric(unlist(lapply(strsplit(as.character(dat[,x]),"\\(|\\)-"),function(x)x[2])))
  tmp <- dat$tt
  if(!is.null(group)){
    tmp <- paste0(tmp,"-",apply(dat[,group, drop=F],1,paste,collapse="-"))
  }
  out <- do.call(rbind,lapply(split(dat,tmp),function(z){
    rr <- z[1,]
    if(stat=="mean"){
      aa <- apply(z[,y,drop=F],2,mean)
    }else aa <- apply(z[,y,drop=F],2,sum)
    rr[,y] <- aa
    rr
  }))
  rownames(out) <- NULL
  out <- out[,colnames(out)!="tt"]
  out
}

#=============================================================================

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#=============================================================================
# get_distance (in meters) between two points P1 and P2 using latitude and longitude
# P1 = c(lat1, long1)    &     P2 = c(lat2, long2)
# P1 = c(40.48800, -87.00600);  P2 = c(40.47927, -86.98905)
# P1=c(40.6892, -74.0442); P2=c(38.8890, -77.0345)
get_distance <- function(P1 , P2){
  P1 <- as.vector(as.matrix(P1))
  P2 <- as.vector(as.matrix(P2))
  R <- 6371E3 # R is earth’s radius
  x1 <- P1[1] * pi /180
  x2 <- P2[1] * pi /180
  dx <- (P2[1]-P1[1]) * pi /180
  dy <- (P2[2]-P1[2]) * pi /180
  a <- (sin(dx/2)^2) + cos(x1)*cos(x2)*(sin(dy/2)^2)
  b <- 2*atan2(sqrt(a), sqrt(1-a))
  R*b
}

#=============================================================================
# group=latlon=latlon2=d.max=NULL; metric="median"; verbose=TRUE; id.col=NULL
distance_to_mean <- function(data, group=NULL, latlon=NULL,
                             latlon2=NULL, d.max=NULL,
                             metric=c("median","mean"),
                             id.col=NULL, verbose=TRUE){

  metric <- match.arg(metric)

  latlon0 <- c("latitude","longitude")

  if(is.null(group)){
    groupvar <- rep("none",nrow(data))
  }else{
    groupvar <- data[,group]
  }
  group_levels <- names(table(groupvar))

  if(is.null(latlon)){
    latlon <- latlon0
  }
  stopifnot(length(latlon)==2)

  if(is.null(id.col)){
    ID <- paste("Row",1:nrow(data))
  }else{
    ID <- data[,id.col]
  }

  if(verbose){
    if(is.null(group)){
      cat("Calculating the ",metric," distance across all data ...\n",sep="")
    }else{
     cat("Calculating the ",metric," distance within each of the ",
         length(group_levels)," level(s) of '",group,"' ...\n",sep="")
    }
  }

  OUT <- matrix(NA, ncol=4, nrow=nrow(data))
  colnames(OUT) <- c(latlon,"d","dMax")

  for(k in 1:length(group_levels)){
    index <- which(groupvar == group_levels[k])
    dat.tmp <- data[index,]
    tt <-  as.matrix(dat.tmp[,latlon])
    md <- as.vector(apply(tt, 2, metric, na.rm=T))

    d0 <- as.vector(apply(tt, 1, function(z)get_distance(z, md))/1000)
    out <- data.frame(tt,d=d0)

    if(!is.null(d.max) & !is.null(latlon2)){
      flag <- any(is.na(d0)) | (length(which(d0 > d.max)) > 0)
      if(flag){
        tt0 <- as.matrix(dat.tmp[,latlon2])
        d02 <- as.vector(apply(tt0, 1, function(z)get_distance(z, md))/1000)
        for(j in 1:nrow(out)){
           a1 <- ifelse(is.na(d0[j]),TRUE,d0[j]>d.max) # WS coordinates is NA or > dmid
           a2 <- ifelse(is.na(d02[j]),FALSE,d02[j]<d.max) # Field coordinates is not NA and < dmid

            if(a1 & a2){
              prefix1 <- paste(paste0(substr(latlon,1,3),"=",out[j,latlon]),collapse=",")
              prefix2 <- paste(paste0(substr(latlon,1,3),"=",tt0[j,]),collapse=",")
              out[j,] <- as.vector(c(tt0[j,],d=d02[j]))
              if(verbose){
                tmp <- ifelse(is.null(group),"",paste0(", ",group,"=",groupvar[index[j]]))
                msg <- paste(" - ",ID[index[j]],tmp,": ",
                    "\n\t","WS coordinates (",prefix1,") were changed to",
                    "\n\t","field coordinates (",prefix2,")\n",sep="")
                cat(msg,"\n")
              }
            }
        }
      }
    }
    out$dMax <- max(out$d, na.rm=TRUE)

    OUT[index,] <- as.matrix(out)
  }

  colnames(OUT)[1:2] <- latlon0
  data.frame(OUT)
}

#=============================================================================
# entry=NULL; latlon=NULL; d.max=NA; all=FALSE; metric="median"
get_coordinates_mean <- function(data, entry=NULL, latlon=NULL,
                                 metric=c("median","mean"),
                                 d.max=NA, all=FALSE)
{
  metric <- match.arg(metric)

  if(is.null(latlon)){
    latlon <- c("latitude","longitude")
  }
  stopifnot(length(latlon)==2)

  stopifnot(is.list(entry))
  if(any(!unlist(lapply(entry, function(x)length(x)==1)))){
    stop("Object 'entry' must be a list of entries with one element ")
  }
  if(any(!names(entry) %in% colnames(data))){
    stop("All 'names(entry)' must be contained in 'colnames(data)'")
  }

  np <- length(entry)
  tt <- vector('list', np*(np+1)/2)
  names0 <- c()
  cont <- 0
  for(i in 1:np){
    for(j in (i):np){
      cont <- cont + 1

      if(i==j){
        names0[cont] <- names(entry)[i]
        index <- which(data[,names(entry)[i]]==entry[[i]])
      }else{
        names0[cont] <- paste0(names(entry)[i],"_",names(entry)[j])
        index <- which(data[,names(entry)[i]]==entry[[i]] & data[,names(entry)[j]]==entry[[j]])
      }
      tt[[cont]] <- data[index, , drop=FALSE]
    }
  }
  names(tt) <- names0

  res <- do.call(rbind,lapply(1:length(tt), function(i){
    tt0 <- tt[[i]]

    if(!is.na(d.max)){
      tmp <- distance_to_mean(tt0, latlon=latlon, metric="median", verbose=FALSE)
      index <- which(tmp$d < d.max)
      n.dropped <- nrow(tt0) - length(index)
      tt0 <- tt0[index, , drop=FALSE]
    }else{
      n.dropped <- 0
    }

    item0 <- lapply(1:length(entry),function(k){
      x <- tt0[,names(entry)[k]]
      unique(x[!is.na(x)])
    })
    names(item0) <- names(entry)

    list0 <- lapply(item0,function(x)ifelse(length(x)==1,x,NA))
    n0 <- lapply(item0,length)
    names(n0) <- paste0("n.",names(n0))

    n.OK <- sum(apply(tt0[,latlon],1,function(x)all(!is.na(x))))

    data.frame(t(apply(tt0[, latlon],2,metric, na.rm=TRUE)),
               n0, n=nrow(tt0), n.dropped, n.OK, list0)
  }))
  rownames(res) <- names(tt)

  if(!all){
    index <- apply(res[,latlon],1,function(x)all(!is.na(x)))
    res <- res[index, , drop=F]
  }

  res
}

#=============================================================================

get_pairwise_distance <- function(data, latlon=NULL, group, metric=c("median","mean"))
{
  metric <- match.arg(metric)

  if(is.null(latlon)){
    latlon <- c("latitude","longitude")
  }
  stopifnot(length(latlon)==2)

  if(!group %in% colnames(data) & length(grep(":",group))>0){
    tmp <- unlist(strsplit(group,":"))
    stopifnot(all(tmp %in% colnames(data)))
    tmp <- data.frame(apply(data[,tmp],1,function(x)paste(as.character(x), collapse=":")))
    names(tmp) <- data
    data <- data.frame(data,tmp, check.names = FALSE)
  }
  dat2 <- do.call(rbind,lapply(split(data, data[,group]), function(x){
    if(metric=="median"){
      md <- apply(x[,latlon], 2, median, na.rm=T)
    }
    if(metric=="mean"){
      md <- apply(x[,latlon], 2, mean, na.rm=T)
    }
    c(md,n=nrow(x))
  }))

  index <- apply(dat2[,latlon],1, function(x)any(is.na(x)))
  dat2 <- dat2[!index,]

  res <- c()
  for(i in 1:(nrow(dat2)-1)){
    for(j in (i+1):nrow(dat2)){
      d <- get_distance(as.vector(dat2[i,latlon]), as.vector(dat2[j,latlon]))/1000
      res <- rbind(res,data.frame(rownames(dat2)[i],rownames(dat2)[j],d))
    }
  }
  colnames(res) <- c(paste0(group,1:2),"distance")
  res
}

#=============================================================================

perform_action <- function(data, actions)
{
  stopifnot(is.list(actions))
  index <- unlist(lapply(actions, function(x) all(c("COLUMN","ENTRY") %in% names(x))))
  actions <- actions[index]

  if(length(actions) == 0){
    stop("No actions to be performed")
  }

  out <- data[]
  for(j in 1:length(actions))
  {
    action <- actions[[j]]
    index <- which(data[,action$COLUMN] %in% action$ENTRY)

    if(length(index) > 0)
    {
      vv <- names(action)[names(action) %in% colnames(data)]
      if(length(vv)>0){
        for(k in 1:length(vv)){
          out[index, vv[k]] <- action[[vv[k]]]
        }
        a1 <- paste0(" - '",action$COLUMN,"=",paste(action$ENTRY,collapse=","),"'")
        a2 <- paste0("(n=",length(index),")")
        a3 <- paste(paste0(names(action[vv]),"=",action[vv]), collapse="|")
        msg <- paste0(a1," ", a2,":",
              ifelse(!is.null(action$MSG),paste0("\n\t",action$MSG),""),
              "\n\t",a3,"\n")
        cat(msg,"\n")
      }else{
        cat(" - No action was performed for 'actions[[",j,"]]'\n", sep="")
      }
    }else{
      cat("Not all item(s) '",paste(action$ENTRY,collapse=","),
          "' were not found in column '",action$COLUMN,"'\n",sep="")
    }
  }
  return(out)
}

#=============================================================================
# color="year_loc"; point.size=point.shape=1; show.mean=FALSE; show.all=TRUE; latlon=title=NULL
plot_US <- function(data, color="year_loc", pathToShp,
          point.size=1, point.shape=1, show.mean=FALSE,
          show.lines=FALSE,
          show.all=TRUE, latlon=NULL, title=NULL){

  if(is.null(latlon)){
    latlon <- c("latitude","longitude")
  }
  stopifnot(length(latlon)==2)

  # ogrInfo(dsn = pathToShp, layer="states")
  dpt <- readOGR(dsn = pathToShp, layer="states",
                 stringsAsFactors=FALSE, verbose=FALSE)

  lev_color <- unique(data[,color])
  color0 <- gg_color_hue(length(lev_color))
  names(color0) <- lev_color

  mp <- apply(data[,latlon], 2, median)

  if(show.lines  & !show.mean){
    cat("'show.mean' is set to TRUE when 'show.lines=TRUE'\n")
  }

  if(show.all){
    plot(dpt[!dpt$STATE_NAME %in% c("Hawaii","Alaska"),], col="gray95", main=title)
  }else{
    plot(dpt[dpt$STATE_ABBR %in% data$state,], col="gray95", main=title)
  }
  points(data[,latlon[2]], data[,latlon[1]], col=color0[data[,color]],
      pch=point.shape, cex=point.size)
  if(show.lines){
    for(i in 1:nrow(data)){
      lines(c(mp[2],data[i,latlon[2]]), c(mp[1],data[i,latlon[1]]),
            col=color0[data[i,color]])
    }
  }
  if(show.mean){
     points(mp[2], mp[1], col="black", pch=20, cex=1)
  }
  cd <- par("usr")
  legend(x=cd[1]+0.0001*(cd[2]-cd[1]), y=cd[3]+0.55*(cd[4]-cd[3]),
        legend=names(color0), fill=color0,bg="gray100", cex=0.6)
}

#------------------------------------------------
# Get information on City coordinades for WS from GOOGLE
#------------------------------------------------
get_rev_geocode <- function(data, latlon){
  res <- data[]
  res$tmp1 <- data[,latlon[1]]
  res$tmp2 <- data[,latlon[2]]
  res <- tidygeocoder::reverse_geocode(res, lat=tmp1, long=tmp2,
                          method='arcgis', full_results=TRUE)

  res <- data.frame(res[,c(colnames(data),"City","Region","RegionAbbr")])
  res$city2 <- gsub("Town of |City of ","",res$City)
  res$state2 <- res[,grep("RegionAbbr",colnames(res))]
  res <- res[,-grep("^City|RegionAbbr",colnames(res))]
  res$city_state2 <- paste0(res$city2,",",res$state2)
  res
}

#------------------------------------------------
# Next two functions are used to search for
# similar names.
# 'Find string' will return the element in x2 that is
# the same (except for some missing characters) to x1
# Example: x1='computer',  x2=c('compter','compoter')
#         find_string(x1,x2) retuns 'compter'
#------------------------------------------------
find_string <- function(x1,x2){
  xx1 <- strsplit(tolower(x1),"")[[1]]
  xx1 <- unique_split(xx1)
  cont <- 0
  res <- c()
  for(k in 1:length(x2)){
    xx2 <- unique_split(strsplit(tolower(x2[k]),"")[[1]])
    if(length(xx1) < length(xx2)){
     tmp <- match(xx1, xx2)
    }else{
     tmp <- match(xx2, xx1)
    }
    if(all(!is.na(tmp))){
      if(all(diff(tmp)>0)){
        cont <- cont +1
        res[cont] <- x2[k]
      }
    }
  }
  res[res!=x1]
}

unique_split <- function(x){
  out <- x[]
  tmp <- table(x)
  for(j in 1:length(x)){
    tt <- as.vector(tmp[x[j]])
    out[x==x[j]] <- paste0(x[x==x[j]],1:tt)
  }
  out
}
