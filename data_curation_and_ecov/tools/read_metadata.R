

# factors_list=latlon_list=latlon_field_list=NULL
# dMax = 15; location_changes=NULL; verbose=TRUE
read_metadata <- function(metadata_folder,
                          agronomic_folder,
                          factors_list=NULL,
                          latlon_list=NULL,
                          latlon_field_list=NULL,
                          dMax = 15,
                          location_changes=NULL,
                          verbose=TRUE){

  verbose <- as.numeric(verbose)

  files <- list.files(metadata_folder)
  stopifnot(length(files)>0)

  years <- unlist(lapply(strsplit(files,"_"),function(x)x[2]))
  if(any(nchar(years)!=4) | any(is.na(as.numeric(years)))){
      stop("Year information in '.csv' files names must be of the YYYY format")
  }

  if(verbose>0){
    message("-> Reading ",length(files)," metadata files:")
    message("   ",paste(years,collapse=','))
  }
  DATA <- vector("list", length(files))
  for(i in 1:length(files)){
    DATA[[i]] <- data.table::fread(paste0(metadata_folder,'/',files[i]),
                                   stringsAsFactors=FALSE,
                                   data.table=FALSE)

    # lower case all column names, remove character #### from colnames
    colnames(DATA[[i]]) <- tolower(gsub("#+","",colnames(DATA[[i]])))
    DATA[[i]]$year <- years[i]
  }
  names(DATA) <- years

  # Rename columns
  if(is.null(factors_list)){
    factors_list <- list(
      location =  c('experiment'),
      city =      c('city'),
      type  =     c('^type$'),
      treatment = c('treatment'),
      farm =      c('farm'),
      field =     c('^field$'),
      irrigated = c('^irrigated'),
      wsn  =      c('wsn|ws_sn|^weather.*serial.*number')
    )
  }

  if(is.null(latlon_list)){
    latlon_list <- list(
      latitude  =  c('^lat$|ws.*lat|^weather.*latitude'),
      longitude  =  c('^long$|ws.*lon|^weather.*longitude')
    )
  }

  if(is.null(latlon_field_list)){
    latlon_field_list <- list(
      'latitude_*'  =  c('corner.*.[1-4]*lat|latitude.*corner.*[1-4]'),
      'longitude_*' =  c('corner.*[1-4].*lon|longitude.*corner.*[1-4]')
    )
  }

  # Reading data from all years ...
  for(i in 1:length(DATA)){
    colnames(DATA[[i]]) <- Rename(colnames(DATA[[i]]),
                                   factors_list, replace='all',
                                   verbose=FALSE)

    colnames(DATA[[i]]) <- Rename(colnames(DATA[[i]]),
                                    latlon_list, replace='all',
                                    verbose=FALSE)

    colnames(DATA[[i]]) <- Rename(colnames(DATA[[i]]),
                                   latlon_field_list, replace='all',
                                   verbose=FALSE)
  }

  # Merge datasets from different years, insert NA when the column is not found
  factors <- names(factors_list)
  if(length(grep("year",factors))==0)  factors <- c("year",factors)
  latlon <- names(latlon_list)
  latlon_field <- unlist(lapply(names(latlon_field_list),function(x)
                       lapply(1:4, function(j)gsub("\\*",j,x))))

  if(length(grep("location",tolower(factors)))!=1L){
     stop("One 'location' name must be provided in 'factors_list'")
  }
  if(length(grep("type",tolower(factors)))!=1L){
     stop("One 'type' name must be provided in 'factors_list'")
  }
  if(length(grep("irrigated",tolower(factors)))!=1L){
     stop("One 'irrigated' name must be provided in 'factors_list'")
  }

  new_names <- c(factors,latlon, latlon_field)
  META <- do.call(rbind, lapply(DATA, function(x){
     missing <- new_names[!new_names %in% colnames(x)]
     if(length(missing) > 0){
       tmp <- data.frame(matrix(NA,ncol=length(missing),nrow=nrow(x)))
       colnames(tmp) <- missing
       if("type" %in% missing){tmp$type <- "hybrid"}
       x <- data.frame(x,tmp)
     }else{
        x <- data.frame(x)
     }
     x[,new_names]
   }))
  rownames(META) <- NULL
  #head(META)

  #------------------------------------------------
  irri_name <- colnames(META)[grep("irrigated",tolower(colnames(META)))]
  location_name <- colnames(META)[grep("location",tolower(colnames(META)))]
  year_name <- colnames(META)[grep("year",tolower(colnames(META)))]
  city_name <- colnames(META)[grep("city",tolower(colnames(META)))]

  META$irrigated <- META[,irri_name]
  META$location <- META[,location_name]
  META$year <- META[,year_name]
  META$city <- META[,city_name]

  #------------------------------------------------
  # Getting information on irrigation ...
  #------------------------------------------------
  if(verbose>0){ message("-> Getting information on irrigation ...") }
  irri_values <- c("irrigation","fertigation")

  META$irrigated <- tolower(META$irrigated)
  files <- list.files(agronomic_folder)
  for(i in 1:length(years)){
    index <- which(META[,year_name]==years[i])
    irri <- META[index,'irrigated']
    if(verbose>1){ message("  [",i,"]\tYear ",years[i],":") }
    if(any(is.na(irri))){  # No IRRIGATED column was found
      if(verbose>1){
       message("   No irrigation information was found in Metadata")
      }
      filename <- grep(years[i],files,value=T)
      if(length(filename)>0){
        META[index[is.na(irri)],'irrigated'] <- "no"
        dat <- data.table::fread(paste0(agronomic_folder,'/',filename),data.table=FALSE)
        colnames(dat) <- tolower(colnames(dat))
        colnames(dat) <- Rename(colnames(dat),list(location="location|experiment"),
                                 replace="all", verbose=FALSE)
        colnames(dat) <- Rename(colnames(dat),list(application="application.*or.*treatment"),
                                replace="all", verbose=FALSE)
        trt <- tolower(dat[,"application"])
        locs <- unique(dat[trt %in% irri_values,"location"])
        if(length(locs)>0){
          if(verbose>1){
            message("   Information on ",paste(irri_values,collapse="|")," for ",
                  length(locs)," locations was found in Agronomic data")
          }

          ii <- META[index,"location"] %in% locs
          META[index[ii],'irrigated'] <- "yes"
        }
      }
    }else{
      META[index[grep("no",irri)],'irrigated'] <- "no"
      META[index[irri == ""],'irrigated'] <- NA
      tmp <- unique(META[index[which(META[index,'irrigated']=="yes")],"location"])
      tmp <- tmp[tmp != ""]
      if(verbose>1){ message("   Irrigation was found in Metadata for ",length(tmp)," locations") }
      tmp <- unique(META[index[which(is.na(META[index,'irrigated']))],"location"])
      if((verbose>1) &(length(tmp)>0)){
        message("   ",length(tmp)," locations with UNKNOWN irrigation: ",
                paste(tmp,collapse=","))
      }
    }
  }

  # Searching for experiments from hybrids only
  type_name <- colnames(META)[grep("type",tolower(colnames(META)))]
  META[,type_name] <- tolower(META[,type_name])
  index <- META[,type_name] == "hybrid"
  if(any(index)){
    if(verbose>0){
      message("-> Found metadata for ",length(which(index))," year-locations with 'hybrids' experiments...")
    }
    META <- META[which(index),]
  }

  # Standardize location names
  STATES <- data.frame(state=state.abb, name=state.name, region=state.region)

  META$location <- gsub(" ","",META$location)
  META$location2 <- gsub("G2F","",META$location)
  META$state <- substr(META$location2,1,2)

  # Keep only US states
  index <- grep("^[A-Z]{2}H",META$location2)  # Format XXH
  drop <- index[!substr(META$location2[index],1,2)%in%STATES$state]
  if(length(drop)>0){
    META <- META[-drop,]
  }

  tt <- split(META,paste(META$year,META$state))
  drop <- c()
  for(i in 1:length(tt)){
    tt0 <- tt[[i]]
    for(j in 1:nrow(tt0)){
      tmp <- regexpr(paste0("^",tt0$state[1],"H[1-9]+"), tt0$location2[j])
      tt0$location2[j] <- ifelse(tmp[1]>0,regmatches(tt0$location2[j],tmp),tt0$location2[j])
    }

    tmp <- grep(paste0(tt0$state[1],"H[1-9]+"),tt0$location2,value=T)
    flag <- ifelse(length(tmp)==nrow(tt0), all(tmp == tt0$location2),FALSE)
    if(flag){
      if(!tt0$state[1] %in% STATES$state){
        if(verbose>0){
         message("   [",i,"]\tYear ",tt0$year[1],": Location(s): ",paste(tt0$location2, collapse=","),
                 " are not US locations")
        }
        drop <- c(drop, which(META$year==tt0[1,'year'] & META$location%in%tt0[,'location']))
      }
    }else{
      index <- which(!tt0$location2%in%tmp)
      loc0 <- tt0$location2[index]
      if((length(index)==1) & (tt0$state[1] %in% STATES$state)){
        if(length(grep(paste0(tt0$state[1],"[A-Z][1-9]"),loc0)) == 1){
          loc1 <- gsub(paste0(tt0$state[1],"[A-Z]"),paste0(tt0$state[1],"H"),loc0)
        }else{
          loc1 <- paste0(tt0$state[1],"H1")
        }
        index0 <- which(META$year==tt0[index,'year'] & META$location==tt0[index,'location'])
        META[index0,'location'] <- loc1
        if(verbose>0){
         message("   [",i,"]\tYear ",tt0[index,'year'],":",
                 " Location name '",tt0[index,'location'],"' replaced to '",loc1,"'")
        }
      }else{
        if(verbose>0){
         message("   [",i,"]\tYear ",tt0$year[1],": Location(s): ",paste(tt0$location2, collapse=","),
                 " are not US locations with format XXH1, XXH2, etc")
        }
        drop <- c(drop, which(META$year==tt0[1,'year'] & META$location%in%tt0[index,'location']))
      }
    }
  }
  if(length(drop)>0){
    if(verbose>1){
     message("   New location names can be set by passing argument 'location_changes' as")
     message("     location_changes <- list(old_name1='new_name1', old_name2='new_name2') ")
    }
  }


  if(!is.null(location_changes)){
    if(verbose){ message("-> Some entries in 'location' column were replaced:") }
    state <- substr(META$location,1,2)
    META$location <- Rename(META$location, location_changes,
                                   multiple.patterns=TRUE,
                                   multiple.existing=TRUE,
                                   replace="all",
                                   verbose=TRUE)
  }

  # Standardizing some location and city names ...
  META$year_loc <- paste0(META$year,"-",META$location)
  META$state <- substr(META$location,1,2)

  # City information
  META$city <- unlist(lapply(strsplit(META$city,","),function(x)x[1]))
  META$city <- gsub("-"," ",META$city)
  META$city <- capitalize(META$city)
  META$city_state <- paste0(META$city,",",META$state)

  # Remove duplicated rows
  index <- paste0(META$year_loc,"_",META$treatment,"_",META$city)
  tt <- table(index); tt2 <- names(tt[tt>1])
  if(length(tt2)>0){
    drop <- c()
    for(i in 1:length(tt2)){
      drop <- c(drop, which(index == tt2[i])[-1])
    }
    META <- META[-drop, ]
    index <- index[-drop]
  }
  stopifnot(nrow(META) == length(unique(index)))

  # Obtaining average field coordinates
  META$latitude_field <- apply(META[,grep("latitude_",latlon_field,value=T)],1,
    function(x) mean(as.numeric(gsub("\n","",x)),na.rm=T)
  )
  META$longitude_field <- apply(META[,grep("longitude_",latlon_field,value=T)],1,
    function(x) mean(as.numeric(gsub("\n","",x)),na.rm=T)
  )
  META <- META[, -grep("latitude_[1-4]",colnames(META))]
  META <- META[, -grep("longitude_[1-4]",colnames(META))]
  latlon_field <- c("latitude_field","longitude_field")

  #------------------------------------------------
  # Get median latlon for unique cities
  #------------------------------------------------
  names0 <- c("city","state","city_state",latlon)
  index1 <- apply(META[,latlon],1,function(x)all(!is.na(x)))
  index2 <- apply(META[,latlon_field],1,function(x)all(!is.na(x)))
  index <- which((!index1)&index2)
  if(length(index)>0){
    tmp <- META[index,latlon_field]; colnames(tmp) <- latlon
    META[index,latlon] <- tmp
    META[index,latlon_field] <- NA # Set latlon_field to NA
  }
  index1 <- apply(META[,latlon],1,function(x)all(!is.na(x)))
  META0 <- META[index1, names0, drop=FALSE]
  META0 <- do.call(rbind,lapply(split(META0, META0$city_state),function(x){
    tt <- do.call(cbind,apply(x[,latlon],2,as.numeric,simplify=F))
    data.frame(x[1,names0[!names0%in%latlon],drop=F],
               n=nrow(tt), t(apply(tt,2,median)))
  }))
  rownames(META0) <- NULL
  #head(META0)

  if(verbose>0){ message("-> Obtaining correct city names from latlon using ArcGIS ...") }
  INFO <- get_rev_geocode(META0,latlon)
  #head(INFO)

  META$city_new <- META$city
  META$status <- NA
  for(i in 1:nrow(META)){
    meta <- META[i,]
    cont <- 0
    prefix <- paste0("   [",i,"]\t",meta$year_loc,":")

    a1 <- paste0("[",paste(sprintf("%.4f",meta[latlon]),collapse=","),"] ")
    a2 <- paste0("[",paste(sprintf("%.4f",meta[latlon_field]),collapse=","),"] ")
    a1 <- ifelse(verbose<2,"",a1)
    a2 <- ifelse(verbose<2,"",a2)

    if(all(!is.na(meta[latlon]))){
      META$status[i] <- "OK"
      INFO$d1 <- rep(NA, nrow(INFO))
      INFO$d2 <- rep(NA, nrow(INFO))
      for(j in 1:nrow(INFO)){
        INFO$d1[j] <- get_distance(meta[latlon], INFO[j,latlon])/1000
        INFO$d2[j] <- get_distance(meta[latlon_field], INFO[j,latlon])/1000
      }
      info <- INFO[INFO$city_state==meta$city_state,]

      # See if it can be changed for another city with similar name
      tmp <- find_string(meta$city,INFO$city)
      info2 <- INFO[INFO$city%in%tmp,]
      flag <- ifelse(nrow(info2)==0, FALSE,
                     info2[which.max(info2$n),"n"]>info$n)
      if(flag){
        if(info2$state == info2$state2){ # META$city found by ArcGIS within state
          META[i,"city_new"] <- info2$city
          info <- info2
          if(verbose>0){
           message(ifelse(cont>0,"\t\t",prefix)," City name ",meta$city_state,
                   " replaced to ",info2$city_state)
          }
          cont <- cont + 1
        }else{
          META$status[i] <- "CHECK"
          #message("   City ",meta$city_state,": SUSPICIOUS!!!!!!")
        }
      }

      a3 <- paste0("[",paste(sprintf("%.4f",info[latlon]),collapse=","),"] ")
      a3 <- ifelse(verbose<2,"",a3)

      if(info$city_state==info$city_state2){ # META$city is the same found by ArcGIS
        if(info$d1 <= dMax){
          # message("   City ",meta$city_state,": OK")
        }else{
          if(ifelse(is.na(info$d2),FALSE,info$d2<dMax)){  # Check field coordinates
            info2 <- get_rev_geocode(meta, latlon_field)  # Geocode for latlon field
            if((meta$city != info2$city2) & (meta$state==info2$state2) & (info2$city2!="")){
              META[i,"city_new"] <- info2$city2
              if(verbose>0){
                message(ifelse(cont>0,"\t\t",prefix)," City name ",meta$city_state,
                        " replaced to ",info2$city_state2)
              }
              cont <- cont + 1
            }
            META[i,latlon] <- META[i,latlon_field]
            if(verbose>0){
              message(ifelse(cont>0,"\t\t",prefix)," City ",meta$city_state,": latlon ",a1,
                      "changed to field latlon ",a2)
            }
            cont <- cont + 1
          }else{  # Or set to the median value
            META[i,latlon] <- info[latlon]
            if(verbose>0){
              message(ifelse(cont>0,"\t\t",prefix)," City ",meta$city_state,": latlon ",a1,
                      "changed to the median ",a3)
            }
            cont <- cont + 1
          }
        }
      }else{
        if(meta$state==info$state2){
          flag_change <- TRUE
          if(info$d1 > dMax){
            if(ifelse(is.na(info$d2),FALSE,info$d2<dMax)){  # Check field coordinates
              info2 <- get_rev_geocode(meta, latlon_field)  # Geocode for latlon field
              if((meta$city != info2$city2) & (meta$state==info2$state2) & (info2$city2!="")){
                META[i,"city_new"] <- info2$city2
                if(verbose>0){
                  message(ifelse(cont>0,"\t\t",prefix)," City name ",meta$city_state,
                          " replaced to ",info2$city_state2)
                }
                cont <- cont + 1
              }
              flag_change <- FALSE  # do not change later because we use now field coordinates
              META[i,latlon] <- META[i,latlon_field]
              if(verbose>0){
                message(ifelse(cont>0,"\t\t",prefix)," City ",meta$city_state,": latlon ",a1,
                        "changed to field latlon ",a2)
              }
              cont <- cont + 1
            }else{  # Or set to the median value
              META[i,latlon] <- info[latlon]
              if(verbose>0){
               message(ifelse(cont>0,"\t\t",prefix)," City ",meta$city_state,": latlon ",a1,
                       "changed to the median ",a3)
              }
              cont <- cont + 1
            }
          }

          if((info$city2 != "") & flag_change){
            META[i,"city_new"] <- info$city2
            if(verbose>0){
              message(ifelse(cont>0,"\t\t",prefix)," City name ",meta$city_state," replaced to ",info$city_state2)
            }
            cont <- cont + 1
          }

        }else{
          # Check field coordinates first
          flagOK <- FALSE
          if(all(!is.na(meta[latlon_field]))){
            info2 <- get_rev_geocode(meta, latlon_field)  # Geocode for latlon field
            if(info$state == info2$state2){
              META[i,latlon] <- META[i,latlon_field]
              if(verbose>0){
                message(ifelse(cont>0,"\t\t",prefix)," City ",meta$city_state,": latlon ",
                        a1,"changed to field latlon ",a2)
              }
              cont <- cont + 1
              if((meta$city != info2$city2) & (info2$city2!="")){
                META[i,"city_new"] <- info2$city2
                if(verbose>0){
                    message(ifelse(cont>0,"\t\t",prefix)," City name ",meta$city_state,
                            " replaced to ",info2$city_state2)
                }
                cont <- cont + 1
              }
              flagOK <- TRUE
            }
          }
          if(!flagOK){
            tmp <- tidygeocoder::geocode(meta[,c(names0,"city_state")], city=city, state=state, method='osm',full_results=TRUE)
            info2 <- get_rev_geocode(tmp, c("lat","long"))
            if(info2$city_state == info2$city_state2){
              META[i,latlon] <- info2[,c("lat","long")]
              if(verbose>0){
                message(ifelse(cont>0,"\t\t",prefix)," City name ",meta$city_state,
                        " latlon replaced by ArcGIS latlon")
              }
              cont <- cont + 1
              flagOK <- TRUE
            }
          }

          if(!flagOK){
            META$status[i] <- "CHECK"
            #message("   City ",meta$city_state," SUSPICIOUS!!!!!!")
          }
        }
      }
    }else{
      # See if city name must be changed
      if(!is.na(meta$city)){
        info <- INFO[INFO$city_state==meta$city_state,]
        if(info$state == info$state2){
          if((meta$city != info$city2) & (info$city2!="")){
            META[i,"city_new"] <- info$city2
            if(verbose>0){
             message(ifelse(cont>0,"\t\t",prefix)," City name ",meta$city_state,
                      " replaced to ",info$city_state2)
            }
            cont <- cont + 1
          }
        }
      }

      META$status[i] <- "MISSING"
      if(verbose>0){
       message(ifelse(cont>0,"\t\t",prefix)," City ",meta$city_state,": MISSING latlon. To be imputed")
      }
      cont <- cont + 1
    }
  }
  # Calculate again city-state
  #table(META$status)
  city_info <- META[which(META$city != META$city_new),c("city","city_new")]

  META$city <- META$city_new
  META$city_state <- paste0(META$city,",",META$state)

  #--------------------------------------------------
  ## Imputing missing coordinates (and) from available mean ...
  #--------------------------------------------------
  if(verbose>0){ message("-> Retrieving missing coordinates from available information ...") }
  index_missing <- which(META$status=="MISSING")

  for(i in 1:length(index_missing)){
    meta <- META[index_missing[i],]
    prefix <- paste0("   [",i,"]\t",meta$year_loc,":")
    # If the city name was changed previously
    if(meta$city %in% city_info[,1]){
      meta[,'city'] <- city_info[match(meta$city, city_info[,1]),2]
    }

    entry <- list(year=meta$year,location=meta$location,
                  city=meta$city, state=meta$state)
    res <- get_coordinates_mean(META, entry=entry, metric="median",all=TRUE)
    res <- res[res[,"n.city"]==1 & res[,"n.OK"]>0,]

    index <- which(rownames(res)%in%c("year_city","city"))
    if(length(index)==0){
      index <- which(rownames(res)=="year_location")
      if(length(index)==0){
        index <- which(rownames(res) == "year_state")
        if(length(index)==0){
          if(length(index)==0){
            index <- which(rownames(res)=="location")
            if(length(index)==0){
              index <- which(rownames(res)=="state")
            }
          }
        }
      }
    }

    city0 <- ifelse(is.na(meta$city), "MISSING", meta$city)
    if(length(index)>0){
      res0 <- res[index[1], , drop=FALSE]
      a0 <- paste0(" (n=",res0[,'n.OK'],") ")
      a1 <- paste0("[",paste(sprintf("%.4f",res0[latlon]),collapse=","),"]")
      a0 <- ifelse(verbose<2,"",a0)
      a1 <- ifelse(verbose<2,"",a1)

      tmp <- ""
      if(is.na(meta$city) & (!is.na(res0$city))){
       meta$city <- res0$city
       city0 <- res0$city
       tmp <- "city and "
      }

      meta[latlon] <- res0[latlon]
      meta$status <- "OK"
      META[index_missing[i],] <- meta
      if(verbose>0){
        message(prefix," City ",city0,": ",tmp,"latlon imputed from ",
                rownames(res0)," median",a0,a1)
      }
    }else{
      if(verbose>0){
        message(prefix," City ",city0,": No imputation can be made")
      }
    }
  }
  table(META$status)
  META$city_state <- paste0(META$city,",",META$state)

  # Select apropiate columns
  META <- META[,c('year_loc',factors,latlon)]

  # Final edits to file
  farm_name <- colnames(META)[grep("farm",tolower(colnames(META)))]
  if(length(farm_name)==1){
    META[,farm_name] <- gsub("\\,"," ",META[,farm_name])
  }
  field_name <- colnames(META)[grep("field",tolower(colnames(META)))]
  if(length(field_name)==1){
    META[,field_name] <- gsub("\\,","/",META[,field_name])
  }

  # Add region information
  META$region <- ifelse(META$latitude<37,"South","North")

  # Change Delaware to South
  META[META$state=="DE","region"] <- "South"

  META
}
