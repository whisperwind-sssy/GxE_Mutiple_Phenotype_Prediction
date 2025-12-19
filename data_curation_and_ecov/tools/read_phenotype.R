
# factors_list=NULL; dates_list=NULL; traits_list=NULL; verbose=TRUE
read_phenotype <- function(pheno_folder,
                           factors_list=NULL,
                           dates_list=NULL,
                           traits_list=NULL,
                           verbose=TRUE){

  # Files that contain phenotype data
  files <- list.files(pheno_folder)
  stopifnot(length(files)>0)

  years <- unlist(lapply(strsplit(files,"_"),function(x)x[2]))
  if(any(nchar(years)!=4) | any(is.na(as.numeric(years)))){
    stop("Year information in '.csv' files names must be of the YYYY format")
  }

  # Read files
  if(verbose){
    message("-> Reading ",length(files)," phenotypic files:")
    message("   ",paste(years,collapse=','))
  }
  DATA <- vector("list", length(files))
  n0 <- 0
  for(i in 1:length(files)){
    DATA[[i]] <- data.table::fread(paste0(pheno_folder,'/', files[i]),
                                   stringsAsFactors=FALSE,
                                   data.table=FALSE)
    colnames(DATA[[i]]) <- tolower(colnames(DATA[[i]]))
    DATA[[i]]$year <- years[i]
    n0 <- n0 + nrow(DATA[[i]])
  }
  names(DATA) <- years
  # lapply(DATA,function(x)x[1:5,1:4])

  # Rename columns
  if(is.null(factors_list)){
    factors_list <- list(
      location =        c('field.*location'),
      genotype =        c('pedigree'),
      source =          c('source'),
      rep =             c('replicate')
    )
  }

 if(is.null(dates_list)){
    dates_list <- list(
      date_plant =      c('date.*planted'),
      date_harvest =    c('date.*harvest'),
      date_anthesis =   c('anthesis.*mm/dd/yy|anthesis.*date'),
      date_silking =    c('silking.*mm/dd/yy|silking.*date')
    )
  }

  if(is.null(traits_list)){
    traits_list <- list(
      plot_area =       c('plot.*area'),
      seed_number =     c('seed'),
      plants_stand =    c('stand.*plants'),
      grain_moisture =  c('grain.*moisture'),
      yield =           c('grain.*yield')
    )
  }

  factors <- names(factors_list)
  if(length(grep("year",factors))==0)  factors <- c("year",factors)
  dates <- names(dates_list)
  traits <- names(traits_list)
  all_names <- c(factors,dates,traits)

  if(length(grep("yield",tolower(traits)))!=1L){
    stop("One 'yield' name must be provided in 'traits_list'")
  }
  if(length(grep("location",tolower(factors)))!=1L){
    stop("One 'location' name must be provided in 'factors_list'")
  }

  if(verbose){
    message("-> Obtaining common columns from files ...")
    message("   ",length(factors)," factors: ",paste0(factors,collapse=", "))
    message("   ",length(dates)," dates: ",paste0(dates,collapse=", "))
    message("   ",length(traits)," traits: ",paste0(traits,collapse=", "))
  }
  for(i in 1:length(DATA)){
   colnames(DATA[[i]]) <- Rename(colnames(DATA[[i]]),
                                 factors_list, replace="all", verbose=FALSE)
   colnames(DATA[[i]]) <- Rename(colnames(DATA[[i]]),
                                 dates_list, replace="all", verbose=FALSE)
   colnames(DATA[[i]]) <- Rename(colnames(DATA[[i]]),
                                 traits_list, replace="all", verbose=FALSE)

   missing_col <- all_names[!all_names %in% colnames(DATA[[i]])]
   if(length(missing_col)>0){
     tmp <- c(ifelse(any(factors%in%missing_col),"factors_list",NA),
              ifelse(any(dates%in%missing_col),"dates_list",NA),
              ifelse(any(traits%in%missing_col),"traits_list",NA))
     tmp <- paste(tmp[!is.na(tmp)], collapse=",")
     stop("Some column names could not be found for year ",names(DATA)[i],":",
          "\n\t",paste(paste0("'",missing_col,"'"), collapse=","),
          "\n\t","Change argument(s) ",tmp," appropriately")
   }
  }

  # merge datasets from different years
  # Convert to character dates and merge all years
  if(verbose){
    message("-> Merging ",n0," records from all ",length(DATA)," years ...")
  }
  pheno <- do.call(rbind, lapply(DATA, function(x){
    out <- x[,all_names]
    for(k in seq_along(dates)){
      # message(out$year[1],": Class ",dates[k],"='",class(out[,dates[k]]),"'")
      if(any(c("IDate","Date") %in% class(out[,dates[k]])) )
      { # format to mm/dd/yyy
        #dt <- cbind(month(out[,dates[k]]),mday(out[,dates[k]]),year(out[,dates[k]]))
        dt <- cbind(lubridate::month(out[,dates[k]]),
                    lubridate::mday(out[,dates[k]]),
                    lubridate::year(out[,dates[k]]))
        new_dt <- apply(dt,1,paste,collapse="/")
        iNA <- which(apply(dt,1,function(z)all(is.na(z))))
        if(length(iNA)>0) new_dt[iNA] <- NA
        out[,dates[k]] <- new_dt
      }else{
        out[,dates[k]] <- as.character(out[,dates[k]])
      }
    }
    out
  }))
  rownames(pheno) <- NULL
  head(pheno)

  # Remove Germany and Ontario, Canada
  location_name <- colnames(pheno)[grep("location",tolower(colnames(pheno)))]
  pheno$location <- pheno[,location_name]

  tmp <- c('GEH','ONH')
  index <- grep(paste0(tmp,collapse="|"),pheno$location)
  if(length(index) >0) pheno <- pheno[-index,]

  # Transform variable types
  for(k in seq_along(factors)){
    pheno[,factors[k]] <- factor(as.character(pheno[,factors[k]]))
  }

  for(k in seq_along(traits)){
    if(!storage.mode(pheno[,traits[k]]) %in% c("double","integer")){
      stop("Trait '",traits[k],"' is not of the numeric type")
    }
    pheno[,traits[k]] <- as.numeric(pheno[,traits[k]])
  }

  # Fix date format
  if(verbose){ message("-> Treating columns with 'date' format ...") }
  for(k in seq_along(dates)){
    pheno[,dates[k]] <- fix_date(pheno[,dates[k]])
  }

  # Remove empty rows
  index <- apply(pheno[,factors],1,function(x){
    sum(is.na(x) | nchar(as.character(x))==0) == length(x)
  })
  if(any(index)){
    if(verbose){
      message("-> Removed ",length(which(index))," fully empty or NA rows with factors ...")
    }
    pheno <- pheno[-which(index),]
  }

  # Remove rows with no date information
  index <- rowSums(is.na(pheno[,dates])) == length(dates)
  if(any(index)){
    if(verbose){
      message("-> Removed ",length(which(index))," records with no date information ...")
    }
    pheno <- pheno[-which(index),]
  }

  # Remove rows with no yield information
  yield_name <- colnames(pheno)[grep("yield",tolower(colnames(pheno)))]
  index <- is.na(pheno[,yield_name] )
  if(any(index)){
    if(verbose){
      message("-> Removed ",length(which(index))," records with NA for yield ...")
    }
    pheno <- pheno[-which(index),]
  }

  # Change yield measure to kg/ha
  pheno[,yield_name] <- pheno[,yield_name] * 62.77
  pheno[,yield_name] <- pheno[,yield_name]/1000    # Tons/ha

  pheno$year_loc <- paste0(pheno$year,"-",pheno$location)
  pheno <- pheno[,c("year_loc",all_names)]

  if(verbose){
    message("-> Output: Phenotypic data with ",nrow(pheno)," rows and ",ncol(pheno)," columns")
    message("-> Done !")
  }
  return(pheno)
}
