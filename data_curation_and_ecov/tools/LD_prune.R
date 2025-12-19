

# Load 'Prune' function from SFSI R-package (Lopez-Cruz & de los Campos)
source("https://raw.githubusercontent.com/MarcooLopez/SFSI/main/R/Prune.R")

#================================================

#================================================
LD_prune <- function(X, MAP, threshold=0.95, n.windows=NULL, window.size=NULL,
                     d.max=NULL, mc.cores=1L, verbose=TRUE)
{
    colnames(MAP) <- toupper(colnames(MAP))
    CHRs <- sort(unique(MAP$CHR))

    if(is.null(colnames(X))){
      colnames(X) <- paste0("S",MAP$CHR,"_",MAP$POS)
    }
    index <- NULL
    for(j in 1:ncol(MAP)){
      if(all(MAP[,j]==colnames(X))) index <- j
    }
    if(is.null(index)){
      MAP <- data.frame(NAME=colnames(X),MAP)
    }else{
      colnames(MAP)[index] <- "NAME"
    }

    if(!is.null(n.windows) & !is.null(window.size)){
      window.size <- NULL
      warning("Input 'window.size' is ignored when also 'n.windows' is provided",
              .immediate=TRUE)
    }

    if("MAF" %in% colnames(MAP)){
      message("SNP with the highest MAF will be selected from each group sharing an R2>threshold")
    }else{
      message("First-appeared SNP will be selected from each group sharing an R2>threshold")
    }

    if(any(is.na(X))){
      message("Correlation is obtained from complete pairwise observations")
    }

    for(chr in CHRs){
      index <- which(MAP$CHR==chr)
      tmp <- sprintf('%.2f',diff(range(MAP$POS[index]))/1E6)
      cat(" - Chromosome ",chr,": size=",tmp," Mb. ",length(index)," SNPs\n",sep="")
    }

    if(!is.null(d.max)) d.max <- as.integer(d.max)

    compApply <- function(i)
    {
      chr <- CHRs[i]
      indexCHR <- which(MAP$CHR==chr)

      if(is.null(n.windows) & is.null(window.size)){
        n.windows0 <- 1L
        window.size0 <- length(indexCHR)
      }else{
        if(is.null(n.windows)){
          n.windows0 <- ceiling(length(indexCHR)/window.size)
          window.size0 <- table(rep(1:n.windows0, each=window.size)[1:length(indexCHR)])

          if(window.size0[n.windows0] < window.size*0.30){ # If the last window is very small
            window.size0[n.windows0-1] <- window.size0[n.windows0-1] + window.size0[n.windows0]
            window.size0 <- window.size0[1:(n.windows0-1)]
            n.windows0 <- n.windows0 - 1
          }

        }else{
          if(is.null(window.size)){
            n.windows0 <- n.windows
            window.size0 <- table(rep(1:n.windows,ceiling(length(indexCHR)/n.windows))[1:length(indexCHR)])
          }
        }
      }
      stopifnot(length(window.size0) == n.windows0)
      stopifnot(sum(window.size0) == length(indexCHR))
      window.size.cum <- c(0, cumsum(window.size0))

      id0 <- c()
      out0  <- vector('list', length(threshold))
      nConn0 <- vector('list', length(threshold))
      #names(out0) <- names(nConn0) <- threshold
      for(k in 1:n.windows0)
      {
        tmp <- c(window.size.cum[k] + 1, window.size.cum[k+1])
        index0 <- indexCHR[seq(tmp[1], tmp[2])]
        if(any(is.na(index0))){ stop("Marker matrix could not be subset at window ",k) }

        if(any(is.na(X[,index0]))){
          R0 <- cor(X[,index0], use="pairwise.complete.obs")^2
        }else{
          R0 <- (crossprod(scale(X[,index0]))/(nrow(X) - 1))^2
        }

        if(!is.null(d.max)){
          # D0 <- as.matrix(dist(MAP[index0,"POS",drop=F], method='manhattan'))
          D0 <- sapply(MAP$POS[index0],function(x)abs(x-MAP$POS[index0]))
        }else {
          D0 <- NULL
        }

        if("MAF" %in% colnames(MAP)){
          MAF0 <- MAP$MAF[index0]
        }else{
          MAF0 <- NULL
        }
        for(tr in seq_along(threshold)){
          ans <- Prune(R0, threshold=threshold[tr], D=D0, d.max=d.max, MAF=MAF0, verbose=FALSE)
          ##tt <- match(ans$prune.in,colnames(R))
          ##A <- (R[tt,tt] > threshold) & (D[tt,tt] <= d.max); diag(A) <- FALSE
          ##range(colSums(A))
          out0[[tr]] <- rbind(out0[[tr]], data.frame(CHR=chr,window=k,NAME=ans$prune.in))
          nConn0[[tr]] <- c(nConn0[[tr]], ans$nConn)

          message("    CHR=",sprintf('%2d',chr),". Window ",ifelse(n.windows0==1,NA,k)," [",
                  sprintf('%3d',tmp[1]),",",tmp[2],"]: ",sprintf('%5d',length(ans$prune.in)),
                  " SNPs selected with R2=",threshold[tr])
        }
        id0 <- c(id0, MAP[index0,"NAME"])
      }
      for(tr in seq_along(threshold)){
        out0[[tr]] <- out0[[tr]][order(match(out0[[tr]]$NAME, MAP$NAME[indexCHR])),]
        rownames(out0[[tr]]) <- NULL
        if(n.windows0 == 1L){
          out0[[tr]]$window <- NA
        }
      }

      return(list(pruneIn=out0, id=id0, nConn=nConn0))
    }

    if(mc.cores == 1L){
       res = lapply(X=seq_along(CHRs),FUN=compApply)
    }else{
       res = parallel::mclapply(X=seq_along(CHRs),FUN=compApply,mc.cores=mc.cores)
    }
    if(any(unlist(lapply(res, function(x)x$id)) != MAP$NAME)){
      stop("Something went wrong during the prunning procedure failed")
    }

    pruneIn <- lapply(seq_along(threshold), function(tr){
      do.call(rbind, lapply(res, function(x)x$pruneIn[[tr]]))
    })

    nConn <- lapply(seq_along(threshold), function(tr){
      unlist(lapply(res, function(x)x$nConn[[tr]]))
    })
    names(pruneIn) <- names(nConn) <- threshold

    tmp <- unlist(lapply(pruneIn, nrow))
    for(i in 1:length(tmp)){
     message(" - ",sprintf('%5d',tmp[i])," of ",ncol(X)," SNPs selected across the ",
             length(CHRs)," CHRs for R2=",names(tmp)[i])
    }

    return(list(pruneIn=pruneIn, nConn=nConn))
}
