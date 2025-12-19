

significance_code <- function(x){
  #ifelse(x<=0.001,"***",ifelse(x<=0.01,"**",ifelse(x<=0.05,"*",ifelse(x<=0.1,".","ns"))))
  ifelse(x>0.1,"ns",ifelse(x>0.05,".",ifelse(x>0.01,"*",ifelse(x>0.001,"**","***"))))
}

get_HI <- function(temp, RH){
  # To test:  temp=90; RH=100
  # temp: air temperature (F)
  # RH: relative humidity (percentage)
  HI <- -42.379 + 2.04901523*temp + 10.14333127*RH - 0.22475541*temp*RH -
         (6.83783E-3)*(temp^2) - (5.481717E-2)*(RH^2) + (1.22874E-3)*(temp^2)*RH +
          (8.5282E-4)*temp*(RH^2) - (1.99E-6)*(temp^2)*(RH^2)

  as.numeric(HI)
}

CtoF <- function(x){
  (1.8*x)+32
}

FtoC <- function(x){
  (x-32)/1.8
}

#=============================================================================

get_ec_info <- function(EC){

  if(length(dim(EC)) == 2L){
    stopifnot(!is.null(colnames(EC)))
    namesEC <- colnames(EC)
  }else{
    if(is.character(EC)){
      namesEC <- EC
    }else{
      stop("'EC' must be either a matrix a character vector with EC names")
    }
  }

  out <- data.frame(do.call(rbind,lapply(strsplit(namesEC,"_"),function(x){
    if(length(x)==1){
       x <- c(x,NA,NA)
    }else{
      if(length(x)==2){
         x <- c(x,NA)
      }
    }
    if(length(which(x==""))>0) x[which(x=="")] <- NA
    x
  })))
  out <- data.frame(ec=namesEC,out)
  colnames(out) <- c("name","variable","period","layer")
  out$layer <- as.numeric(as.character(out$layer))

  # Add start/end
  if(length(grep("^TT\\(", out$period)) == nrow(out))
  {
    period0 <- gsub("^TT\\(|\\)","",out$period)
    out$start <- as.numeric(period0)*100 - 50
    out$end <- as.numeric(period0)*100 + 50

    plevels <- unique(period0)
    plevels <- plevels[order(as.numeric(plevels))]
    plevels <- paste0("TT(",plevels,")")
    out$period <- factor(as.character(out$period), levels=plevels)

  }else{
    stages <- c(Ger="Germination", Eme="Emergence", EnJ="EndJuvenile",
                Flo="FloralInitiation", Fla="FlagLeaf", Flw="Flowering",
                StG="StartGrainFill", EnG="EndGrainFill", Mat="Maturity",
                Har="HarvestRipe")

    a1 <- substr(out$period,1,3)
    a2 <- substr(out$period,4,6)
    if(!all(c(a1,a2) %in% names(stages))){
      tmp <- c(a1[!a1 %in% names(stages)], a2[!a2 %in% names(stages)])
      message("Some period names could not be found: ",paste(unique(tmp),collapse=","))
    }
    out$start <- as.vector(stages[match(a1, names(stages))])
    out$end <- as.vector(stages[match(a2, names(stages))])

    plevels <- paste0(names(stages)[-length(stages)],names(stages)[-1])
    out$period <- factor(as.character(out$period), levels=plevels)
  }

  return(out)
}

#=============================================================================

aggregate_ec <- function(EC, by=c("layer","period"), stat=c("sum","mean"))
{
  by <- match.arg(by)
  stat <- match.arg(stat)

  EC_info <- get_ec_info(colnames(EC))
  infoX <- split(EC_info, EC_info$variable)
  perioded <- unlist(lapply(infoX,function(x) 1-mean(is.na(x$period)) ))
  layered <- unlist(lapply(infoX,function(x) 1-mean(is.na(x$layer)) ))
  stopifnot(perioded %in% c(0,1))
  stopifnot(layered %in% c(0,1))
  perioded <- perioded==1
  layered <- layered==1

  ecovs <- unique(EC_info$variable)

  if(length(unique(perioded)) > 1L){
    stop("Some variables are perioded and others not. All variables must be of the same type")
  }
  isperioded <- ifelse(all(perioded),TRUE,ifelse(all(!perioded),FALSE,NA))

  if(by=="period"){
    if(any(layered)){
      message(sum(layered)," of ",length(ecovs)," ECOVs were found layered. ",sum(!layered),
              " ECOVs were aggregated across periods and")
      message(sum(layered)," ECOVs were aggregated across periods within layer")
    }else{
      message("All ",length(ecovs)," ECOVs were aggregated across periods")
    }
  }
  if(by=="layer"){
    if(any(layered)){
      message(sum(!layered)," of ",length(ecovs)," ECOVs are not layered. Only ",sum(layered),
              " ECOVs were aggregated across layers",ifelse(isperioded," (within period)",""))
    }else{
      message(" No ECOVs is layered. No action is performed")
    }
  }

  for(k in 1:length(ecovs)){
    ecov <- ecovs[k]
    index <- EC_info$variable == ecov
    EC_info0 <- EC_info[index,]

    if(by == "layer"){
      if(layered[ecov]){
        if(isperioded){
          w0 <- do.call(cbind,lapply(split(EC_info0,EC_info0$period),function(x){
            apply(EC[,x$name],1,stat)
          }))
          colnames(w0) <- paste0(ecov,"_",colnames(w0))
        }else{
          w0 <- t(t(apply(EC[,index],1,stat)))
          colnames(w0) <- ecov
        }
      }else{
        w0 <- EC[,index,drop=FALSE]
      }

    }else{ # by period
      if(layered[ecov]){
        w0 <- do.call(cbind,lapply(split(EC_info0,EC_info0$layer),function(x){
          apply(EC[,x$name],1,stat)
        }))
        colnames(w0) <- paste0(ecov,"__",colnames(w0))
      }else{
        w0 <- t(t(apply(EC[,index],1,stat)))
        colnames(w0) <- ecov
      }
    }

    if(k==1){
      W <- w0
    }else{
      W <- cbind(W, w0)
    }
  }

  W
}

#=============================================================================

get_SDR <- function(EC, EC_KL, check.small=FALSE){

  EC <- as.matrix(EC)
  EC_info <- get_ec_info(EC)

  periods <- levels(EC_info$period)

  Ws <- EC[,EC_info$variable == "ESW"]

  # Check if Ws is layered
  EC_info0 <- EC_info[EC_info$variable == "ESW",]
  layers <- unique(EC_info0$layer)
  if(any(is.na(layers))){
    stop("EC 'ESW' must be layered")
  }
  layers <- layers[order(layers)]

  Wssum <- matrix(NA, nrow=nrow(Ws),ncol=length(periods))
  dimnames(Wssum) <- list(rownames(Ws),paste0("ESW_",periods))
  for(i in 1:nrow(Ws)){
    EC_KL0 <- EC_KL[EC_KL$year_loc==rownames(Ws)[i],]
    if(!all(layers %in% EC_KL0$layer)){
      stop("Some layers were not found for '",rownames(Ws)[i],"'")
    }
    KL <- EC_KL0[match(layers, EC_KL0$layer),"KL"]

    for(j in seq_along(periods)){
      Wssum[i,j] <- sum(Ws[i,paste0("ESW_",periods[j],"_",layers)]*KL)
    }
  }

  a1 <- EC[,EC_info$variable == "Eo"]
  a2 <- EC[,EC_info$variable == "CoverGreen"]

  # Check some CoverGreen equal to zero and set to min(CG[CG>0])
  if(check.small){
    eps <- sqrt(.Machine$double.eps)
    for(i in 1:nrow(a2)){
      index <- which(a2[i,]<eps)
      if(length(index)>0){
        a2[i, index] <- min(a2[i,-index])
      }
    }
  }

  Wd <- 2*a1[,paste0("Eo_",periods)]*a2[,paste0("CoverGreen_",periods)]

  SDR <- as.matrix(Wssum/Wd)
  colnames(SDR) <- paste0("SDR_",periods)
  SDR

}

#=============================================================================

get_period_ec <- function(EC, variable, period, p.breaks=NULL,
                          bottom=TRUE, stat=c("sum","mean")){
  stat <- match.arg(stat)

  EC_info <- get_ec_info(EC)
  period_levels <- levels(EC_info$period)

  # Get SD
  STAGE <- EC_info[match(period_levels, EC_info$period),c("period","start","end")]

  indexEC <- which(EC_info$variable %in% variable)
  SI <- as.matrix(EC[,indexEC])
  EC_info0 <- EC_info[indexEC,]

  if(any(table(EC_info0$period)>1)){
    message("ECs are aggregated within each period")
  }

  out <- c()
  for(i in 1:nrow(STAGE)){
    stg <- STAGE[i,]
    #tmp <- paste0("SDR_",stg$period)
    #tmp <- unlist(lapply(strsplit(tmp,""),function(x)
    #   paste(ifelse(x %in% c("(",")"), paste0("\\",x),x),collapse="")))
    #sdr <- SDR[,grep(tmp,colnames(SDR),value=T),drop=F]
    si <- SI[,EC_info0$period == stg$period, drop=FALSE]
    si <- apply(si,1,sum)
    tt <- do.call(rbind,lapply(strsplit(names(si),"-"),function(x){
                   c(year=x[1],location=paste(x[2:length(x)],collapse="-"))
    }))
    tmp <- data.frame(x=i,year_loc=names(si), tt,
               stg[rep(1,length(si)),],SI=si)
    out <- rbind(out, tmp)
  }
  #  head(out)
  #  Get the SDR group
  dat_group <- do.call(rbind,lapply(split(out, out$year_loc),function(x){
    index <- which(x$period %in% period)
    if(stat == "mean"){
      SI <- mean(x[index,"SI"])
    }else{
      SI <- sum(x[index,"SI"])
    }
    data.frame(x[1,c("year_loc","year","location")], nEC=length(index), SI=SI)
  }))
  rownames(dat_group) <- NULL

  if(!is.null(p.breaks)){
    #quantile(dat_group$SDR, probs=c(0.2, 0.25, 0.33, 0.5, 0.66, 0.8))
    if(bottom){
      breaks0 <- round(quantile(dat_group$SI, probs=p.breaks),2)
    }else{
      breaks0 <- round(quantile(dat_group$SI, probs=1-p.breaks),2)
      names(breaks0) <- paste0(100*(p.breaks),"%")
    }

    if(length(breaks0)>1){
      labels0 <- unlist(lapply(1:(length(breaks0)-1),function(i)
        paste0("(",breaks0[i],",",breaks0[i+1],"]")))
    }else{
      labels0 <- NULL
    }
    SI_levels <- c(paste0("<=",breaks0[1]),labels0,paste0(">",breaks0[length(breaks0)]))
    rg <- range(dat_group$SI)
    breaks2 <- c(rg[1]-0.01*diff(rg),breaks0,rg[2]+0.01*diff(rg))
    dat_group$SI_level <- cut(dat_group$SI, labels=SI_levels, breaks=breaks2)
    if(any(is.na(dat_group$SI_level))){
      stop("Some entries could not be assigned to a SI level")
    }
    stopifnot(all(levels(dat_group$SI_level)==SI_levels))
  }else{
    breaks0 <- NULL
  }

  list(dat=dat_group, breaks=breaks0)
}


EWAS_fit <- function(pheno, W, PC_adjust=TRUE, nPC=5, alpha=0.05){

  WtW <- tcrossprod(W)
  WtW <- WtW/mean(diag(WtW))
  EVD <- eigen(WtW, symmetric=TRUE)
  indexPC <- which(EVD$values>1E-8)

  if(PC_adjust){
    PCnames = paste0("PC",1:nPC)
    PC <- sweep(EVD$vectors[,indexPC], 2, sqrt(EVD$values[indexPC]),FUN='*')
    PC <- PC[,1:nPC]
    dimnames(PC) <- list(rownames(WtW),PCnames)
  }else{
    PC <- NA
  }

  bg <- BGData(geno=W, pheno=data.frame(pheno,PC))

  nt <- meff(eigen=EVD$values[indexPC], method="galwey")
  out <- c()

  for(k in 1:ncol(pheno)){
    if(PC_adjust){
      form0 = paste(colnames(pheno)[k]," ~ 1+",paste(PCnames,collapse="+"))
    }else{
      form0 = paste(colnames(pheno)[k]," ~ 1")
    }
    fm <- GWAS(formula(form0), data=bg, method=c("lsfit","rayOLS")[1])

    tmp <- data.frame(trait=colnames(pheno)[k],ecov=rownames(fm),fm,
                      bonferroni=alpha/nt, check.names=F)
    out <- rbind(out, tmp, make.row.names=FALSE)
  }
  colnames(out)[grep("Pr\\(",colnames(out))] <- "Pvalue"
  out
}

#--------------------------------------
# Plot heatmap
#--------------------------------------
# PC_adjust=TRUE; nPC=5; xlab="ECOV"; ylab="-log10(Pvalue)";
# alpha=0.05; color_name=NULL; R2_name=expression(R^2)
# expand=0.015; rel_height=0.5; point.size=1
EWAS <- function(pheno, W, PC_adjust=TRUE, nPC=5, ecov_annotation=NULL,
                      alpha = 0.05, color="Type", xlab="ECOV",
                      ylab="-log10(P)", color_name=NULL, traits=colnames(pheno),
                      point.size=1,expand=0.015, plot=TRUE){

  stopifnot(all(rownames(W) == rownames(pheno)))
  ecov_names <- colnames(W)
  W <- scale(W)

  out <- EWAS_fit(pheno, W, PC_adjust=PC_adjust, nPC=nPC, alpha=alpha)

  out$ecov <- factor(as.character(out$ecov), levels=ecov_names)

  if(is.null(ecov_annotation)){
    out$color <- data[,color]
    color_name <- NULL
  }else{
    out$color <- ecov_annotation[match(out$ecov, rownames(ecov_annotation)),1]
    color_name <- colnames(ecov_annotation)[1]
  }
  out$color <- factor(out$color)

  dat0 <- out[out$trait %in% traits,]
  if(is.null(names(traits))){
    dat0$trait <- factor(as.character(dat0$trait), levels=traits)
  }else{
    dat0$trait <- names(traits)[match(dat0$trait,traits)]
    dat0$trait <- factor(as.character(dat0$trait), levels=names(traits))
  }

  if(is.null(ecov_annotation)){
    out <- out[, -grep("^color$",colnames(out))]
  }else{
    colnames(out) <- gsub("^color$",color_name,colnames(out))
  }

  # Get threshold for manhattan plot
  datthr <- do.call(rbind,lapply(split(dat0, dat0$trait), function(x){
    x[1,]
  }))
  datthr <- datthr[!is.na(datthr$bonferroni),]

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  color0 <- gg_color_hue(nlevels(dat0$color))
  tmp <- 21:25
  color0 <- gg_color_hue(nlevels(dat0$color))
  shape0 <- rep(tmp, ceiling(nlevels(dat0$color)/length(tmp)))[seq(nlevels(dat0$color))]
  names(color0) <- levels(dat0$color)

  if(plot){
    plot1 <- ggplot(dat0, aes(ecov, -log10(Pvalue))) +
          geom_point(aes(fill=color,shape=color), size=point.size) +
          theme_bw() + labs(x=xlab, y=ylab, fill=color_name, shape=color_name) +
          scale_x_discrete(expand=c(expand,expand)) +
          geom_hline(data=datthr, aes(yintercept=-log10(bonferroni)), color="blue",
                     linetype="dotted", linewidth=0.4) +
          facet_wrap(~trait, ncol=1, scales="free_y") +
          scale_fill_manual(values=color0) +
          scale_shape_manual(values=shape0) +
          theme(panel.grid.minor.x = element_blank(),
                panel.grid.major.x = element_blank(),
                axis.text.x = element_blank(),
                strip.text.x = element_text(size=10, margin=margin(t=2,b=2))
              ) +
          guides(fill=guide_legend(override.aes=list(size=4)))

  }else{
    plot1=NULL
  }
  return(list(out=out, plot=plot1))
}

#--------------------------------------
# data = dat
# xlab="ECOV"; ylab="-log10(Pvalue)"; color_name=NULL; point.size=1; expand=0.015;
# point.color="gray20"; rel_height=0.05; color_pal=NULL; text.size=3
plot_EWAS <- function(data, W=NULL, ecov_annotation=NULL, xlab="ECOV",
                      color_pal=NULL, alpha=0.5, group='trait',
                      ylab="-log10(P)", select=NULL, annotation_height=0.05,
                      color_name=NULL, point.color = "gray35", point.size=1,
                      text.size=3, expand.x=0.01,expand.y=0.01){

  data$group <- data[,group]
  # Select by group
  if(is.null(select)){
    select <- unique(data$group)
  }
  data <- data[data$group %in% select,]
  if(is.null(names(select))){
    data$group <- factor(as.character(data$group), levels=select)
  }else{
    data$group <- names(select)[match(data$group,select)]
    data$group <- factor(as.character(data$group), levels=names(select))
  }

  # Add stage info
  tmp <- get_ec_info(data$ecov)
  data$variable <- tmp$variable
  data$period <- tmp$period

  ecov_names <- as.character(unique(data$ecov))
  if(is.null(ecov_annotation)){
    data$color <- factor("Chromosome 1")
    color_name <- NULL

  }else{
    stopifnot(all(ecov_names %in% rownames(ecov_annotation)))
    if(!is.factor(ecov_annotation[,1])){
      ecov_annotation[,1] <- as.factor(ecov_annotation[,1])
    }
    ecov_levels <- levels(ecov_annotation[,1])
    data$color <- ecov_annotation[match(data$ecov, rownames(ecov_annotation)),1]
    data$color <- factor(as.character(data$color),
                         levels=ecov_levels[ecov_levels %in% unique(data$color)])
  }

  # get x position ordered by annotation
  tmp <- as.vector(unlist(lapply(split(data, data$color), function(x){
    as.vector(unlist(lapply(split(x,x$variable), function(z){
      tt <- z[match(levels(z$period), z$period),'ecov'] # Order by stage
      tt <- tt[!is.na(tt)]
      # tt <- ecov_names[ecov_names %in% unique(as.character(x$ecov))]
      if(is.null(W)){
        return(tt)
      }else{  # Make clusterization
        D <- SFSI::cov2dist(cor(W[,tt]))
        hc <- hclust(as.dist(D))
        return(tt[hc$order])
      }
    })))
  })))
  data$x <- match(data$ecov, tmp)
  dat0 <- data[order(data$group,data$x),]

  color_levels <- levels(dat0$color)
  datcol <- do.call(rbind,lapply(split(dat0, dat0$color), function(x){
    data.frame(color=x[1,"color"],
               xmin=min(x$x)-0.5, xmid=(min(x$x)+max(x$x))/2,
               xmax=max(x$x)+0.5, y=1)
  }))

  # Get a dataset for periods
  tmp <- split(dat0, dat0$group)
  datstage <- tmp[[which.max(unlist(lapply(tmp, nrow)))]]
  datstage$group <- factor(levels(dat0$group)[nlevels(dat0$group)], levels=levels(dat0$group))

  # Get threshold for manhattan plot
  datthr <- do.call(rbind,lapply(split(dat0, dat0$group), function(x){
    x[1,]
  }))
  datthr <- datthr[!is.na(datthr$bonferroni),]

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  # Color for 'CHR'
  if(is.null(color_pal)){
    color0 <- gg_color_hue(length(color_levels))
  }else{
    color0 <- rep(color_pal, ceiling(length(color_levels)/length(color_pal)))
    color0 <- color0[seq_along(color_levels)]
  }
  names(color0) <- color_levels

  # Add color for stages
  fun_color_range <- colorRampPalette(c("yellow", "red"))
  tmp <- fun_color_range(nlevels(dat0$period))
  names(tmp) <- levels(dat0$period)
  color0 <- c(color0, tmp)

  breaks0 <- c(datcol$xmin,datcol$xmax[nrow(datcol)])
  rg_stg <- range(-log10(datstage$Pvalue))

  # Data for enrichment
  datenr <- do.call(rbind,lapply(split(dat0,dat0$group),function(x){
    x$n0 <- nrow(x)
    x$n1 <- sum(-log10(x$Pvalue) > -log10(x$bonferroni))
    tt <- do.call(rbind,lapply(split(x,x$color),function(z){
      n0i <- sum(-log10(z$Pvalue) > -log10(z$bonferroni))
      if(n0i >0){
        Pr <- phyper(
                q = n0i-1,               # number of red marbles in the draw - 1 (see below)
                m = nrow(z),             # number of red marbles in urn
                n = z$n0[1] - nrow(z),   # number of green marbles in urn
                k = z$n1[1],                  # Number of drawn marbles
                lower.tail=FALSE)        # compute P( X > overlap ), hence the '-1' above
       }else{
         Pr <- NA
       }
       data.frame(z[1,c("group","color")],Pr=Pr,xmid=(min(z$x)+max(z$x))/2)
    }))
    tt
  }))
  datenr <- datenr[!is.na(datenr$Pr),]
  datenr$label <- significance_code(datenr$Pr)
  datenr$size <- ifelse(datenr$label=="ns",5,ifelse(datenr$label==".",6,3.5))
  datenr$vjust <- ifelse(datenr$label=="ns",2,ifelse(datenr$label==".",4,3))
  rownames(datenr) <- NULL

  plot1 <- ggplot(dat0) + theme_bw() +
          labs(x=NULL, y=ylab, fill=color_name) +
          #scale_x_discrete(expand=c(expand,expand)) +
          scale_x_continuous(breaks=datcol$xmid, labels=datcol$color,
                            expand=c(expand.x,expand.x)) +
          scale_y_continuous(expand=expansion(mult = c(expand.y, 0.1))) +
          facet_wrap(~group, ncol=1, scales="free_y")

  if(!is.null(ecov_annotation)){
    tmp <- -0.5*annotation_height*rg_stg[2] -1*expand.y*diff(rg_stg)
    plot1 <- plot1 +
             geom_rect(data=datcol, aes(xmin=xmin,ymin=0,xmax=xmax,ymax=Inf,fill=color),
                       alpha=alpha) +
             geom_tile(data=datstage, aes(x=x,y=tmp,fill=period),
                       color="gray20",height=annotation_height*diff(rg_stg)) +
             #geom_hline(yintercept=tmp, color="orange") +
             #geom_hline(yintercept=-1*expand.y*diff(rg_stg), color="blue") +
             #geom_hline(yintercept=-1*expand.y*rg_stg[2], color="green") +
             coord_cartesian(ylim=c(0,NA), clip="off")
  }

  plot1 <- plot1 +
        geom_point(aes(x=x, y=-log10(Pvalue)), color=point.color, size=point.size, shape=20) +
        geom_hline(data=datthr, aes(yintercept=-log10(bonferroni)), color="blue",
                   linetype="dashed")
  if(sum(datenr$label=="ns")>0){
    plot1 <- plot1 +
        geom_text(data=datenr[datenr$label=="ns",], aes(x=xmid,y=Inf,label=label),size=2.8,
                  color="red", vjust=1)
  }
  if(sum(datenr$label==".")>0){
    plot1 <- plot1 +
        geom_text(data=datenr[datenr$label==".",], aes(x=xmid,y=Inf,label=label),size=6,
                  color="red", vjust=0.3)
  }
  if(sum(!datenr$label%in%c("ns","."))>0){
    plot1 <- plot1 +
        geom_text(data=datenr[!datenr$label%in%c("ns","."),], aes(x=xmid,y=Inf,label=label),size=4.5,
                  color="red", vjust=1.1)
   }

  plot1 + scale_fill_manual(values=color0) +
          theme(panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.text.x = element_text(angle=90, vjust=0.5, hjust=1,
                                    margin=margin(t=105*(annotation_height+expand.y))),
            axis.ticks.x=element_blank(),
            strip.text.x = element_text(size=8, margin=margin(t=1.7,b=1.7)),
            legend.position="none"
          )
}

#--------------------------------------
# Plot enrichment
#--------------------------------------
# data <- out[[1]]
enrichment_plot <- function(data, traits=unique(data$trait),
                  xlab=NULL, ylab=NULL){

  library(ggrepel)
  stopifnot("TYPE" %in% toupper(colnames(data)))
  colnames(data)[grep("TYPE", toupper(colnames(data)))] <- "TYPE"

  dat <- do.call(rbind,lapply(split(data, data$trait),function(x){
    n0 <- sum(-log10(x$bonferroni) < -log10(x$Pvalue))
    tt <- do.call(rbind,lapply(split(x, x$TYPE),function(z){
      n1 <- sum(-log10(z$bonferroni) < -log10(z$Pvalue))
      data.frame( z[1,c("trait","TYPE")],
                  n0=nrow(z),
                  p0=nrow(z)/nrow(x),
                  n1=n1,
                  p1=ifelse(n0>0,n1/n0,0)
                )
    }))
    tt
  }))
  dat$enriched <- dat$p1 > dat$p0

  if(is.null(names(traits))){
    dat$trait <- factor(as.character(dat$trait), levels=traits)
  }else{
    dat$trait <- names(traits)[match(dat$trait,traits)]
    dat$trait <- factor(as.character(dat$trait), levels=names(traits))
  }

  if(is.null(xlab)) xlab <- "Proportion among all ECOVs"
  if(is.null(ylab)) ylab <- "Proportion in the discovery set"
  rg <- range(c(dat$p0,dat$p1))

  pt <- ggplot(dat, aes(p0, p1)) +
   geom_label_repel(aes(label=TYPE, color=enriched), max.overlaps=25,
                    size=2, label.padding=unit(0.18,"lines")) +
   geom_point(aes(fill=enriched), shape=21) +
   lims(x=rg, y=rg) + theme_bw() + #labs(x=NULL,y=NULL) +
   facet_wrap(~trait) +
   labs(x="Proportion among all ECOVs",y="Proportion in the discovery set") +
   geom_abline(slope=1,intercept=0,color="gray50",linetype="dashed") +
   theme(legend.position="none",
         strip.text.x = element_text(size=10, margin=margin(t=2.5,b=2.5))
       )
  pt

}

#--------------------------------------
# W=NULL; ecov_annotation=NULL; xlab="ECOV"; color_pal=NULL; alpha=0.5; group="region"
# ylab="Variance"; select=NULL; hline.color="blue"; vline.color = "gray40";
# text.size=3; expand.x=0.01; expand.y=0.01;annotation_height=0.05
# bar.color=NA; comp_mean=NULL; bar.width=1
plot_varcomp <- function(data, comp_names, W=NULL, ecov_annotation=NULL,
                      xlab="ECOV", color_pal=NULL, alpha=0.5, group="region",
                      comp_mean=NULL,
                      ylab="Variance", select=NULL, annotation_height=0.05,
                      bar.width=1, bar.color=NA, comp_color=NULL,
                      hline.color="blue", vline.color = "gray40",
                      text.size=3, expand.x=0.01, expand.y=0.0){
  # head(data)
  data$group <- data[,group]
  # Select by group
  if(is.null(select)){
    select <- unique(data$group)
  }
  data <- data[data$group %in% select,]
  if(is.null(names(select))){
    data$group <- factor(as.character(data$group), levels=select)
  }else{
    data$group <- names(select)[match(data$group,select)]
    data$group <- factor(as.character(data$group), levels=names(select))
  }

  # Add stage info
  tmp <- get_ec_info(data$ecov)
  data$variable <- tmp$variable
  data$period <- tmp$period

  # Get percentage variance
  stopifnot(all(comp_names %in% colnames(data)))
  if(!is.null(names(comp_names))){
    if(length(unique(names(comp_names))) < length(comp_names)){
      stop("Names of 'comp_names' should be all different from each other")
    }
    colnames(data)[match(comp_names,colnames(data))] <- names(comp_names)
    comp_names <- names(comp_names)
  }
  dat0 <- do.call(rbind,lapply(split(data, data$group),function(x){
      do.call(rbind,lapply(split(x, x$ecov),function(z){
        if(nrow(z)>1){
          stop("More than one observation for ecov '",z[1,'ecov'],"'")
        }
        VAR <- as.vector(as.matrix(z[,comp_names]))
        names(VAR) <- comp_names
        VAR <- VAR[!is.na(VAR)]
        data.frame(z[rep(1,length(VAR)),c("ecov","group","variable","period")],
                   comp=names(VAR),VAR=as.vector(VAR),
                   propVAR=as.vector(VAR/sum(VAR)))
      }))
  }))
  rownames(dat0) <- NULL
  dat0$comp <- factor(as.character(dat0$comp), levels=comp_names)

  ecov_names <- as.character(unique(data$ecov))
  if(is.null(ecov_annotation)){
    dat0$color <- factor("Chromosome 1")

  }else{
    stopifnot(all(ecov_names %in% rownames(ecov_annotation)))
    if(!is.factor(ecov_annotation[,1])){
      ecov_annotation[,1] <- as.factor(ecov_annotation[,1])
    }
    ecov_levels <- levels(ecov_annotation[,1])
    dat0$color <- ecov_annotation[match(dat0$ecov, rownames(ecov_annotation)),1]
    dat0$color <- factor(as.character(dat0$color),
                         levels=ecov_levels[ecov_levels %in% unique(dat0$color)])
  }

  # get x position ordered by annotation
  tmp <- as.vector(unlist(lapply(split(dat0, dat0$color), function(x){
    as.vector(unlist(lapply(split(x,x$variable), function(z){
      tt <- z[match(levels(z$period), z$period),'ecov'] # Order by stage
      tt <- tt[!is.na(tt)]
      # tt <- ecov_names[ecov_names %in% unique(as.character(x$ecov))]
      if(is.null(W)){
        return(tt)
      }else{  # Make clusterization
        D <- SFSI::cov2dist(cor(W[,tt]))
        hc <- hclust(as.dist(D))
        return(tt[hc$order])
      }
    })))
  })))
  dat0$x <- match(dat0$ecov, tmp)
  dat0 <- dat0[order(dat0$group,dat0$x),]

  color_levels <- levels(dat0$color)
  datcol <- do.call(rbind,lapply(split(dat0, dat0$color), function(x){
    data.frame(color=x[1,"color"],
               xmin=min(x$x)-0.5, xmid=(min(x$x)+max(x$x))/2,
               xmax=max(x$x)+0.5, y=1)
  }))
  datcol$group <- factor(levels(dat0$group)[nlevels(dat0$group)], levels=levels(dat0$group))

  # Get a dataset for periods
  datstage <- split(dat0, dat0$group)[[levels(dat0$group)[nlevels(dat0$group)]]]
  tmp <- split(datstage, datstage$comp)
  datstage <- tmp[[which.max(unlist(lapply(tmp, nrow)))]]

  fun_color_range <- colorRampPalette(c("yellow", "red"))

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if(is.null(color_pal)){
    color0 <- fun_color_range(nlevels(dat0$period))
    #color0 <- gg_color_hue(length(color_levels))
    names(color0) <- levels(dat0$period)
  }else{
    tmp <- rep(color_pal, ceiling(length(color_levels)/length(color_pal)))
    color0 <- tmp[seq_along(color_levels)]
    names(color0) <- c(color_levels)
  }
  # Add color for variance components
  if(is.null(comp_color)){
    tmp <- RColorBrewer::brewer.pal(8, "Pastel1")
    tmp <- tmp[seq_along(comp_names)]
  }else{
    if(length(comp_color)<length(comp_names)){
      stop(length(comp_color),' colors were provided for variance components,',
           ' at least ',length(comp_names),' are needed')
    }
    tmp <- comp_color[seq_along(comp_names)]
  }
  names(tmp) <- c(comp_names)
  color0 <- c(color0, tmp)

  datcol$color2 <- color0[datcol$color]

  breaks0 <- c(datcol$xmin,datcol$xmax[nrow(datcol)])

  # Average values
  if(is.null(comp_mean)){
    comp_mean <- lapply(split(dat0$comp,dat0$group),function(x)comp_names[comp_names%in%x][1])
    message("Average is shown for the first component level within group")
  }
  stopifnot(all(levels(dat0$group) %in% names(comp_mean)))

  datmean <- do.call(rbind,lapply(split(dat0, dat0$group),function(x){
    tt <- comp_mean[[x$group[1]]]
    mm <- tapply(X=x$propVAR, INDEX=as.character(x$comp), FUN=mean)
    index <- match(tt, names(mm))
    if(length(index)==0){
      stop("Component(s) ",paste(tt,collapse=",")," were not found in group '",x$group[1],"'")
    }
    data.frame(x[1,"group",drop=F], propVAR=sum(mm[index]))
  }))

  dat0$comp <- factor(as.character(dat0$comp), levels=rev(levels(dat0$comp)))

  plot1 <- ggplot(dat0) + theme_bw() +
          geom_col(aes(x=x, y=propVAR, fill=comp), color=bar.color, width=bar.width) +
          labs(x=NULL, y=ylab, fill=NULL) +
          scale_x_continuous(breaks=datcol$xmid, labels=datcol$color,
                            expand=c(expand.x,expand.x)) +
          scale_y_continuous(expand=c(expand.y,expand.y)) +
          facet_wrap(~group, ncol=1, scales="free_y")

  if(!is.null(ecov_annotation)){
    if(is.null(color_pal)){
      plot1 <- plot1 +
               geom_tile(data=datstage, aes(x=x,y=-0.5*(annotation_height+expand.y),fill=period),
                         color="gray20",height=annotation_height)
    }else{
      plot1 <- plot1 +
               geom_rect(data=datcol, aes(xmin=xmin,ymin=-1*annotation_height,xmax=xmax,ymax=0,fill=color))
    }
    plot1 <- plot1 + coord_cartesian(ylim=c(0,1), clip="off") +
                     geom_vline(xintercept=breaks0, color=vline.color)
  }

  plot1 <- plot1 +
           geom_hline(data=datmean, aes(yintercept=propVAR),
                      linetype="dashed", color=hline.color) +
           scale_fill_manual(values=color0, breaks=comp_names) +
           theme(panel.grid.minor.x = element_blank(),
               panel.grid.major.x = element_blank(),
               axis.text.x = element_text(angle=90, vjust=0.5, hjust=1,
                                          margin=margin(t=110*(annotation_height+expand.y))),
               strip.text.x = element_text(size=10, margin=margin(t=2,b=2)),
               legend.position="bottom", legend.justification="right",
               axis.ticks.x=element_blank(),
               legend.margin=margin(0,0,0,0),
               legend.box.margin=margin(t=-10)
             )
  plot1

}


plot_EWAS_old <- function(data, W=NULL, ecov_annotation=NULL, xlab="ECOV",
                      color_pal=NULL, alpha=0.5,
                      ylab="-log10(Pvalue)", traits=unique(data$trait),
                      color_name=NULL, point.color = "gray20", point.size=1,
                      text.size=3, expand=0.015,
                      rel_height=0.05){

  ecov_names <- as.character(unique(data$ecov))
  if(is.null(ecov_annotation)){
    data$color <- factor("Chromosome 1")
    color_name <- NULL

  }else{
    stopifnot(all(ecov_names %in% rownames(ecov_annotation)))
    if(!is.factor(ecov_annotation[,1])){
      ecov_annotation[,1] <- as.factor(ecov_annotation[,1])
    }
    data$color <- ecov_annotation[match(data$ecov, rownames(ecov_annotation)),1]
    color_name <- colnames(ecov_annotation)[1]
    if(!is.factor(ecov_annotation[,1])){
      ecov_annotation[,1] <- as.factor(ecov_annotation[,1])
    }
  }

  dat0 <- data[data$trait %in% traits,]
  if(is.null(names(traits))){
    dat0$trait <- factor(as.character(dat0$trait), levels=traits)
  }else{
    dat0$trait <- names(traits)[match(dat0$trait,traits)]
    dat0$trait <- factor(as.character(dat0$trait), levels=names(traits))
  }

  # get x position ordered by annotation
  tmp <- as.vector(unlist(lapply(split(dat0, dat0$color), function(x){
    tt <- ecov_names[ecov_names %in% unique(as.character(x$ecov))]
    if(is.null(W)){
      return(tt)
    }else{  # Make clusterization
      D <- SFSI::cov2dist(cor(W[,tt]))
      hc <- hclust(as.dist(D))
      return(tt[hc$order])
    }
  })))
  dat0$x <- match(dat0$ecov, tmp)
  dat0 <- dat0[order(dat0$trait,dat0$x),]

  color_levels <- levels(dat0$color)
  datcol <- do.call(rbind,lapply(split(dat0, dat0$color), function(x){
    #tt <- dat0[dat0$color==color_levels[k],]
    #data.frame(min=min(tt$x)-0.5, mid=(min(tt$x)+max(tt$x))/2,max=max(tt$x)+0.5)
    data.frame(color=x[1,"color"],
               xmin=min(x$x)-0.5, xmid=(min(x$x)+max(x$x))/2,
               xmax=max(x$x)+0.5, y=1)
  }))
  dattext <- datcol[]
  dattext$trait <- factor(levels(dat0$trait)[nlevels(dat0$trait)], levels=levels(dat0$trait))

  # Get threshold for manhattan plot
  datthr <- do.call(rbind,lapply(split(dat0, dat0$trait), function(x){
    x[1,]
  }))
  datthr <- datthr[!is.na(datthr$bonferroni),]

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if(is.null(color_pal)){
    color0 <- gg_color_hue(length(color_levels))
  }else{
    color0 <- rep(color_pal, ceiling(length(color_levels)/length(color_pal)))
    color0 <- color0[seq_along(color_levels)]
  }
  names(color0) <- color_levels

  breaks0 <- c(datcol$xmin,datcol$xmax[nrow(datcol)])

  plot1 <- ggplot(dat0) + theme_bw() +
          labs(x=NULL, y=ylab, fill=color_name) +
          #scale_x_discrete(expand=c(expand,expand)) +
          scale_x_continuous(breaks=datcol$xmid, labels=datcol$color,
                            expand=c(expand,expand)) +
          facet_wrap(~trait, ncol=1, scales="free_y")

  if(!is.null(ecov_annotation)){
    plot1 <- plot1 +
             geom_rect(data=datcol, aes(xmin=xmin,ymin=-Inf,xmax=xmax,ymax=Inf,fill=color),
                       alpha=0.5)
  }

  plot1 <- plot1 +
        geom_point(aes(x=x, y=-log10(Pvalue)), color=point.color, size=point.size, shape=20) +
        geom_hline(data=datthr, aes(yintercept=-log10(bonferroni)), color="blue",
                   linetype="dashed") +
        scale_fill_manual(values=color0) +
         theme(panel.grid.minor.x = element_blank(),
               panel.grid.major.x = element_blank(),
               axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
               strip.text.x = element_text(size=10, margin=margin(t=2,b=2)),
               legend.position="none"
             )

  # Plot spectrum bar
  plotbar <- ggplot(datcol, aes(y=y)) +
    geom_rect(aes(xmin=xmin,ymin=-Inf,xmax=xmax,ymax=Inf,fill=color),
              alpha=0.5) +
    scale_x_continuous(breaks=NULL, expand=c(expand,expand),
               sec.axis = sec_axis(~., breaks=breaks0)) +
    scale_y_continuous(breaks=c(1.0),labels=c(color_name),expand=c(0,0), limits=c(0.8, 1.2)) +
    scale_fill_manual(values=color0) +
    geom_text(aes(x=xmid, label=color), size=text.size, angle=90) +
    labs(x=NULL,y=NULL) +
    theme(panel.background=element_rect(color="gray100",fill="gray100"),
          axis.text.x=element_blank(),
          axis.text.y=element_text(face="bold"),
          #element_text(size=5.5,angle=90, vjust=0.5, hjust=1)
          #axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin=unit(c(1,5,1,5), "pt"),
          legend.position="none")

    # arrange into grid and align
    #cowplot::plot_grid(plot1, plotbar, align='v', axis='lr', ncol=1,
    #         rel_heights=c(1,rel_height))
    plot1
}

#--------------------------------------

cartToPolar=function(x,y){
  r=sqrt(x^2+y^2)
  theta=atan(y/x)
  if((x<0 & y>0) | (x<0 & y<=0)) theta <- theta+pi
  if(x>0 & y<0) theta <- theta+2*pi

  return(c('r'=r,'theta'=theta))
}

polarToCart=function(r,theta){
  y=r*sin(theta)
  x=r*cos(theta)
  return(c('x'=x,'y'=y))
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}

# rotate_point(c(1,1),10)
rotate_point <- function(xy, rotate=45){
    tmp <- cartToPolar(xy[1],xy[2])
    r <- tmp['r']
    tmp['theta'] <- tmp['theta'] + rotate*pi/180

    xy2 <- polarToCart(tmp['r'],tmp['theta'])

    dat <- circleFun(c(0,0), 2*r, npoints=100)

    ggplot2::ggplot(dat,ggplot2::aes(x,y)) + ggplot2::theme_classic() +
    ggplot2::geom_path() +
    ggplot2::geom_segment(x=0, y=0, xend=xy[1], yend=xy[2], color="gray60") +
    ggplot2::geom_segment(x=0, y=0, xend=xy2[1], yend=xy2[2], color="blue") +
    ggplot2::geom_point(x=0,y=0, size=1, color="black")
}

#------------------------------------
# Rotate xy a given angle in degrees
#------------------------------------
rotate_xy <- function(xy, angle){

  if(length(dim(xy)) != 2L){
    if(length(xy)==2){
      xy <- matrix(xy, ncol=2)
    }else{
      stop("'xy' should be a matrix with 2 columns")
    }
  }
  n <- nrow(xy)
  tmp <- do.call(rbind,lapply(1:n, function(i) as.vector(cartToPolar(xy[i,1],xy[i,2]))))

  res <- vector("list",length(angle))
  for(k in 1:length(angle)){
    theta <- tmp[,2] + angle[k]*pi/180  # add to theta
    newxy <- do.call(rbind,lapply(1:n,function(i)
              polarToCart(as.vector(tmp[i,1]),theta[i])))
    if(n==1){
      res[[k]] <- drop(newxy)
    }else{
      res[[k]] <- newxy
    }
  }
  if(length(angle)==1) res <- res[[1]]
  res
}

get_rotated_axis <- function(r=1, side=1, angle){
  stopifnot(side %in% c(1,2))
  if(side == 1L){
    xx <- matrix(c(-r,0,r,0),ncol=2,byrow=TRUE)
  }else{
    xx <- matrix(c(0,-r,0,r),ncol=2,byrow=TRUE)
  }
  # rotate_point(xx[1,],angle)
  do.call(rbind,lapply(1:nrow(xx),function(i){
    ax <- rotate_xy(xx[i,], angle)
    tmp <- as.vector(cartToPolar(ax[1],ax[2]))
    if(side==1){
      x2 <- (r^2)/tan(tmp[2])^2
    }else{
      x2 <- (r^2)*tan(tmp[2])^2
    }
    tmp[1] <- sqrt(x2 + r^2)
    polarToCart(tmp[1],tmp[2])
  }))
}
#------------------------------------

rotate_PC <- function(PC, latitude, longitude, signPC=NULL){

  stopifnot(length(latitude) == nrow(PC))
  stopifnot(length(longitude) == nrow(PC))

  # Make a grid search for rotation
  grid = seq(0, 360, length=1000)

  stopifnot(ncol(PC)>1L)
  PC <- PC[,1:2]

  if(is.null(signPC)){
    M <- rbind(c(1,1), c(-1,-1), c(-1,1), c(1,-1))
  }else{
    M <- t(signPC)
  }

  B <- matrix(NA, nrow=length(grid), ncol=nrow(M))
  for(k in 1:nrow(M)){

    PC0 <- sweep(PC, 2, M[k,],FUN='*')

    b <- unlist(lapply(seq_along(grid), function(i){
      xy <- rotate_xy(xy=PC0, angle=grid[i])
      #c(coef(lm(xy[,1]~longitude))[2], coef(lm(xy[,2]~latitude))[2])
      #c(abs(cor(longitude,xy[,1])), abs(cor(latitude,xy[,2])))
      cor(longitude,xy[,1]) + cor(latitude,xy[,2])
    }))

    B[,k] <- b
  }

  # Select the PC signs with highest regression coef
  index <- which.max(apply(B,2,max))
  if(is.null(signPC)){
    if(index <= 2){ # cases (1,1) and (-1,-1)
      index <- ifelse(grid[which.max(B[,1])] < grid[which.max(B[,2])],1,2)
    }else{     # cases (1,-1) and (-1,1)
      index <- ifelse(grid[which.max(B[,3])] < grid[which.max(B[,4])],3,4)
    }
  }

  signPC <- as.vector(M[index,])
  PC <- sweep(PC, 2, signPC, FUN='*') # Change appropiate axis signs

  # Center to range [-1, +1]
  moveAxis <- function(A){
    mm <- apply(A, 2, function(x) (max(x)+min(x))/2)
    A <- sweep(A, 2, mm, FUN="-")
    apply(A, 2, function(x)x/max(abs(x)))
  }
  # plot(grid, rowSums(B[[index]]))

  shift <- grid[which.max(B[,index])]

  xy <- rotate_xy(PC, shift)  # Get the rotated PC's
  PC <- moveAxis(PC)
  xy <- moveAxis(xy)

  # Rotated axis corners: x: (-1,0) & (1,0).  y: (0,1) & (0,-1)
  # rotate_point(c(-1,0),88)
  x.axis <- get_rotated_axis(r=1, side=1, shift)
  y.axis <- get_rotated_axis(r=1, side=2, shift)

  list(PC=PC, xy=xy, signPC=signPC, shift=shift, x.axis=x.axis, y.axis=y.axis)
}
