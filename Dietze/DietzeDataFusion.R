#---------------- Load plot and tree ring data. -------------------------------------------------------#
setwd("C:/Users/eth_admin/Documents/GitHub/DBH/Dietze")
rm(list=ls())


## Read tree data
trees <- read.csv("treecores2014.csv")

##' @title Clean_Tucson
##' @name Clean_Tucson
##' @description  tree core QAQC
Clean_Tucson <- function(file){
  lines = scan(file,character(),sep="\n")
  split = strsplit(lines," ")
  tags   = NULL
  decade = NULL
  for(i in 1:length(split)){
    tags[i]   = split[[i]][1]
    decade[i] = split[[i]][2]
  }
  utags = unique(tags)
  newfile = paste(file,".COR.txt",sep="")
  if(file.exists(newfile)) file.remove(newfile)
  for(tag in utags){
    rows = rev(which(tags == tag))
    keep = 1
    for(i in 1:length(rows)){
      if(rows[i] - rows[keep] >= -1){
        keep = i
      } else {
        break
      }
    }
    keep = min(keep,length(rows))
    rows = rev(rows[1:keep])
    append = file.exists(newfile)
    write(lines[rows],newfile,append=append)
  }
  return(newfile)
}

##' @title Read_Tucson
##' @name Read_Tucson
##' @export
##' @description wrapper around read.tucson that loads a whole directory of tree ring files 
##' and calls a 'clean' function that removes redundant records 
##' (WinDendro can sometimes create duplicate records when editing)
Read_Tucson <- function(folder){
  
  require(dplR)
  
  filenames <- dir(folder,pattern =  "TXT",full.names=TRUE)
  filenames <- c(filenames,dir(folder,pattern =  "rwl",full.names=TRUE))
  filenames <- c(filenames,dir(folder,pattern = "rw",full.names=TRUE))
  corrected = grep(pattern="COR.txt",x=filenames)
  if(length(corrected)>0){
    filenames = filenames[-corrected]
  }
  filedata <- list()
  for (file in filenames){
    file = Clean_Tucson(file)
    filedata[[file]] <- read.tucson(file, header = FALSE)
  }
  
  return(filedata)
}


## Read tree ring data
rings <- Read_Tucson("C:/Users/eth_admin/Documents/GitHub/DBH/Dietze/Tucson")

###### helper function matching up two data types
matchInventoryRings <- function(trees,rings,extractor="TreeCode",nyears=30,coredOnly=TRUE){
  
  ## build tree codes
  id.build = function(x){do.call(paste0("to.",extractor),x)}
  names(trees) = toupper(names(trees))
  tree.ID = id.build(list(SITE=trees$SITE,PLOT=trees$PLOT,SUB=trees$SUB,TAG=trees$TAG))
  
  ## build tree ring codes
  if(is.list(rings)){
    ring.file <- rep(names(rings),times=sapply(rings,ncol))
    rings <- combine.rwl(rings)
  }
  ring.ID <- names(rings)
  id.extract = function(x){do.call(paste0("from.",extractor),list(x=x))}
  ring.info <- id.extract(ring.ID)
  
  ## matching up data sets by tree
  mch = match(tree.ID,ring.ID)
  cored = apply(!is.na(trees[,grep("DATE_CORE_COLLECT",names(trees))]),1,any)
  unmatched = which(cored & is.na(mch))
  write.table(tree.ID[unmatched],file="unmatched.txt")
  mch[duplicated(mch)] <- NA  ## if there's multiple stems, match the first
  
  ## combine data into one table
  combined = cbind(trees,t(as.matrix(rings))[mch,-(nyears-1):0 + nrow(rings)])
  if(coredOnly==TRUE){
    combined = combined[!is.na(combined$"2000"),]
  }
  return(combined)
}


## Match observations & format for JAGS
combined <- matchInventoryRings(trees,rings,extractor="Tag",nyears=36,coredOnly=FALSE) #WARNINGS

#### helper function formatting data for JAGS
##' @title Format ring & plot data for JAGA
##' @name  buildJAGSdata_InventoryRings
##' @param combined  object returned from matchInventoryRings. Matrix with both increment and plot data
##' @param inc.unit.conv  conversion factor from loaded increments to cm (of radius)
##' @return list
##' @author Michael Dietze
##' @description builds the JAGS data object for the tree ring / inventory fusion code
##' also sets all the priors
##' @export
buildJAGSdata_InventoryRings <- function(combined,inc.unit.conv = 0.1){
  
  ## pull out growth to a matrix, convert to cm of diameter
  y = as.matrix(combined[,!is.na(as.numeric(colnames(combined)))])*inc.unit.conv*2
  time = as.numeric(colnames(y))
  
  ## pull out diameter to a matrix
  DBH.cols = grep("DBH",colnames(combined))
  DBH = as.matrix(combined[,DBH.cols])
  class(DBH) <- "numeric"
  z = matrix(NA,nrow(y),ncol(y))
  DBH.years = as.numeric(sub("DBH","",colnames(combined)[DBH.cols]))
  DBH.years = ifelse(DBH.years < 20,DBH.years+2000,DBH.years+1900)
  z[,which(time %in% DBH.years)] = DBH
  
  ## if present, pull out mortality and recruitment
  COND.cols = grep("COND",colnames(combined))
  if(length(COND.cols)>0){
    COND = as.matrix(combined[,COND.cols])
    w = matrix(NA,nrow(y),ncol(y))
    COND.years = as.numeric(sub("COND","",colnames(combined)[COND.cols]))
    COND.years = ifelse(COND.years < 20,COND.years+2000,ifelse(COND.years < 100,COND.years+1900,COND.years))
    w[,which(time %in% COND.years)] = COND
    ## convert COND matrix to numeric 0/1
    w[w=="L"] = 1
    w[w=='miss'] = NA
    for(i in 1:nrow(w)){
      ## recruitment
      r = which(w[i,] %in% 'R')
      if(length(r)>0){
        w[i,seq_len(r-1)] = 0  ## should really set last census, not last year, to 0, and then put NAs in between  ****
        w[i,r] = 1
      }
      ## mortality
      m = which(w[i,] %in% 'M')
      if(length(m)>0){
        w[i,m:ncol(w)] = 0
      }
      ## known to be alive
      l = which(w[i,] %in% '1')
      if(length(l)>1){
        
      }
    }
    class(w) <- "numeric"
    ## Fill in COND matrix
    
    
  } else {
    w = matrix(1,nrow(y),ncol(y))
  }
  
  ## build data object for JAGS
  n = nrow(y)
  data <- list(y=y[1:n,],z = z[1:n,],ni=n,nt=ncol(y),x_ic=1,tau_ic=0.0001,
               a_dbh=16,r_dbh=8,a_inc=0.001,r_inc=1,a_add=1,r_add=1,time=time)
  
  return(data)
  
}


data <- buildJAGSdata_InventoryRings(combined) #WARNINGS
status.end()

#---------------- Load plot and tree ring data. -------------------------------------------------------#
status.start("TREE RING MODEL")
## Tree Ring model
n.iter = 3000

##### DBH + tree-ring data fusion
##### JAGS model
##' @name InventoryGrowthFusion
##' @title InventoryGrowthFusion
##' @description this code fuses forest inventory data with tree growth data (tree ring or dendrometer band)
##' for the same plots. Code is a rewrite of Clark et al 2007 Ecol Appl into JAGS
##' 
##' @param data  list of data inputs
##' @param random = whether or not to include random effects
##' @note Requires JAGS
##' @return an mcmc.list object
##' @export
InventoryGrowthFusion <- function(data,n.iter,random=TRUE){
  require(rjags)
  
  TreeDataFusionMV = "
  model{
  ### Loop over all individuals
  for(i in 1:ni){
  
  #### Data Model: DBH
  for(t in 1:nt){
  z[i,t] ~ dnorm(x[i,t],tau_dbh)
  }
  
  #### Data Model: growth
  for(t in 2:nt){
  inc[i,t] <- x[i,t]-x[i,t-1]
  y[i,t] ~ dnorm(inc[i,t],tau_inc)
  }
  
  #### Process Model
  for(t in 2:nt){
  Dnew[i,t] <- x[i,t-1] + mu
  x[i,t]~dnorm(Dnew[i,t],tau_add)
  }
  
  x[i,1] ~ dnorm(x_ic,tau_ic)
  }  ## end loop over individuals
  
  #### Priors
  tau_dbh ~ dgamma(a_dbh,r_dbh)
  tau_inc ~ dgamma(a_inc,r_inc)
  tau_add ~ dgamma(a_add,r_add)
  mu ~ dnorm(0.5,0.5)
  }"

  if(random==TRUE){
    ## version with tree and year random effects
    TreeDataFusionMV = "
    model{
    
    ### Loop over all individuals
    for(i in 1:ni){
    
    #### Data Model: DBH
    for(t in 1:nt){
    z[i,t] ~ dnorm(x[i,t],tau_dbh)
    }
    
    #### Data Model: growth
    for(t in 2:nt){
    inc[i,t] <- x[i,t]-x[i,t-1]
    y[i,t] ~ dnorm(inc[i,t],tau_inc)
    }
    
    #### Process Model
    for(t in 2:nt){
    Dnew[i,t] <- x[i,t-1] + mu + ind[i] + year[t]
    x[i,t]~dnorm(Dnew[i,t],tau_add)
    }
    
    ## individual effects
    ind[i] ~ dnorm(0,tau_ind)
    
    ## initial condition
    x[i,1] ~ dnorm(x_ic,tau_ic)
    }  ## end loop over individuals
    
    ## year effects
    for(t in 1:nt){
    year[t] ~ dnorm(0,tau_yr)
    }
    
    
    #### Priors
    tau_dbh ~ dgamma(a_dbh,r_dbh)
    tau_inc ~ dgamma(a_inc,r_inc)
    tau_add ~ dgamma(a_add,r_add)
    tau_ind ~ dgamma(1,0.1)
    tau_yr  ~ dgamma(1,0.1)
    mu ~ dnorm(0.5,0.5)
    
    }"
}
    
  ## state variable initial condition
  z0 = t(apply(data$y,1,function(y){-rev(cumsum(rev(y)))})) + data$z[,ncol(data$z)] 
  
  ## JAGS initial conditions
  nchain = 3
  init <- list()
  for(i in 1:nchain){
    y.samp = sample(data$y,length(data$y),replace=TRUE)
    init[[i]] <- list(x = z0,tau_add=runif(1,1,5)/var(diff(y.samp),na.rm=TRUE),
                      tau_dbh=1,tau_inc=1500,tau_ind=50,tau_yr=100,ind=rep(0,data$ni),year=rep(0,data$nt))
  }
  
  ## compile JAGS model
  j.model   <- jags.model (file = textConnection(TreeDataFusionMV),
                           data = data,
                           inits = init,
                           n.chains = 3)
  ## burn-in
  jags.out   <- coda.samples (model = j.model,
                              variable.names = c("tau_add","tau_dbh","tau_inc","mu","tau_ind","tau_yr"),
                              n.iter = min(n.iter,2000))
  plot(jags.out)
  
  ## run MCMC
  jags.out   <- coda.samples (model = j.model,
                              variable.names = c("x","tau_add","tau_dbh","tau_inc","mu",
                                                 "tau_ind","tau_yr","ind","year"),
                              n.iter = n.iter)
  
  return(jags.out)



jags.out = InventoryGrowthFusion(data,n.iter=n.iter) #WARNINGS
save(trees,rings,combined,data,jags.out,
     file=file.path(settings$outdir,"treering.Rdata"))

##### diagnostics
##' @name InventoryGrowthFusionDiagnostics
##' @title InventoryGrowthFusionDiagnostics
##' @param jags.out output mcmc.list from InventoryGrowthFusion
##' @param combined  data output from matchInventoryRings
##' @author Michael Dietze
##' @export 
InventoryGrowthFusionDiagnostics <- function(jags.out,combined){
  
  #### Diagnostic plots
  
  ### DBH
  #par(mfrow=c(3,2))
  layout(matrix(1:8,4,2,byrow=TRUE))
  out <- as.matrix(jags.out)
  x.cols = which(substr(colnames(out),1,1)=="x")
  ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))
  ci.names = parse.MatrixNames(colnames(ci),numeric=TRUE)
  
  smp = sample.int(data$ni,min(8,data$ni))
  for(i in smp){
    sel = which(ci.names$row == i)
    rng = c(range(ci[,sel],na.rm=TRUE),range(data$z[i,],na.rm=TRUE))
    plot(data$time,ci[2,sel],type='n',ylim=range(rng),ylab="DBH (cm)",main=i)
    ciEnvelope(data$time,ci[1,sel],ci[3,sel],col="lightBlue")
    points(data$time,data$z[i,],pch="+",cex=1.5)
    #   lines(data$time,z0[i,],lty=2)
    
    ## growth
    sel = which(ci.names$row == i)
    inc.mcmc = apply(out[,x.cols[sel]],1,diff)
    inc.ci = apply(inc.mcmc,1,quantile,c(0.025,0.5,0.975))*5
    #inc.names = parse.MatrixNames(colnames(ci),numeric=TRUE)
    
    plot(data$time[-1],inc.ci[2,],type='n',ylim=range(inc.ci,na.rm=TRUE),ylab="Ring Increment (mm)")
    ciEnvelope(data$time[-1],inc.ci[1,],inc.ci[3,],col="lightBlue")
    points(data$time,data$y[i,]*5,pch="+",cex=1.5,type='b',lty=2)
  }
  
  if(FALSE){
    ##check a DBH
    plot(out[,which(colnames(out)=="x[3,31]")])
    abline(h=z[3,31],col=2,lwd=2)
    hist(out[,which(colnames(out)=="x[3,31]")])
    abline(v=z[3,31],col=2,lwd=2)
  }
  
  ## process model
  vars = (1:ncol(out))[-c(which(substr(colnames(out),1,1)=="x"),grep("tau",colnames(out)),
                          grep("year",colnames(out)),grep("ind",colnames(out)))]
  par(mfrow=c(1,1))
  for(i in vars){
    hist(out[,i],main=colnames(out)[i])
  }
  if(length(vars)>1) pairs(out[,vars])
  
  ## Standard Deviations
  #layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
  par(mfrow=c(2,3))
  prec = out[,grep("tau",colnames(out))]
  for(i in 1:ncol(prec)){
    hist(1/sqrt(prec[,i]),main=colnames(prec)[i])
  }
  cor(prec)
  #  pairs(prec)
  
  
  par(mfrow=c(1,1))
  ### YEAR
  year.cols = grep("year",colnames(out))
  if(length(year.cols>0)){
    ci.yr <- apply(out[,year.cols],2,quantile,c(0.025,0.5,0.975))
    plot(data$time,ci.yr[2,],type='n',ylim=range(ci.yr,na.rm=TRUE),ylab="Year Effect")
    ciEnvelope(data$time,ci.yr[1,],ci.yr[3,],col="lightBlue")
    lines(data$time,ci.yr[2,],lty=1,lwd=2)
    abline(h=0,lty=2)
  }
  
  ### INDIV
  ind.cols= which(substr(colnames(out),1,3)=="ind")
  if(length(ind.cols)>0){
    boxplot(out[,ind.cols],horizontal=TRUE,outline=FALSE,col=as.factor(combined$PLOT))
    abline(v=0,lty=2)
    tapply(apply(out[,ind.cols],2,mean),combined$PLOT,mean)
    table(combined$PLOT)
    
    spp = combined$SPP
    #    boxplot(out[order(spp),ind.cols],horizontal=TRUE,outline=FALSE,col=spp[order(spp)])
    boxplot(out[,ind.cols],horizontal=TRUE,outline=FALSE,col=spp)
    abline(v=0,lty=2)
    spp.code = levels(spp)[table(spp)>0]
    legend("bottomright",legend=rev(spp.code),col=rev(which(table(spp)>0)),lwd=4)
    tapply(apply(out[,ind.cols],2,mean),combined$SPP,mean)
  }
}
########### NEXT STEPS: ############
#what explain's the year effects? climate
#what explain's the individual effects? size, species, canopy position, plot -> landscape

pdf(file.path(settings$outdir,"treering.Diagnostics.pdf"))
InventoryGrowthFusionDiagnostics(jags.out,combined)
dev.off()
status.end()