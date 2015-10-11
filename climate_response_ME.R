###############################################################
## Climate sensitivity of TRW in Monument Canyon (Jemez Mts)
###############################################################

### modified from Flurin Babst
### Sept 2015

### load packages and functions

library(dplR)
library(bootRes)
library(pspline)

monthly.correlations <- function(treeringdata,monthlydata,startyear,endyear) {
temp<- ts.union(treeringdata,lag(monthlydata,-1),monthlydata)
colnames(temp) <- c("treeringdata",paste("p",colnames(monthlydata),sep=""),colnames(monthlydata))
temp <- window(temp,startyear,endyear)
seasons <- cbind(rowMeans(temp[,2:13]),rowMeans(temp[,7:9]),rowMeans(temp[,5:10]), rowMeans(temp[,13:17]), rowMeans(temp[,10:11]),rowMeans(temp[,13:15]),rowMeans(temp[,14:16]),rowMeans(temp[,15:17]),rowMeans(temp[,16:17]),rowMeans(temp[,16:18]),rowMeans(temp[,17:18]),rowMeans(temp[,18:19]),rowMeans(temp[,17:22]),rowMeans(temp[,18:20]),rowMeans(temp[,19:21]),rowMeans(temp[,19:20]),rowMeans(temp[,20:21]),rowMeans(temp[,21:22]),rowMeans(temp[,22:23])) # finish to include seasonal means,etc
colnames(seasons) <- c("pYear","pJJA","pAMJJAS", "pD_cA","pSO","DJF", "JFM","FMA","MA","MAM","AM","MJ","AMJJAS","MJJ","JJA","JJ","JA","AS","SO")
seasons <- ts(seasons,start=start(temp)[1])
temp2 <- ts.union(temp,seasons)
colnames(temp2) <- c(colnames(temp),colnames(seasons))
correlations <- cor(temp2[,1],temp2[,2:ncol(temp2)],use="pairwise.complete.obs")
return(correlations)
}


## load detr()
source("detr.R")
source("allfunctionsdavid.v13.R")

### load data

# tree-ring data + detrend

# Kent Holsinger - 2015-10-11
# setwd("C:/Users/mekevans/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/treerings/FalkHolsingerProject/MCN2015Data/MCNIncrements/plots")

setwd("plot-data")

files <- list.files(pattern = "rwl")
data.tr <- vector(mode = "list", length = length(files)); data.tr.det <- data.tr
data.tr.det.mean <- ts(matrix(data = NA, ncol = length(files), nrow = 65), start = 1950, frequency = 1); colnames(data.tr.det.mean) <- files


for(i in 1:length(files)){
	a <- read.rwl(files[i], header = FALSE) # header = T misses data from 1950s for first core
	data.tr[[i]] <- ts(a[,2:ncol(a)], start = as.numeric(rownames(a)[1]), frequency = 1)
	data.tr.det[[i]] <- detr(data.tr[[i]],10)
	data.tr.det.mean[,i] <- apply(data.tr.det[[i]],1,mean,na.rm = T)
	}


# climate data

# Kent Holsinger - 2015-10-11
# setwd("C:/Users/mekevans/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/treerings/FalkHolsingerProject/extractclimate")
setwd("..")

prism <- read.table("PRISM_MCN.txt", header = T)      #starts Jan 1981
vpd <- ts(read.table("vpd_MCC_cru1901_1913.txt", header = T), start = 1901, frequency = 12)
vpd81 <- window(vpd, start = 1981, end = 2014, extend = T)

data.climate <- cbind(prism, c(as.vector(vpd81), rep(NA,11)))
colnames(data.climate) <- c(colnames(prism), "vpd")


months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
prc <- data.climate[,3]; dim(prc) <- c(12,length(prc)/12); prc.ts <- ts(t(prc), start = 1981, frequency = 1); colnames(prc.ts) <- months
Tmean <- data.climate[,4]; dim(Tmean) <- c(12,length(Tmean)/12); Tmean.ts <- ts(t(Tmean), start = 1981, frequency = 1); colnames(Tmean.ts) <- months
Tmax <- data.climate[,5]; dim(Tmax) <- c(12,length(Tmax)/12); Tmax.ts <- ts(t(Tmax), start = 1981, frequency = 1); colnames(Tmax.ts) <- months
vpdM <- data.climate[,6]; dim(vpdM) <- c(12,length(vpdM)/12); vpdM.ts <- ts(t(vpdM), start = 1981, frequency = 1); colnames(vpdM.ts) <- months



### climate correlations

correls.prc <- vector(mode = "list", length = length(files))
correls.Tmean <- correls.prc; correls.Tmax <- correls.prc; correls.vpdM <- correls.prc

for(i in 1:length(correls.prc)){
	correls.prc[[i]] <- monthly.correlations(data.tr.det.mean[,i], prc.ts, 1981, 2014)
	correls.Tmean[[i]] <- monthly.correlations(data.tr.det.mean[,i], Tmean.ts, 1981, 2014)
	correls.Tmax[[i]] <- monthly.correlations(data.tr.det.mean[,i], Tmax.ts, 1981, 2014)
	correls.vpdM[[i]] <- monthly.correlations(data.tr.det.mean[,i], vpdM.ts, 1981, 2014)
}

# Kent Holsinger - 2015-10-11
# all.correls.prc <- rbind(correls.prc[[1]],correls.prc[[2]],correls.prc[[3]],correls.prc[[4]],correls.prc[[5]], correls.prc[[6]]); rownames(all.correls.prc) <- files
# all.correls.Tmean <- rbind(correls.Tmean[[1]],correls.Tmean[[2]],correls.Tmean[[3]],correls.Tmean[[4]],correls.Tmean[[5]], correls.Tmean[[6]]); rownames(all.correls.Tmean) <- files
# all.correls.Tmax <- rbind(correls.Tmax[[1]],correls.Tmax[[2]],correls.Tmax[[3]],correls.Tmax[[4]],correls.Tmax[[5]], correls.Tmax[[6]]); rownames(all.correls.Tmax) <- files
# all.correls.vpdM <- rbind(correls.vpdM[[1]],correls.vpdM[[2]],correls.vpdM[[3]],correls.vpdM[[4]],correls.vpdM[[5]], correls.vpdM[[6]]); rownames(all.correls.vpdM) <- files
all.correls.prc <- matrix(unlist(correls.prc), nrow=length(correls.prc), byrow=TRUE); colnames(all.correls.prc[[1]]); rownames(all.correls.prc) <- files
all.correls.Tmean <- matrix(unlist(correls.Tmean), nrow=length(correls.prc), byrow=TRUE); colnames(all.correls.prc[[1]]); rownames(all.correls.prc) <- files
all.correls.Tmax <- matrix(unlist(correls.Tmax), nrow=length(correls.prc), byrow=TRUE); colnames(all.correls.prc[[1]]); rownames(all.correls.prc) <- files
all.correls.vpdM <- matrix(unlist(correls.vpdM), nrow=length(correls.prc), byrow=TRUE); colnames(all.correls.prc[[1]]); rownames(all.correls.prc) <- files


### climate response functions

crf <- vector(mode = "list", length = length(files))
crf.tp <- crf

for(i in 1:length(crf)){
	forcrf.tr <- data.frame(window(data.tr.det.mean, start = 1981, end = 2013)[,i]); rownames(forcrf.tr) <- 1981:2013; colnames(forcrf.tr) <- "TRW"
	forcrf.climate <- data.climate[1:396,-5]
	forcrf.tp.climate <- data.climate[1:396,-c(5,6)]
	crf[[i]] <- dcc(forcrf.tr, forcrf.climate, start = 1, end = 9)
	crf.tp[[i]] <- dcc(forcrf.tr, forcrf.tp.climate, start = -4, end = 9)
}


### write data

# Kent Holsinger - 2015-10-11
# setwd("/Users/Flurin/Documents/post_doc_Tucson/Jemez_Margaret/climate_correlations")
setwd("climate-correlations")

write.table(all.correls.prc, file = "monthly_seasonal_correlations_TRW_precip_1981_2014.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(all.correls.Tmean, file = "monthly_seasonal_correlations_TRW_Tmean_1981_2014.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(all.correls.Tmax, file = "monthly_seasonal_correlations_TRW_Tmax_1981_2014.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(all.correls.vpdM, file = "monthly_seasonal_correlations_TRW_VPD_1981_2014.txt", sep = "\t", col.names = T, row.names = T, quote = F)

# Kent Holsinger - 2015-10-11
# setwd("/Users/Flurin/Documents/post_doc_Tucson/Jemez_Margaret/climate_response_functions")
setwd("../climate-response-functions")

for(i in 1:length(crf)){
	write.table(crf[[i]], file = paste("climate_response_function_1981_2014_", files[i], sep = ""), sep = "\t", col.names = T, row.names = T, quote = F)
	write.table(crf.tp[[i]], file = paste("climate_response_function_TPlag_1981_2014_", files[i], sep = ""), sep = "\t", col.names = T, row.names = T, quote = F)
}



### plots

# Kent Holsinger - 2015-10-11
# setwd("/Users/Flurin/Documents/post_doc_Tucson/Jemez_Margaret/plots")
setwd("../plots")

nameofpicturefile = "climate_correlations.pdf"
pdf(file = nameofpicturefile, width = 15, height = 6)

barplot(all.correls.prc[1,], las = 3, ylim = c(-0.6,0.6), ylab = "correlation", main = "precipitation")
barplot(all.correls.Tmean[1,], las = 3, ylim = c(-0.6,0.6), ylab = "correlation", main = "mean T")
barplot(all.correls.Tmax[1,], las = 3, ylim = c(-0.6,0.6), ylab = "correlation", main = "max T")
barplot(all.correls.vpdM[1,], las = 3, ylim = c(-0.6,0.6), ylab = "correlation", main = "VPD")

dev.off()


nameofpicturefile = "climate_response_functions.pdf"
pdf(file = nameofpicturefile, width = 10, height = 6)

for(i in 1:length(crf)){dcplot(crf[[i]])}

dev.off()

nameofpicturefile = "climate_response_functions_TPlag.pdf"
pdf(file = nameofpicturefile, width = 10, height = 6)

for(i in 1:length(crf.tp)){dcplot(crf.tp[[i]])}

dev.off()















