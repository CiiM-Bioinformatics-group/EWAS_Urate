# We started with the 300BCG DNA methylation data after QC 
# We will get the DNA methylation data after trimming at each time point (day0, day14, day90)
library(minfi)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)

M.val <- readRDS("input/mVals_filtered.rds")

removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}


system.time(OutlierResults<-removeOutliers(M.val))
M.val<-OutlierResults[[1]]

save(M.val, file="input/mVals_filtered_trimmed.Rdata")
save(Log, file="input/mVals_filtered_trimmed_log.Rdata")

load("input/pdv1.rdata")
load("input/pdv2.rdata")
load("input/pdv3.rdata")

m1.trim <- M.val[,which(colnames(M.val) %in% rownames(pdv1))] # 283 samples
m2.trim <- M.val[,which(colnames(M.val) %in% rownames(pdv2))] # 283 samples
m3.trim <- M.val[,which(colnames(M.val) %in% rownames(pdv3))] # 282 samples

save(m1.trim, file = "input/m1.trim.rdata")
save(m2.trim, file = "input/m2.trim.rdata")
save(m3.trim, file = "input/m3.trim.rdata")


