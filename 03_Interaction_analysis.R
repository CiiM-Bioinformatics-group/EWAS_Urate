# Interaction analysis
# Input: DNA methylation data, Urate value at each time point
# Output: Summary statistic of the interaction analysis
library(data.table)
library(MASS) 
library(sandwich) 
library(lmtest) 
library(parallel) 
library(R.utils)
library(openxlsx)
library(qqman)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


load("input/m1.trim.rdata")
load("input/uv1.rdata")
load("input/pdv1.rdata")

mv1.t <- t(m1.trim)

sum(pdv1$id == uv1$id)
sum(pdv1$id == substr(rownames(mv1.t),1,9))


RLMtest = function(meth_matrix, methcol, urate, age, gender, plate, CD8, CD4, NK, B, Mono, Neu) {mod = try(rlm(meth_matrix[, methcol] ~ urate+age+gender+plate+CD8+CD4+NK+B+Mono+Neu+urate*gender, maxit=200))
	cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))

if (class(cf)=="try-error") {
  bad <- as.numeric(rep(NA, 3))
  names(bad)<- c("Estimate", "Std. Error", "Pr(>|z|)")
  bad
}
else{
  cf[20, c("Estimate", "Std. Error", "Pr(>|z|)")]
}
}

res <- lapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest, meth_matrix=mv1.t, urate=uv1$con, age=pdv1$age, gender=pdv1$gender, plate=pdv1$Sample_Plate, CD8=pdv1$CD8T, CD4=pdv1$CD4T, NK=pdv1$NK, B=pdv1$Bcell, Mono=pdv1$Mono, Neu=pdv1$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_,4))
setattr(res, "names", make.names(names(res), unique=TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result<-data.table(result)
result[, probeID := probelistnamesB]
setnames(result, c("BETA","SE", "P_VAL", "probeID")) # rename columns
setcolorder(result, c("probeID","BETA","SE", "P_VAL"))
result$padj <- p.adjust(result$P_VAL, method ="BH")

lambda = median(qchisq(as.numeric(as.character(result$P_VAL)),df=1,lower.tail = F),na.rm=T)/qchisq(0.5,1)
lambda 

write.xlsx(result, file = "output/v1_mtrimmed_interaction.rlm.xlsx")

