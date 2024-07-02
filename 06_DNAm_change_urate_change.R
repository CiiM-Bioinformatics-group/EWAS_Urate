# Epigenome-wide association between DNA methylation change and Urate change
# We need to input DNA methylation change and urate change
# We will get the summary statistics of the association analysis

library(readxl)
library(stringr)
library(data.table)
library(MASS) 
library(sandwich) 
library(lmtest) 
library(parallel) 
library(R.utils)
library(openxlsx)
library(dplyr)

load("/input/mv2_v1_regress.rdata")
mv1.t <- t(mv2_v1_regress)

load("input/uv1.rdata")
load("input/uv2.rdata")
load("input/pdv1.rdata")

uv2_minus_uv1 <- data.frame(uv2 = uv2$con, uv1=uv1$con)
uv2_minus_uv1$minus <- uv2_minus_uv1$uv2-uv2_minus_uv1$uv1
rownames(uv2_minus_uv1) <- uv1$id
pdv1$Sample_Plate <- droplevels(pdv1$Sample_Plate)
str(pdv1)

sum(pdv1$id == rownames(uv2_minus_uv1))
sum(pdv1$id == rownames(mv1.t))
pdv1$age <- as.numeric(pdv1$age)

mv1.t <- as.data.frame(mv1.t)

RLMtest = function(meth_matrix, methcol, urate, age, sex) {mod = try(rlm(meth_matrix[, methcol] ~ urate+age+sex, maxit=200))
	cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))

if (class(cf)=="try-error") {
  bad <- as.numeric(rep(NA, 3))
  names(bad)<- c("Estimate", "Std. Error", "Pr(>|z|)")
  bad
}
else{
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
}
}

res <- mclapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest, meth_matrix=mv1.t, urate=uv2_minus_uv1$minus, age=pdv1$age, sex = pdv1$gender)

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
result$lambda <- lambda 
save(result, file = "output/methy_regress_urate_change/all_v2_v1_methy_regress.rdata")










