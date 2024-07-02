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

load("/input/mv1.t.rata")
load("/input/pdv1.rdata")

phe.use <- pdv1 %>% select(plate, CD8T, CD4T, NK, Bcell, Mono, Neu)
sum(substr(rownames(phe.use),1,9) == substr(rownames(mv1.t), 2, 10))

LMtest = function(meth_matrix, methcol, plate, CD8, CD4, NK, B, Mono, Neu) {
	mod = lm(meth_matrix[, methcol] ~ plate+CD8+CD4+NK+B+Mono+Neu)
	mod$residuals	
	}

res <- lapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), LMtest, meth_matrix=mv1.t, plate=pdv1$Sample_Plate, CD8=pdv1$CD8T, CD4=pdv1$CD4T, NK=pdv1$NK, B=pdv1$Bcell, Mono=pdv1$Mono, Neu=pdv1$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_,283))
setattr(res, "names", make.names(names(res), unique=TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result<-data.table(result)
colnames(result) <- rownames(mv1.t)
result[, probeID := probelistnamesB]

t1_raw_regress <- result
save(t1_raw_regress, file = "/input/t1_raw_regress.rdata")

#### Time 2 
load("/input/mv2.rdata")
load("/input/pdv2.rdata")

phe.use <- pdv2 %>% select(plate, CD8T, CD4T, NK, Bcell, Mono, Neu)
sum(substr(rownames(phe.use),1,9) == substr(colnames(mv2), 2, 10))
mv2.t <- t(mv2)
sum(substr(rownames(phe.use),1,9) == substr(rownames(mv2.t), 2, 10))

res <- lapply(setNames(seq_len(ncol(mv2.t)), dimnames(mv2.t)[[2]]), LMtest, meth_matrix=mv2.t, plate=pdv2$Sample_Plate, CD8=pdv2$CD8T, CD4=pdv2$CD4T, NK=pdv2$NK, B=pdv2$Bcell, Mono=pdv2$Mono, Neu=pdv2$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_,283))
setattr(res, "names", make.names(names(res), unique=TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result<-data.table(result)
colnames(result) <- rownames(mv1.t)
result[, probeID := probelistnamesB]

t2_raw_regress <- result
save(t2_raw_regress, file = "/input/t2_raw_regress.rdata")

#### Time 3 
load("/input/mv3.rdata")
load("/input/pdv3.rdata")

phe.use <- pdv3 %>% select(plate, CD8T, CD4T, NK, Bcell, Mono, Neu)
sum(substr(rownames(phe.use),1,9) == substr(colnames(mv3), 2, 10))
mv3.t <- t(mv3)
sum(substr(rownames(phe.use),1,9) == substr(rownames(mv1.t), 2, 10))

LMtest = function(meth_matrix, methcol, plate, CD8, CD4, NK, B, Mono, Neu) {
	mod = lm(meth_matrix[, methcol] ~ plate+CD8+CD4+NK+B+Mono+Neu)
	mod$residuals	
	}

res <- lapply(setNames(seq_len(ncol(mv3.t)), dimnames(mv3.t)[[2]]), LMtest, meth_matrix=mv3.t, plate=pdv3$Sample_Plate, CD8=pdv3$CD8T, CD4=pdv3$CD4T, NK=pdv3$NK, B=pdv3$Bcell, Mono=pdv3$Mono, Neu=pdv3$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_,283))
setattr(res, "names", make.names(names(res), unique=TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result<-data.table(result)
colnames(result) <- rownames(mv1.t)
result[, probeID := probelistnamesB]

t3_raw_regress <- result
save(t3_raw_regress, file = "/input/t3_raw_regress.rdata")


sum(substr(colnames(t1_raw_regress),2,10)== substr(colnames(t2_raw_regress),2,10))
sum(t1_raw_regress$probeID == t2_raw_regress$probeID)
t1_raw_regress <- data.frame(t1_raw_regress)
t2_raw_regress <- data.frame(t2_raw_regress)
rownames(t1_raw_regress) <- t1_raw_regress$probeID
rownames(t2_raw_regress) <- t2_raw_regress$probeID
t1_raw_regress <- t1_raw_regress[,-282]
t2_raw_regress <- t2_raw_regress[,-282]
mv2_v1_regress <- t2_raw_regress - t1_raw_regress
save(mv2_v1_regress, file = "/input/mv2_v1_regress.rdata")

colnames(t1_raw_regress) <- substr(colnames(t1_raw_regress),2,10)
colnames(t2_raw_regress) <- substr(colnames(t2_raw_regress),2,10)
colnames(t3_raw_regress) <- substr(colnames(t3_raw_regress),2,10)

intersect(colnames(t1_raw_regress), colnames(t3_raw_regress))
intersect(colnames(t2_raw_regress), colnames(t3_raw_regress))

t1_raw_regress_2 <- t1_raw_regress[,which(colnames(t1_raw_regress) %in% colnames(t3_raw_regress))]
t2_raw_regress_2 <- t2_raw_regress[,which(colnames(t2_raw_regress) %in% colnames(t3_raw_regress))]
rownames(t3_raw_regress) <- t3_raw_regress$robeID

t3_raw_regress <- t3_raw_regress[,-281]

sum(colnames(t1_raw_regress_2) == colnames(t3_raw_regress))
sum(rownames(t1_raw_regress_2) == rownames(t3_raw_regress))
sum(colnames(t2_raw_regress_2) == colnames(t3_raw_regress))
sum(rownames(t2_raw_regress_2) == rownames(t3_raw_regress))

mv3_v1_regress <- t3_raw_regress - t1_raw_regress_2
save(mv3_v1_regress, file = "/input/mv3_v1_regress.rdata")

mv3_v2_regress <- t3_raw_regress - t2_raw_regress_2
save(mv3_v2_regress, file = "/input/mv3_v2_regress.rdata")
















