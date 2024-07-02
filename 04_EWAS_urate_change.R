# Urate change EWAS
# We will start with the input of urate level from two time points, and baseline DNA methylation M value
# We will get the summary statistics of the EWAS of urate change. 
library(data.table)
library(MASS) 
library(sandwich) 
library(lmtest) 
library(parallel) 
library(R.utils)
library(openxlsx)
library(dplyr)
library(qqman)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

load("input/m1.trim.rdata")
load("input/uv1.rdata")
load("input/uv2.rdata")
load("input/pdv1.rdata")

mv1.t <- t(m1.trim)

uv2_2 <- uv2[which(uv2$name %in% uv1$name),]
sum(uv2_2$name == uv1$name)
uv2_minus_uv1 <- data.frame(uv2_2 = uv2_2$con, uv1=uv1$con)
uv2_minus_uv1$minus <- uv2_minus_uv1$uv2_2-uv2_minus_uv1$uv1

uv2_minus_uv1$minus <- qnorm((rank(uv2_minus_uv1$minus, na.last="keep") - 0.5) / sum(!is.na(uv2_minus_uv1$minus)))

pdv1$Sample_Plate <- droplevels(pdv1$Sample_Plate)
str(pdv1)

RLMtest = function(meth_matrix, methcol, urate, age, gender, plate, CD8, CD4, NK, B, Mono, Neu) {mod = try(rlm(urate ~ meth_matrix[, methcol]+age+gender+plate+CD8+CD4+NK+B+Mono+Neu, maxit=200))
if(class(mod) == "try-error"){
  print(paste("error thrown by column", methcol))
  invisible(rep(NA, 3))
}else cf = coeftest(mod, vcov=vcovHC(mod, type="HC0"))
cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
}

res <- mclapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest, meth_matrix=mv1.t, urate=uv2_minus_uv1$minus, age=pdv1$age, gender=pdv1$gender, plate=pdv1$Sample_Plate, CD8=pdv1$CD8T, CD4=pdv1$CD4T, NK=pdv1$NK, B=pdv1$Bcell, Mono=pdv1$Mono, Neu=pdv1$Neu)

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

write.xlsx(result, file = "output/all_v2_minus_v1_mval_trimmed_urate_rank.xlsx")

ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annotdata<-ann850k[c( "chr","pos","UCSC_RefGene_Name","Relation_to_Island")]
result <- data.frame(result)
rownames(result) <- result$probeID
Ml<-data.frame(merge(result,annotdata,by="row.names"))
Ml$BP<-as.integer(Ml$pos)
Ml$SNP<-as.character(Ml$Row.names)
Ml$CHR<-as.numeric(gsub("[^0-9]", "", (Ml$chr)))
Ml$P<-Ml$P_VAL
Ml1<-Ml[,c("CHR","BP","P","SNP")]

sum(is.na(Ml1))/4
Ml1 <- na.omit(Ml1)

pdf("output/qqplot_all_v2_minus_v1_mval_trimmed_urate_rank.pdf")
qq(Ml1$P, main = "Q-Q plot of v2 minus v1(mval_trimmed)")
dev.off()

pdf("output/manhan_all_v2_minus_v1_mval_trimmed_urate_rank.pdf")
manhattan(Ml1, main = "Manhattan Plot", cex = 0.6, cex.axis = 0.9, col = c("lightblue","royalblue"), suggestiveline = F)
dev.off()

lambda = median(qchisq(as.numeric(as.character(Ml1$P)),df=1,lower.tail = F),na.rm=T)/qchisq(0.5,1)
lambda 




