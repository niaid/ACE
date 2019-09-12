# see Ca/Co counts in start data
data <- read.table("plink_start/TB_GWAS.fam")
names(data) <- c("fam", "ind", "father", "mother", "sex", "affected")

#plot assoc1

data <- read.table("assoc1/TB_GWAS.assoc.fisher", head = T)
names(data)
summary(data$P)
min(data$P)
#if this returns 0 can do. looks like it was 1.047e-92
min(data$P[data$P != 0])
require(qqman)
manhattan(data, col = c("orange", "blue"), ylim = c(0, 95))

#make qq plot also
qq(data$P)

#plot pca1

data <- read.table("pca1/TB_GWAS.eigenvec", head = T)
##plot PC1 vs PC2
plot(data$PC1, data$PC2)

#it would be useful to color the points by Ca/Co status

pheno <- read.table("qc1/TB_GWAS.fam")
head(pheno)
names(pheno) <- c("fam", "ind", "father", "mother", "sex", "affected")
head(pheno)

case <- pheno$affected == 2
control <- pheno$affected == 1
plot(data$PC1, data$PC2, type = "n", xlab = "PC1", ylab = "PC2", main = "TB pca")
points(data$PC1[control], data$PC2[control], col = rgb(0,0,1, 0.3), pch = 20)
points(data$PC1[case], data$PC2[case], col = rgb(1,0,0, 0.3), pch = 20)
abline(v = 0.008)

#subset to just PC1 > 0.008

keep_data <- subset(data, PC1 > 0.008, select=c(FID, IID))
write.table(keep_data, "ancestry_filt1/samples_to_keep.txt", quote = F, row.names = F, col.names = F)

# go back to plink script to subset


## plot pca2
data <- read.table("pca2/TB_GWAS.eigenvec", head = T)
pheno <- read.table("ancestry_filt1/TB_GWAS.fam")
names(pheno) <- c("fam", "ind", "father", "mother", "sex", "affected")
case <- pheno$affected == 2
control <- pheno$affected == 1
plot(data$PC1, data$PC2, type = "n", xlab = "PC1", ylab = "PC2", main = "TB pca after ancestry filtering")
points(data$PC1[control], data$PC2[control], col = rgb(0,0,1, 0.3), pch = 20)
points(data$PC1[case], data$PC2[case], col = rgb(1,0,0, 0.3), pch = 20)

## still substructure in pca plot. Now try 1:1 match using optmatch

# optmatch requires 0 for control and 1 for cases
data$CaCo <- pheno$affected - 1
table(data$CaCo)

require(optmatch)
pm <- pairmatch(CaCo ~ PC1 + PC2 + PC3, controls = 1, data = data, remove.unmatchables = T)
data2 <- data[which(!is.na(pm)),]
table(data2$CaCo)

#plot after matching
case <- data2$CaCo == 1
control <- data2$CaCo == 0
plot(data2$PC1, data2$PC2, type = "n", xlab = "PC1", ylab = "PC2", main = "TB pca after matching")
points(data2$PC1[control], data2$PC2[control], col = rgb(0,0,1, 0.3), pch = 20)
points(data2$PC1[case], data2$PC2[case], col = rgb(1,0,0, 0.3), pch = 20)

#output samples to subset
head(data2)
head(data2[,1:2])

write.table(data2[,1:2], "ancestry_filt2/samples_to_keep.txt", quote = F, row.names = F, col.names = F)

# go back to plink to subset

# plot assoc3

data <- read.table("assoc3/TB_GWAS.assoc.fisher", head = T)
names(data)
summary(data$P)
min(data$P)
#if this returns 0 can do. looks like it was 2.787e-42
min(data$P[data$P != 0])

require(qqman)
manhattan(data, col = c("green3", "gold", "red2"), ylim = c(0, 45))

#make qq plot also
qq(data$P)
