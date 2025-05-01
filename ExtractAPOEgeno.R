apoe_w5_geno <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Covar_Pheno/APOE/W5_APOE_genotypes.GT.FORMAT",header=F,sep="\t")
apoe_w8_geno <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Covar_Pheno/APOE/W8_APOE_genotypes.GT.FORMAT",header=F,sep="\t")
apoe_w5_geno_t <- as.data.frame(t(apoe_w5_geno))
apoe_w8_geno_t <- as.data.frame(t(apoe_w8_geno))

# Remove the first row
apoe_w5_geno_t <- apoe_w5_geno_t[-1, ]
apoe_w8_geno_t <- apoe_w8_geno_t[-1, ]
# Set the second row as column names
colnames(apoe_w5_geno_t) <- apoe_w5_geno_t[1, ]
colnames(apoe_w8_geno_t) <- apoe_w8_geno_t[1, ]
# Remove the now redundant second row
apoe_w5_geno_t <- apoe_w5_geno_t[-1, ]
apoe_w8_geno_t <- apoe_w8_geno_t[-1, ]

colnames(apoe_w5_geno_t) <- c("IID","rs429358","rs7412")
colnames(apoe_w8_geno_t) <- c("IID","rs429358","rs7412")

library(tidyr)
apoe_w5_geno_t <- separate(apoe_w5_geno_t, rs429358, into = c("rs429358_1", "rs429358_2"), sep = "\\|")
apoe_w8_geno_t <- separate(apoe_w8_geno_t, rs429358, into = c("rs429358_1", "rs429358_2"), sep = "\\|")
apoe_w5_geno_t <- separate(apoe_w5_geno_t, rs7412, into = c("rs7412_1", "rs7412_2"), sep = "\\|")
apoe_w8_geno_t <- separate(apoe_w8_geno_t, rs7412, into = c("rs7412_1", "rs7412_2"), sep = "\\|")

#rs429... ref=0=T
#rs7412 ref=0=C
apoe_geno <- rbind(apoe_w5_geno_t, apoe_w8_geno_t)
apoe_geno$hap1 <- ifelse(apoe_geno$rs429358_1 == 1 & apoe_geno$rs7412_1 == 0, 1, 0)
apoe_geno$hap2 <- ifelse(apoe_geno$rs429358_2 == 1 & apoe_geno$rs7412_2== 0, 1, 0)
apoe_geno$e4 <- apoe_geno$hap1+apoe_geno$hap2

apoe_for_plink <- data.frame(0,apoe_geno$IID,apoe_geno$e4)
colnames(apoe_for_plink) <- c("FID","IID","APOE4")
write.table(apoe_for_plink,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_status_pheno_data.txt",col.names = TRUE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')


