library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library(data.table)
library(rsq)

#Read in biomarker phenotype
library(readxl)
biomarker <- read_excel("C:/Users/luwan/Desktop/UGA/KY_lab/Regression/Dementia_W8_n436.xlsx", na = "#NULL!")
#biomarker_clean <- biomarker[complete.cases(biomarker[c("AgePCSCw8", "pTau181_W8", "Abeta40W8","Abeta42W8", "GFAPW8", "NFlightW8")]), ]
biomarker_clean <- biomarker[complete.cases(biomarker["AgePCSCw8"]), ]

#ancestry proportion
admix4results <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Admixture/Gibbons_flipped_admix_allchr.4.Q",header=F,sep=" ")
admix4IID <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Admixture/Gibbons_flipped_admix_allchr.fam",header=F, col.names=c('V0','Individual', 'V2','V3','V4','V5'),sep="\t")
admix4proportion <- cbind(admix4results, admix4IID$Individual)
colnames(admix4proportion) <- c("pop1", "pop2", "pop3", "pop4", "IID")

#sex
W5sex <- read_excel("C:/Users/luwan/Desktop/UGA/KY_lab/Expectation_check_b_g/w5_expectation.xlsx")
W8sex <- read_excel("C:/Users/luwan/Desktop/UGA/KY_lab/Expectation_check_b_g/w8_expectation.xlsx")
colnames(W5sex)[colnames(W5sex) == "Relationship_W5"] <- "Relationship"
colnames(W8sex)[colnames(W8sex) == "Relationship_W8"] <- "Relationship"
sexinfo <- rbind(W5sex, W8sex)

#combine data
admix_sex <- merge(admix4proportion, sexinfo[c("PID", "IID", "gsex")], by.x = "IID", by.y = "IID")
admix_sex_age_biomarker <- merge(admix_sex, biomarker_clean[c("PID", "AgePCSCw8", "pTau181_W8", "Abeta40W8","Abeta42W8", "GFAPW8", "NFlightW8")], by.x = "PID", by.y = "PID")
df_combo <- admix_sex_age_biomarker
df_combo$gsex <- factor(df_combo$gsex)

#W5 phenotype data
biomarker_w5 <- read_excel("C:/Users/luwan/Desktop/UGA/KY_lab/Regression/Dementia_W5_n559.xls", na = "#NULL!")
W5_age <- read_excel("C:/Users/luwan/Desktop/UGA/KY_lab/Regression/PCSC-age5.xlsx", na = "#NULL!")
biomarker_w5_age <- merge(biomarker_w5, W5_age[c("PID", "AgePCSCw5.5")], by.x = "PID", by.y = "PID")
colnames(biomarker_w5_age)[colnames(biomarker_w5_age) == "AgePCSCw5.5"] <- "Age"
admix_sex_age_biomarker_W5 <- merge(admix_sex, biomarker_w5_age[c("PID", "Age", "pTau181_W5", "Abeta40W5","Abeta42W5", "GFAPW5", "NFlightW5")], by.x = "PID", by.y = "PID")







#Remove NAs
admix_sex_age_biomarker_W5_clean <- admix_sex_age_biomarker_W5[complete.cases(admix_sex_age_biomarker_W5[, c("GFAPW5")]), ]
sum(is.na(admix_sex_age_biomarker_W5$GFAPW5))
sum(is.na(admix_sex_age_biomarker_W5_clean$NFlightW5))
admix_sex_age_biomarker_W5_clean <- merge(admix_sex_age_biomarker_W5_clean, sexinfo[c("PID", "Race")], by.x = "PID", by.y = "PID")
admix_sex_age_biomarker_W5_clean <- admix_sex_age_biomarker_W5_clean[!duplicated(admix_sex_age_biomarker_W5_clean$PID), ]

#admix_sex_age_biomarker_clean <- admix_sex_age_biomarker[complete.cases(admix_sex_age_biomarker[, c("GFAPW8")]), ]
admix_sex_age_biomarker_clean <- merge(admix_sex_age_biomarker, sexinfo[c("PID", "Race")], by.x = "PID", by.y = "PID")
#admix_sex_age_biomarker_clean <- admix_sex_age_biomarker_clean[!duplicated(admix_sex_age_biomarker_clean$PID), ]







#change name
mergeW8 <- df_combo
mergeW8 <- setNames(mergeW8, c("PID", "IID", "pop1","AFR","pop3","pop4", "gsex","Age","ptau181","Abeta40","Abeta42","GFAP","NFlight"))
mergeW5 <- admix_sex_age_biomarker_W5
mergeW5 <- setNames(mergeW5, c("PID", "IID", "pop1","AFR","pop3","pop4", "gsex","Age","ptau181","Abeta40","Abeta42","GFAP","NFlight"))


# Filter df2 to only include rows where column2 is not in df1's column2
filtered_mergeW5 <- anti_join(mergeW5, mergeW8, by = "IID")

# Merge df1 with the filtered df2
merged_W5W8_dup <- rbind(mergeW8, filtered_mergeW5)
merged_W5W8_rel <- merged_W5W8_dup[!duplicated(merged_W5W8_dup$PID), ]

#write.table(merged_W5W8_rel$IID, file = "C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/combined_indiv_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


##remove related individuals
rm_list <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/W5W8_min_related_list_to_remove.txt",header=F,sep="\t")
merged_W5W8_allAFR <- merged_W5W8_rel[!merged_W5W8_rel$IID %in% rm_list$V2, ]
#Remove missing biomarker
merged_W5W8_allAFR_nomissing <- merged_W5W8_allAFR[complete.cases(merged_W5W8_allAFR[, c("GFAP")]), ]
sum(is.na(merged_W5W8_allAFR_nomissing$ptau181))



#Count Race
merged_W5W8_allAFR_nomissing_race <- merge(merged_W5W8_allAFR_nomissing, sexinfo[c("PID", "Race")], by.x = "PID", by.y = "PID")
merged_W5W8_allAFR_nomissing_race <- merged_W5W8_allAFR_nomissing_race[!duplicated(merged_W5W8_allAFR_nomissing_race$PID), ]
table(merged_W5W8_race$Race)
mean(merged_W5W8_allAFR_nomissing_race$ptau181,na.rm =TRUE)


# #Filter out AFR<0.1
# merged_W5W8 <- merged_W5W8_allAFR[merged_W5W8_allAFR$AFR >= 0.1, ]
# table(merged_W5W8_race$Race)
##add race and count individuals
# merged_W5W8_race <- merge(merged_W5W8, sexinfo[c("PID", "Race")], by.x = "PID", by.y = "PID")
# merged_W5W8_race <- merged_W5W8_race[!duplicated(merged_W5W8_race$PID), ]
# table(merged_W5W8_race$Race)

# #Save total list for PCA
# all_keep_indiv_unrel_0.9 <- data.frame(0, merged_W5W8$IID)
# write.table(all_keep_indiv_unrel_0.9, file = "C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/all_keep_indiv_unrel_0.9.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
all_keep_indiv_unrel <- data.frame(0, merged_W5W8_allAFR_nomissing$IID)
write.table(all_keep_indiv_unrel, file = "C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/all_keep_indiv_unrel.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

merged_W5W8 <- merged_W5W8_allAFR_nomissing_race

#Log transformation
merged_W5W8$log_pTau181 <- log(merged_W5W8$ptau181)
merged_W5W8$log_GFAP <- log(merged_W5W8$GFAP)
merged_W5W8$log_NFlight <- log(merged_W5W8$NFlight)

#Categorize Abeta40W8_age and Abeta42W8_age
merged_W5W8$Abeta40_cat <- ifelse(merged_W5W8$Abeta40 < 3.59, "0", "1")
merged_W5W8$Abeta42_cat <- ifelse(merged_W5W8$Abeta42 < 1.38, "0", "1")
merged_W5W8$Abeta40_cat <- factor(merged_W5W8$Abeta40_cat,ordered = F,levels=c("0","1"))
merged_W5W8$Abeta42_cat <- factor(merged_W5W8$Abeta42_cat,ordered = F,levels=c("0","1"))


#Check phenotype
library(ggplot2)
pTau181_hist <- ggplot(merged_W5W8, aes(x = ptau-181)) + 
  geom_histogram(binwidth = 1, fill = "black", color = "black") +
  labs(x = "pTau181", y = "# of Individuals")
Abeta40_hist <- ggplot(merged_W5W8, aes(x = Abeta40)) + 
  geom_histogram(binwidth = 1, fill = "black", color = "black") +
  labs(x = "Abeta40", y = "# of Individuals")
Abeta42_hist <- ggplot(merged_W5W8, aes(x = Abeta42)) + 
  geom_histogram(fill = "black", color = "black") +
  labs(x = "Abeta42", y = "# of Individuals")
GFAP_hist <- ggplot(merged_W5W8, aes(x = GFAP)) + 
  geom_histogram(binwidth = 1, fill = "black", color = "black") +
  labs(x = "GFAP", y = "# of Individuals")
NFlight_hist <- ggplot(merged_W5W8, aes(x = NFlight)) + 
  geom_histogram(binwidth = 1, fill = "black", color = "black") +
  labs(x = "NFL", y = "# of Individuals")
# hist_plot <- pTau181_hist + Abeta40_hist + Abeta42_hist + GFAP_hist + NFlight_hist +
hist_plot <- pTau181_hist + GFAP_hist + NFlight_hist +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "Biomarker Levels at Sample Collection", tag_levels = 'A')
print(hist_plot)

cat_Abeta40 <- ggplot(merged_W5W8, aes(x = Abeta40_cat)) + 
  geom_bar() +
  geom_text(stat = 'count', aes(label = ..count.., y = ..count..), vjust = -0.5)+
  labs(x = "Abeta40", y = "# of Individuals")
cat_Abeta42 <- ggplot(merged_W5W8, aes(x = Abeta42_cat)) + 
  geom_bar() +
  geom_text(stat = 'count', aes(label = ..count.., y = ..count..), vjust = -0.5)+
  labs(x = "Abeta42", y = "# of Individuals")
cat_biomarker <- cat_Abeta40 + cat_Abeta42 +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "Categorized Biomarker", tag_levels = 'A')
print(cat_biomarker)

#Check age
check_Age <- function(df,pheno,lab_pheno="pheno_name"){
  age_vs_bio <-
    ggplot(df, aes(x = Age, y = pheno)) + 
      geom_point()+
      labs(x = "Age", y = lab_pheno)
  return(age_vs_bio)}
ptau_age <- check_Age(merged_W5W8,merged_W5W8$ptau181,lab_pheno="pTau181")
Abeta40_age <- check_Age(merged_W5W8,merged_W5W8$Abeta40,lab_pheno="Abeta40")
Abeta42_age <- check_Age(merged_W5W8,merged_W5W8$Abeta42,lab_pheno="Abeta42")
GFAP_age <- check_Age(merged_W5W8,merged_W5W8$GFAP,lab_pheno="GFAP")
NFlight_age <- check_Age(merged_W5W8,merged_W5W8$NFlight,lab_pheno="NFlight")
# age_vs_biomarker <- ptau_age + Abeta40_age + Abeta42_age + GFAP_age + NFlight_age +
age_vs_biomarker <- ptau_age + GFAP_age + NFlight_age +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "Age vs Biomarker", tag_levels = 'A')
print(age_vs_biomarker)

#Check biomarker distribution
pTau181_AFR <- ggplot(merged_W5W8, aes(x = AFR, y = ptau181)) + 
  geom_point()+
  labs(x = "AFR ancestry proportion", y = "pTau181")
Abeta40_AFR <- ggplot(merged_W5W8, aes(x = AFR, y = Abeta40)) + 
  geom_point()+
  labs(x = "AFR ancestry proportion", y = "Abeta40")
Abeta42_AFR <- ggplot(merged_W5W8, aes(x = AFR, y = Abeta42)) + 
  geom_point()+
  labs(x = "AFR ancestry proportion", y = "Abeta42")
GFAP_AFR <- ggplot(merged_W5W8, aes(x = AFR, y = GFAP)) + 
  geom_point()+
  labs(x = "AFR ancestry proportion", y = "GFAP")
NFlight_AFR <- ggplot(merged_W5W8, aes(x = AFR, y = NFlight)) + 
  geom_point()+
  labs(x = "AFR ancestry proportion", y = "NFL")

# AFR_vs_biomarker <- pTau181_AFR + Abeta40_AFR + Abeta42_AFR + GFAP_AFR + NFlight_AFR +
AFR_vs_biomarker <- pTau181_AFR + GFAP_AFR + NFlight_AFR +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "AFR ancestry proportion vs Biomarker", tag_levels = 'A')
print(AFR_vs_biomarker)

#check log transformation
library(ggplot2)
log_pTau181_hist <- ggplot(merged_W5W8, aes(x = log_pTau181)) + 
  geom_histogram(binwidth = 0.2, fill = "black", color = "black")+
  labs(x = "log(pTau181)", y = "# of Individuals")
log_GFAP_hist <- ggplot(merged_W5W8, aes(x = log_GFAP)) + 
  geom_histogram(binwidth = 0.2, fill = "black", color = "black")+
  labs(x = "log(GFAP)", y = "# of Individuals")
log_NFlight_hist <- ggplot(merged_W5W8, aes(x = log_NFlight)) + 
  geom_histogram(binwidth = 0.2, fill = "black", color = "black")+
  labs(x = "log(NFL)", y = "# of Individuals")

library(patchwork)
log_biomarker <- log_pTau181_hist + log_GFAP_hist + log_NFlight_hist +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "Log Transformed Biomarker", tag_levels = 'A')
print(log_biomarker)

#transformed vs AFR
check_log_AFR <- function(df, pheno, lab_pheno="pheno_name"){
  log_bio_AFR <-
    ggplot(df, aes(x = AFR, y = pheno)) + 
    geom_point() +
    labs(x = "AFR ancestry proportion", y = lab_pheno) +
    theme(text = element_text(size = 16),  # Adjusts overall text size
          axis.title = element_text(size = 14),  # Specific size for axis titles
          axis.text = element_text(size = 12))  # Specific size for axis text
  return(log_bio_AFR)
}
log_ptau_AFR <- check_log_AFR(merged_W5W8,merged_W5W8$log_pTau181,lab_pheno="log(p-Tau181)")
log_GFAP_AFR <- check_log_AFR(merged_W5W8,merged_W5W8$log_GFAP,lab_pheno="log(GFAP)")
log_NFlight_AFR <- check_log_AFR(merged_W5W8,merged_W5W8$log_NFlight,lab_pheno="log(NFL)")
log_vs_biomarker <- log_ptau_AFR + log_GFAP_AFR + log_NFlight_AFR + 
  plot_layout(guides = 'collect') + # Collects all legends into one
  plot_annotation(title = "AFR vs Log-transformed Biomarker", 
                  tag_levels = 'A',
                  theme = theme(plot.title = element_text(size = 20)))  # Increase and bold the title size
print(log_vs_biomarker)





##Try combined model
# model1: biomarker ~ ancestry% + sex + age
header <- c("sampleSize","phenotype","population","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","pop_rsq","sex_rsq","age_rsq")
write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/RegressionResult/combinedW5W8_unrel_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
for (i in 15:17){
  sampleSize <- sum(!is.na(merged_W5W8[,i]))
  fitModel <- glm(unlist(merged_W5W8[,i])~AFR+gsex+Age,data=merged_W5W8)#,family = "binomial")
  ml_summary <- summary(fitModel)
  pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
  sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
  age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
  rsqValue <- rsq.partial(fitModel, adj = TRUE)
  pop_rsq <- rsqValue$partial.rsq[1]
  sex_rsq <- rsqValue$partial.rsq[2]
  age_rsq <- rsqValue$partial.rsq[3]
  result <- t(as.data.frame(c(sampleSize,colnames(merged_W5W8[i]),"AFR",pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,pop_rsq,sex_rsq,age_rsq)))
  # write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/combinedW5W8_unrel_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/RegressionResult/combinedW5W8_unrel_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
}


## Model for categorized
# model1: biomarker ~ ancestry% + sex + age
header <- c("sampleSize","phenotype","population","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue")
write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/regression/Results/combinedW5W8_unrel_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
for (i in 17:18){
  sampleSize <- sum(!is.na(merged_W5W8[,i]))
  fitModel <- glm(unlist(merged_W5W8[,i])~AFR+gsex+Age,data=merged_W5W8,family = "binomial")
  ml_summary <- summary(fitModel)
  pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
  sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
  age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
  result <- t(as.data.frame(c(sampleSize,colnames(merged_W5W8[i]),"AFR",pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue)))
  write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/combinedW5W8_unrel_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
}





#Wilcox test
#Log transformation
merged_W5W8_allAFR$log_pTau181 <- log(merged_W5W8_allAFR$ptau181)
merged_W5W8_allAFR$log_GFAP <- log(merged_W5W8_allAFR$GFAP)
merged_W5W8_allAFR$log_NFlight <- log(merged_W5W8_allAFR$NFlight)

#Categorize Abeta40W8_age and Abeta42W8_age
merged_W5W8_allAFR$Abeta40_cat <- ifelse(merged_W5W8_allAFR$Abeta40 < 3.59, "0", "1")
merged_W5W8_allAFR$Abeta42_cat <- ifelse(merged_W5W8_allAFR$Abeta42 < 1.38, "0", "1")
merged_W5W8_allAFR$Abeta40_cat <- factor(merged_W5W8_allAFR$Abeta40_cat,ordered = F,levels=c("0","1"))
merged_W5W8_allAFR$Abeta42_cat <- factor(merged_W5W8_allAFR$Abeta42_cat,ordered = F,levels=c("0","1"))

sorted_EUR <- merged_W5W8_allAFR[order(-merged_W5W8_allAFR$pop1), ]  # Sorts the DataFrame in descending order based on 'target_column'
sorted_AFR <- merged_W5W8_allAFR[order(-merged_W5W8_allAFR$AFR), ]
top90_EUR <- sorted_EUR[sorted_EUR$pop1 > 0.9, ]
top90_EUR$pop <- "EUR"
top90_AFR <- sorted_AFR[sorted_AFR$AFR > 0.9, ]
top90_AFR$pop <- "AFR"
combined_top90 <- rbind(top90_EUR,top90_AFR)
wilcox.test(top90_EUR$ptau181, top90_AFR$ptau181, paired = FALSE)
wilcox.test(top90_EUR$Abeta40, top90_AFR$Abeta40, paired = FALSE)
wilcox.test(top90_EUR$Abeta42, top90_AFR$Abeta42, paired = FALSE)
wilcox.test(top90_EUR$GFAP, top90_AFR$GFAP, paired = FALSE)
wilcox.test(top90_EUR$NFlight, top90_AFR$NFlight, paired = FALSE)

#Plot Beeswarm
#install.packages("beeswarm")
library(beeswarm)

# Bee swarm plot by group
pTau181_bee <- beeswarm(ptau181 ~ pop,
           data = combined_top90,
           pch = 19, 
           pwcol = gsex)
Abeta40_bee <-  beeswarm(Abeta40 ~ pop,
           data = combined_top90,
           pch = 19, 
           pwcol = gsex)
Abeta42_bee <-  beeswarm(Abeta42 ~ pop,
           data = combined_top90,
           pch = 19, 
           pwcol = gsex)
GFAPW8_bee <-  beeswarm(GFAP ~ pop,
           data = combined_top90,
           pch = 19, 
           pwcol = gsex)
NFl_bee <-  beeswarm(NFlight ~ pop,
           data = combined_top90,
           pch = 19, 
           pwcol = gsex)



library(ggstatsplot)
biomarker_list <- c("pTau181","Abeta40","Abeta42","GFAP","NFlight")
wilcox.ptau <- wilcox.test(ptau181 ~ pop, data = combined_top90)
W.ptau <- wilcox.ptau$statistic
p_value.ptau <- wilcox.ptau$p.value
n.ptau <- length(combined_top90$ptau181)

wilcox.GFAP <- wilcox.test(GFAP ~ pop, data = combined_top90)
W.GFAP <- wilcox.GFAP$statistic
p_value.GFAP <- wilcox.GFAP$p.value
n.GFAP <- length(combined_top90$GFAP)
# Add the custom subtitle with only W, p, CI, and n
custom.GFAP <- paste("W =", W.GFAP, ", p =", round(p_value.GFAP, 2), ", n =", n.GFAP)

wilcox.NFlight <- wilcox.test(NFlight ~ pop, data = combined_top90)
W.NFlight <- wilcox.NFlight$statistic
p_value.NFlight <- wilcox.NFlight$p.value
n.NFlight <- length(combined_top90$NFlight)
# Add the custom subtitle with only W, p, CI, and n
custom.ptau <- bquote(italic(W) == .(W.ptau) * "," ~ italic(p) == .(round(p_value.ptau, 3)) * "," ~ italic(n) == .(n.ptau))
custom.GFAP <- bquote(italic(W) == .(W.GFAP) * "," ~ italic(p) == .(round(p_value.GFAP, 3)) * "," ~ italic(n) == .(n.GFAP))
custom.NFlight <- bquote(italic(W) == .(W.NFlight) * "," ~ italic(p) == .(round(p_value.NFlight, 3)) * "," ~ italic(n) == .(n.NFlight))


# plot with statistical results
wilcox_pTau <- ggbetweenstats( # independent samples
  data = combined_top90,
  x = pop,
  y = ptau181,
  plot.type = "box", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE, # remove median,
  results.subtitle = FALSE
  #ggplot.component = list(geom_point(aes(color = gsex), position = position_jitterdodge())) # Add dots colored by sex
  )+ 
  theme(
    axis.title.x = element_text(size = 14), # Adjust x-axis title font size
    axis.title.y = element_text(size = 14), # Adjust y-axis title font size
    axis.text = element_text(size = 13),     # Adjust axis text font size
    plot.subtitle = element_text(size = 13) # Adjust subtitle font size
  ) +
  labs(
    x = "Population",                      # Customize x-axis title
    y = "pTau-181 level"                   # Customize y-axis title
  )+
  ggtitle(label = NULL, subtitle = custom.ptau)
# wilcox_abeta40 <- ggbetweenstats( # independent samples
#   data = combined_top90,
#   x = pop,
#   y = Abeta40,
#   plot.type = "box", # for boxplot
#   type = "nonparametric", # for wilcoxon
#   centrality.plotting = FALSE # remove median
# )
# wilcox_abeta42 <- ggbetweenstats( # independent samples
#   data = combined_top90,
#   x = pop,
#   y = Abeta42,
#   plot.type = "box", # for boxplot
#   type = "nonparametric", # for wilcoxon
#   centrality.plotting = FALSE # remove median
# )
wilcox_GFAP <- ggbetweenstats( # independent samples
  data = combined_top90,
  x = pop,
  y = GFAP,
  plot.type = "box", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE, # remove median,
  results.subtitle = FALSE
  #ggplot.component = list(geom_point(aes(color = gsex), position = position_jitterdodge())) # Add dots colored by sex
)+ 
  theme(
    axis.title.x = element_text(size = 14), # Adjust x-axis title font size
    axis.title.y = element_text(size = 14), # Adjust y-axis title font size
    axis.text = element_text(size = 13),     # Adjust axis text font size
    plot.subtitle = element_text(size = 13) # Adjust subtitle font size
  ) +
  labs(
    x = "Population",                      # Customize x-axis title
    y = "GFAP level"                   # Customize y-axis title
  )+
  ggtitle(label = NULL, subtitle = custom.GFAP)


wilcox_NFL <- ggbetweenstats( # independent samples
  data = combined_top90,
  x = pop,
  y = NFlight,
  plot.type = "box", # for boxplot
  type = "nonparametric", # for wilcoxon
  centrality.plotting = FALSE, # remove median,
  results.subtitle = FALSE
  #ggplot.component = list(geom_point(aes(color = gsex), position = position_jitterdodge())) # Add dots colored by sex
)+ 
  theme(
    axis.title.x = element_text(size = 14), # Adjust x-axis title font size
    axis.title.y = element_text(size = 14), # Adjust y-axis title font size
    axis.text = element_text(size = 13),     # Adjust axis text font size
    plot.subtitle = element_text(size = 13) # Adjust subtitle font size
  ) +
  labs(
    x = "Population",                      # Customize x-axis title
    y = "NFL level"                   # Customize y-axis title
  )+
  ggtitle(label = NULL, subtitle = custom.NFlight)

# wilcox_plot <- wilcox_pTau + wilcox_abeta40 + wilcox_abeta42 +  wilcox_GFAP+ wilcox_NFL + 
wilcox_plot <- wilcox_pTau + wilcox_GFAP+ wilcox_NFL + 
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "Wilcoxon Test", tag_levels = 'A')&
  theme(
    plot.tag = element_text(face = "bold") # Make plot tags bold
  )
print(wilcox_plot)
# 
# pdf("C:/Users/luwan/Desktop/UGA/KY_lab/wilcox.pdf", width = 15, height = 10)
# print(wilcox_plot)
# dev.off()


#Factor population
combined_top90$pop_cat <- ifelse(combined_top90$pop == "EUR", "0", "1")
combined_top90$pop_cat <- factor(combined_top90$pop_cat,ordered = F,levels=c("0","1"))
# model1: biomarker ~ ancestry + sex + age
Model1_AFREUR <- function(df, num1, num2, fileinput){
  header <- c("sampleSize","phenotype","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","pop_rsq","sex_rsq","age_rsq")
  write.table(t(as.data.frame(header)),file=fileinput,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  for (i in num1:num2){
    sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
    fitModel <- glm(unlist(df[,i])~pop_cat+gsex+Age,data=df)
    ml_summary <- summary(fitModel)
    pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
    sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
    age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
    rsqValue <- rsq.partial(fitModel, adj = TRUE)
    pop_rsq <- rsqValue$partial.rsq[1]
    sex_rsq <- rsqValue$partial.rsq[2]
    age_rsq <- rsqValue$partial.rsq[3]
    result <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(df[i]),pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,pop_rsq,sex_rsq,age_rsq))))
    write.table(result,file=fileinput,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
    pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
  }
}
Model1binary_AFREUR <- function(df, num1, num2, fileinput){
  header <- c("sampleSize","phenotype","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","pop_rsq","sex_rsq","age_rsq")
  write.table(t(as.data.frame(header)),file=fileinput,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  for (i in num1:num2){
    sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
    fitModel <- glm(unlist(df[,i])~pop_cat+gsex+Age,data=df, family = "binomial")
    ml_summary <- summary(fitModel)
    pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
    sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
    age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
    result <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(df[i]),pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue))))
    write.table(result,file=fileinput,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
    pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
  }
}
Model1_AFREUR(combined_top90, 14, 16, "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W5W8unrel_2sided_regression.txt")
Model1binary_AFREUR(combined_top90, 17, 18, "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W5W8unrel_2sided_regression.txt")



##save df for SNP association and APOE4 status analysis
#pheno
pheno_data <- data.frame(0,merged_W5W8$IID,merged_W5W8$ptau181,merged_W5W8$Abeta40,merged_W5W8$Abeta42,merged_W5W8$GFAP,merged_W5W8$NFlight)
colnames(pheno_data) <- c("FID","IID", "ptau","ab40","ab42","gfap","nfl")
#covar
PC1_20<-read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/PCA/W5W8_PCA.eigenvec",header=FALSE,sep=" ")
covar_data <- merge(merged_W5W8[c("IID", "Age", "gsex")],PC1_20, by.x = "IID", by.y = "V2")
covar_data <- covar_data[, c("V1", setdiff(names(covar_data), "V1"))]
colnames(covar_data) <- c("FID","IID", "age","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")

write.table(pheno_data,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/biomarker_data_for_SNPassociation.txt",col.names = TRUE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
write.table(covar_data,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/covar_PC_for_SNPassociation.txt",col.names = TRUE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')

#Plot PC
PC12 <- ggplot(covar_data, aes(x = PC1, y = PC2)) + 
  geom_point()
  labs(x = "PC1", y = "PC2")
PC34 <- ggplot(covar_data, aes(x = PC3, y = PC4)) + 
  geom_point()
  labs(x = "PC3", y = "PC4")
PC56 <- ggplot(covar_data, aes(x = PC5, y = PC6)) + 
  geom_point()
  labs(x = "PC5", y = "PC6")
PC78 <- ggplot(covar_data, aes(x = PC7, y = PC8)) + 
  geom_point()
PC910 <- ggplot(covar_data, aes(x = PC9, y = PC10)) + 
  geom_point()
PC1112 <- ggplot(covar_data, aes(x = PC11, y = PC12)) + 
  geom_point()
PC1314 <- ggplot(covar_data, aes(x = PC13, y = PC14)) + 
  geom_point()
PC1516 <- ggplot(covar_data, aes(x = PC15, y = PC16)) + 
  geom_point()
PC1718 <- ggplot(covar_data, aes(x = PC17, y = PC18)) + 
  geom_point()

PC_plot <- PC12 + PC34 + PC56 + PC78 + PC910 + PC1112 + PC1314 + PC1516 + PC1718 +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "PCA Plot for Participants", tag_levels = 'A')
print(PC_plot)




### APOE4 regression
APOE4_data <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_status_pheno_data.txt",header=T,sep="\t")
APOE4_data_biomarker <- merge(merged_W5W8, APOE4_data[c("IID", "APOE4")], by.x = "IID", by.y = "IID")
header <- c("sampleSize","phenotype","APOE_coef","APOE_se","APOE_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","APOE_rsq","sex_rsq","age_rsq")
write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_regression_clean.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
for (i in 17:18){
  sampleSize <- sum(!is.na(APOE4_data_biomarker[,i]) & !is.nan(APOE4_data_biomarker[,i]))
  fitModel <- glm(unlist(APOE4_data_biomarker[,i])~APOE4+gsex+Age,data=APOE4_data_biomarker,family="binomial")
  ml_summary <- summary(fitModel)
  APOE_coef <- ml_summary$coefficients[2,1]; APOE_se <- ml_summary$coefficients[2,2]; APOE_pvalue <- ml_summary$coefficients[2,4]
  sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
  age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
  result <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(APOE4_data_biomarker[i]),APOE_coef,APOE_se,APOE_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,APOE_rsq,sex_rsq,age_rsq))))
  write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  APOE_coef="NA"; APOE_se="NA";APOE_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";APOE_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
}
header <- c("sampleSize","phenotype","APOE_coef","APOE_se","APOE_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","APOE_rsq","sex_rsq","age_rsq")
write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_regression_clean.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
for (i in 15:17){
  sampleSize <- sum(!is.na(APOE4_data_biomarker[,i]) & !is.nan(APOE4_data_biomarker[,i]))
  fitModel <- glm(unlist(APOE4_data_biomarker[,i])~APOE4+gsex+Age,data=APOE4_data_biomarker)
  ml_summary <- summary(fitModel)
  APOE_coef <- ml_summary$coefficients[2,1]; APOE_se <- ml_summary$coefficients[2,2]; APOE_pvalue <- ml_summary$coefficients[2,4]
  sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
  age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
  rsqValue <- rsq.partial(fitModel, adj = TRUE)
  APOE_rsq <- rsqValue$partial.rsq[1]
  sex_rsq <- rsqValue$partial.rsq[2]
  age_rsq <- rsqValue$partial.rsq[3]
  result <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(APOE4_data_biomarker[i]),APOE_coef,APOE_se,APOE_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,APOE_rsq,sex_rsq,age_rsq))))
  write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_regression_clean.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  APOE_coef="NA"; APOE_se="NA";APOE_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";APOE_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
}






# #####try age*AFR% interaction term (expect ptau to be significant)
# header <- c("sampleSize","phenotype","APOE_coef","APOE_se","APOE_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","APOE_rsq","sex_rsq","age_rsq")
# write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
# for (i in 17:18){
#   sampleSize <- sum(!is.na(APOE4_data_biomarker[,i]) & !is.nan(APOE4_data_biomarker[,i]))
#   fitModel <- glm(unlist(APOE4_data_biomarker[,i])~APOE4+gsex+Age,data=APOE4_data_biomarker,family="binomial")
#   ml_summary <- summary(fitModel)
#   APOE_coef <- ml_summary$coefficients[2,1]; APOE_se <- ml_summary$coefficients[2,2]; APOE_pvalue <- ml_summary$coefficients[2,4]
#   sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
#   age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
#   result <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(APOE4_data_biomarker[i]),APOE_coef,APOE_se,APOE_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,APOE_rsq,sex_rsq,age_rsq))))
#   write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
#   APOE_coef="NA"; APOE_se="NA";APOE_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";APOE_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
# }
# 
# header <- c("sampleSize","phenotype","population","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","axa_coeff", "axa_se", "axa_pvalue", "pop_rsq","sex_rsq","age_rsq", "axa_rsq")
# write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/AFR_age_interaction_model.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
# for (i in 14:16){
#   sampleSize <- sum(!is.na(merged_W5W8[,i]) & !is.nan(merged_W5W8[,i]))
#   fitModel <- glm(unlist(merged_W5W8[,i])~AFR+gsex+Age+AFR*Age,data=merged_W5W8)
#   ml_summary <- summary(fitModel)
#   pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
#   sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
#   age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
#   axa_coef <- ml_summary$coefficients[5,1]; axa_se <- ml_summary$coefficients[5,2]; axa_pvalue <- ml_summary$coefficients[5,4]
#   rsqValue <- rsq.partial(fitModel, adj = TRUE)
#   pop_rsq <- rsqValue$partial.rsq[1]
#   sex_rsq <- rsqValue$partial.rsq[2]
#   age_rsq <- rsqValue$partial.rsq[3]
#   axa_rsq <- rsqValue$partial.rsq[4]
#   result <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(merged_W5W8[i]),"AFR",pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,axa_coef,axa_se,axa_pvalue,pop_rsq,sex_rsq,age_rsq,axa_rsq))))
#   write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/AFR_age_interaction_model.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
# }
# 
# 
# 
# 
# #try extreme interaction term
# #separate populations
# combined_sorted_1stquintile <- merged_W5W8[merged_W5W8$AFR < 0.75, ]  # Sorts the DataFrame in descending order based on 'target_column'
# combined_sorted_5thquintile <- merged_W5W8[merged_W5W8$AFR > 0.887, ]
# #combine the two groups and do extreme analysis
# combined_sorted_1stquintile$pop <- "Q1"  # Sorts the DataFrame in descending order based on 'target_column'
# combined_sorted_5thquintile$pop <- "Q5"
# combined_sorted_quintile <- rbind(combined_sorted_1stquintile,combined_sorted_5thquintile)
# combined_sorted_quintile$pop_cat <- ifelse(combined_sorted_quintile$pop == "Q1", "0", "1")
# combined_sorted_quintile$pop_cat <- factor(combined_sorted_quintile$pop_cat,ordered = F,levels=c("0","1"))
# #extrame regression
# df=combined_sorted_quintile
# i=9
# sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
# fitModel <- glm(unlist(df[,i])~pop_cat+gsex+Age+pop_cat*Age,data=df)
# ml_summary <- summary(fitModel)
# pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
# sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
# age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
# axa_coef <- ml_summary$coefficients[5,1]; axa_se <- ml_summary$coefficients[5,2]; axa_pvalue <- ml_summary$coefficients[5,4]
# rsqValue <- rsq.partial(fitModel, adj = TRUE)
# pop_rsq <- rsqValue$partial.rsq[1]
# sex_rsq <- rsqValue$partial.rsq[2]
# age_rsq <- rsqValue$partial.rsq[3]
# axa_rsq <- rsqValue$partial.rsq[4]
# result2 <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(df[i]),"AFR",pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,axa_coef,axa_se,axa_pvalue,pop_rsq,sex_rsq,age_rsq,axa_rsq))))
# colnames(result2) <- c("sampleSize","phenotype","population","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","axa_coeff", "axa_se", "axa_pvalue", "pop_rsq","sex_rsq","age_rsq", "axa_rsq")
# 
# #Try age in different groups
# df=combined_sorted_1stquintile
# sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
# fitModel <- glm(unlist(df[,i])~gsex+Age,data=df)
# ml_summary <- summary(fitModel)
# sex_coef <- ml_summary$coefficients[2,1]; sex_se <- ml_summary$coefficients[2,2]; sex_pvalue <- ml_summary$coefficients[2,4]
# age_coef <- ml_summary$coefficients[3,1]; age_se <- ml_summary$coefficients[3,2]; age_pvalue <- ml_summary$coefficients[3,4]
# rsqValue <- rsq.partial(fitModel, adj = TRUE)
# sex_rsq <- rsqValue$partial.rsq[1]
# age_rsq <- rsqValue$partial.rsq[2]
# result3 <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(df[i]),sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,sex_rsq,age_rsq))))
# colnames(result3) <- c("sampleSize","phenotype","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","sex_rsq","age_rsq")
# 
# df=combined_sorted_5thquintile
# sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
# fitModel <- glm(unlist(df[,i])~gsex+Age,data=df)
# ml_summary <- summary(fitModel)
# sex_coef <- ml_summary$coefficients[2,1]; sex_se <- ml_summary$coefficients[2,2]; sex_pvalue <- ml_summary$coefficients[2,4]
# age_coef <- ml_summary$coefficients[3,1]; age_se <- ml_summary$coefficients[3,2]; age_pvalue <- ml_summary$coefficients[3,4]
# rsqValue <- rsq.partial(fitModel, adj = TRUE)
# sex_rsq <- rsqValue$partial.rsq[1]
# age_rsq <- rsqValue$partial.rsq[2]
# result4 <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(df[i]),sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,sex_rsq,age_rsq))))
# colnames(result4) <- c("sampleSize","phenotype","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","sex_rsq","age_rsq")
# combined_quitile1st_5th_result <- rbind(result3,result4)
# 
# 
# #check cutoffs
# check_cutoff_for_age_beta <- function(cutoff){
#   combined_sorted_1stquintile <- merged_W5W8[merged_W5W8$AFR < cutoff, ]  # Sorts the DataFrame in descending order based on 'target_column'
#   combined_sorted_5thquintile <- merged_W5W8[merged_W5W8$AFR > cutoff, ]
#   #combine the two groups and do extreme analysis
#   #model
#   sampleSize <- sum(!is.na(combined_sorted_1stquintile[,9]) & !is.nan(combined_sorted_1stquintile[,9]))
#   fitModel <- glm(unlist(combined_sorted_1stquintile[,9])~gsex+Age,data=combined_sorted_1stquintile)
#   ml_summary <- summary(fitModel)
#   sex_coef <- ml_summary$coefficients[2,1]; sex_se <- ml_summary$coefficients[2,2]; sex_pvalue <- ml_summary$coefficients[2,4]
#   age_coef <- ml_summary$coefficients[3,1]; age_se <- ml_summary$coefficients[3,2]; age_pvalue <- ml_summary$coefficients[3,4]
#   rsqValue <- rsq.partial(fitModel, adj = TRUE)
#   sex_rsq <- rsqValue$partial.rsq[1]
#   age_rsq <- rsqValue$partial.rsq[2]
#   result1 <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(combined_sorted_1stquintile[9]),cutoff,age_coef,age_se,age_pvalue,age_rsq))))
#   colnames(result1) <- c("sampleSize","phenotype","cutoff","age_coef","age_se","age_pvalue","age_rsq")
#   
#   sampleSize <- sum(!is.na(combined_sorted_5thquintile[,9]) & !is.nan(combined_sorted_5thquintile[,9]))
#   fitModel <- glm(unlist(combined_sorted_5thquintile[,9])~gsex+Age,data=combined_sorted_5thquintile)
#   ml_summary <- summary(fitModel)
#   sex_coef <- ml_summary$coefficients[2,1]; sex_se <- ml_summary$coefficients[2,2]; sex_pvalue <- ml_summary$coefficients[2,4]
#   age_coef <- ml_summary$coefficients[3,1]; age_se <- ml_summary$coefficients[3,2]; age_pvalue <- ml_summary$coefficients[3,4]
#   rsqValue <- rsq.partial(fitModel, adj = TRUE)
#   sex_rsq <- rsqValue$partial.rsq[1]
#   age_rsq <- rsqValue$partial.rsq[2]
#   result2 <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(combined_sorted_5thquintile[9]),cutoff,age_coef,age_se,age_pvalue,age_rsq))))
#   colnames(result2) <- c("sampleSize","phenotype","cutoff","age_coef","age_se","age_pvalue","age_rsq")
#   result<-rbind(result1,result2)
#   return(result)
# }
# 
# cutoffs <- c(0.5,0.6,0.7,0.8,0.9)
# for (j in cutoffs){
#   result_name <- paste0("result", gsub("\\.", "", as.character(j)))  # Create a name like result05, result06, etc.
#   assign(result_name, data.frame(check_cutoff_for_age_beta(j)))  # Assign the dataframe to the dynamically created name
# }
# result05<-data.frame(check_cutoff_for_age_beta(0.5))
# result06<-data.frame(check_cutoff_for_age_beta(0.6))
# result07<-data.frame(check_cutoff_for_age_beta(0.7))
# result08<-data.frame(check_cutoff_for_age_beta(0.8))
# result09<-data.frame(check_cutoff_for_age_beta(0.9))
# 
# result00 <- rbind(result05,result06,result07,result08,result09)


##PGS
pTauPGS <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/PGS/pgs_ptauvariants_W5W8.profile",header=T,sep="")
abeta40PGS <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/PGS/pgs_abeta40variants_W5W8.profile",header=T,sep="")
abeta42PGS <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/PGS/pgs_abeta42variants_W5W8.profile",header=T,sep="")
GFAPPGS <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/PGS/pgs_GFAPvariants_W5W8.profile",header=T,sep="")
NFLPGS <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/PGS/pgs_NFLvariants_W5W8.profile",header=T,sep="")
colnames(pTauPGS)[colnames(pTauPGS) == "SCORESUM"] <- "pTau_PGS"
colnames(abeta40PGS)[colnames(abeta40PGS) == "SCORESUM"] <- "abeta40_PGS"
colnames(abeta42PGS)[colnames(abeta42PGS) == "SCORESUM"] <- "abeta42_PGS"
colnames(GFAPPGS)[colnames(GFAPPGS) == "SCORESUM"] <- "GFAP_PGS"
colnames(NFLPGS)[colnames(NFLPGS) == "SCORESUM"] <- "NFL_PGS"

#conbine
df_combo_pgs <- merge(merged_W5W8, pTauPGS[c("IID", "pTau_PGS")], by.x = "IID", by.y = "IID")
df_combo_pgs <- merge(df_combo_pgs, abeta40PGS[c("IID", "abeta40_PGS")], by.x = "IID", by.y = "IID")
df_combo_pgs <- merge(df_combo_pgs, abeta42PGS[c("IID", "abeta42_PGS")], by.x = "IID", by.y = "IID")
df_combo_pgs <- merge(df_combo_pgs, GFAPPGS[c("IID", "GFAP_PGS")], by.x = "IID", by.y = "IID")
df_combo_pgs <- merge(df_combo_pgs, NFLPGS[c("IID", "NFL_PGS")], by.x = "IID", by.y = "IID")


#Check biomarker distribution
pTau181_W8_PGS <- ggplot(df_combo_pgs, aes(x = pTau_PGS, y = ptau181)) + 
  geom_point()+
  labs(x = "PGS", y = "pTau181")
# Abeta40W8_PGS <- ggplot(df_combo_pgs, aes(x = abeta40_PGS, y = Abeta40W8)) + 
#   geom_point()+
#   labs(x = "PGS", y = "Abeta40")
# Abeta42W8_PGS <- ggplot(df_combo_pgs, aes(x = abeta42_PGS, y = Abeta42W8)) + 
#   geom_point()+
#   labs(x = "PGS", y = "Abeta42")
GFAPW8_PGS <- ggplot(df_combo_pgs, aes(x = GFAP_PGS, y = GFAP)) + 
  geom_point()+
  labs(x = "PGS", y = "GFAP")
NFlightW8_PGS <- ggplot(df_combo_pgs, aes(x = NFL_PGS, y = NFlight)) + 
  geom_point()+
  labs(x = "PGS", y = "NFL")

# PGS_vs_biomarker <- pTau181_W8_PGS + Abeta40W8_PGS + Abeta42W8_PGS + GFAPW8_PGS + NFlightW8_PGS +
PGS_vs_biomarker <- pTau181_W8_PGS + GFAPW8_PGS + NFlightW8_PGS +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "PGS vs Biomarkers", tag_levels = 'A')
print(PGS_vs_biomarker)


#correlation
pTau_PGS_correlation <- cor.test(df_combo_pgs$pTau_PGS, df_combo_pgs$ptau181, method = "pearson", use = "complete.obs")
# abeta40_PGS_correlation <- cor.test(df_combo_pgs$abeta40_PGS, df_combo_pgs$Abeta40W8, method = "pearson", use = "complete.obs")
# abeta42_PGS_correlation <- cor.test(df_combo_pgs$abeta42_PGS, df_combo_pgs$Abeta42W8, method = "pearson", use = "complete.obs")
GFAP_PGS_correlation <- cor.test(df_combo_pgs$GFAP_PGS, df_combo_pgs$GFAP, method = "pearson", use = "complete.obs")
NFL_PGS_correlation <- cor.test(df_combo_pgs$NFL_PGS, df_combo_pgs$NFlight, method = "pearson", use = "complete.obs")
print(pTau_PGS_correlation)
# print(abeta40_PGS_correlation)
# print(abeta42_PGS_correlation)
print(GFAP_PGS_correlation)
print(NFL_PGS_correlation)











#Plot gender difference
ggplot(merged_W5W8, aes(x = gsex, y = log_GFAP, fill = gsex)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Biomarker Levels by Gender",
       x = "Gender",
       y = "Biomarker Level") +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal()

#Plot gender difference
library(ggsignif)
pTausex <- ggplot(merged_W5W8, aes(x = gsex, y = log_pTau181, fill = gsex)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(x = "Biological sex",
       y = "log(p-Tau181)") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = c("Women", "Men")) + 
  theme_minimal()+
  theme(legend.position = "none",
        text = element_text(size = 12),  # Changes global text size
        axis.title = element_text(size = 14),  # Specific size for axis titles
        axis.text = element_text(size = 12))  # Specific size for axis text
buffer1 = max(merged_W5W8$log_pTau181, na.rm = TRUE) * 0.08  # 5% above the max value
ptausex_sig <- pTausex + geom_signif(comparisons = list(c("0", "1")), annotations="*", y_position = max(merged_W5W8$log_pTau181, na.rm = TRUE) + buffer1, tip_length = 0.02)

GFAPsex <- ggplot(merged_W5W8, aes(x = gsex, y = log_GFAP, fill = gsex)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(x = "Biological sex",
       y = "log(GFAP)") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = c("Women", "Men")) + 
  theme_minimal()+
  theme(legend.position = "none",
        text = element_text(size = 12),  # Changes global text size
        axis.title = element_text(size = 14),  # Specific size for axis titles
        axis.text = element_text(size = 12))  # Specific size for axis text
buffer = max(merged_W5W8$log_GFAP, na.rm = TRUE) * 0.06  # 5% above the max value
GFAPsex_sig <- GFAPsex + geom_signif(comparisons = list(c("0", "1")), annotations="**", y_position = max(merged_W5W8$log_GFAP, na.rm = TRUE) + buffer, tip_length = 0.02)
NFLsex <- ggplot(merged_W5W8, aes(x = gsex, y = log_NFlight, fill = gsex)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(x = "Biological sex",
       y = "log(NFL)") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_x_discrete(labels = c("Women", "Men")) + 
  theme_minimal()+
  theme(legend.position = "none",
        text = element_text(size = 12),  # Changes global text size
        axis.title = element_text(size = 14),  # Specific size for axis titles
        axis.text = element_text(size = 12))  # Specific size for axis text
#install.packages("ggbeeswarm")
# library(ggbeeswarm)
# ggplot(merged_W5W8, aes(x = gsex, y = log_GFAP, color = gsex)) +
#   geom_quasirandom(groupOnX = TRUE) +  # This creates the beeswarm effect
#   labs(title = "Distribution of Biomarker Levels by Gender",
#        x = "Gender",
#        y = "Biomarker Level") +
#   theme_minimal() +
#   scale_color_brewer(palette = "Set1")
# PGS_vs_biomarker <- pTau181_W8_PGS + Abeta40W8_PGS + Abeta42W8_PGS + GFAPW8_PGS + NFlightW8_PGS +
sex_vs_biomarker <- ptausex_sig + GFAPsex_sig + NFLsex +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(
    title = "Biomarker Levels by Biological Sex",
    tag_levels = 'A',
    theme = theme(plot.title = element_text(size = 20)) # Set the title text size here
  )
print(sex_vs_biomarker)

