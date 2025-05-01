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
W5_W8_diff_rel <- merge(df_combo, biomarker_w5[c("PID", "pTau181_W5", "Abeta40W5","Abeta42W5","GFAPW5","NFlightW5")], by.x = "PID", by.y = "PID")
#unrel_list <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/combined_unrel_AFR0.9_list.txt",header=F, sep="\t")
unrel_list <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/all_keep_indiv_unrel.txt",header=F, sep=" ")



##subtraction method
W5_W8_diff <- W5_W8_diff_rel[W5_W8_diff_rel$IID %in% unrel_list$V2, ]
W5_W8_diff$pTau181_diff = W5_W8_diff$pTau181_W8 - W5_W8_diff$pTau181_W5
W5_W8_diff$GFAP_diff = W5_W8_diff$GFAPW8 - W5_W8_diff$GFAPW5
W5_W8_diff$NFlight_diff = W5_W8_diff$NFlightW8 - W5_W8_diff$NFlightW5

#Check phenotype
library(ggplot2)
pTau181_W85_hist <- ggplot(W5_W8_diff, aes(x = pTau181_diff)) + 
  geom_histogram(binwidth = 1, fill = "black", color = "black") +
  labs(x = "pTau181", y = "# of Individuals")
Abeta40W85_hist <- ggplot(W5_W8_diff, aes(x = Abeta40_diff)) + 
  geom_histogram(binwidth = 1, fill = "black", color = "black") +
  labs(x = "Abeta40", y = "# of Individuals")
Abeta42W85_hist <- ggplot(W5_W8_diff, aes(x = Abeta42_diff)) + 
  geom_histogram(fill = "black", color = "black") +
  labs(x = "Abeta42", y = "# of Individuals")
GFAPW85_hist <- ggplot(W5_W8_diff, aes(x = GFAP_diff)) + 
  geom_histogram(binwidth = 1, fill = "black", color = "black") +
  labs(x = "GFAP", y = "# of Individuals")
NFlightW85_hist <- ggplot(W5_W8_diff, aes(x = NFlight_diff)) + 
  geom_histogram(binwidth = 1, fill = "black", color = "black") +
  labs(x = "NFL", y = "# of Individuals")
#install.packages("patchwork")
library(patchwork)
W85_hist_plot <- pTau181_W85_hist + GFAPW85_hist + NFlightW85_hist +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "W8 - W5", tag_levels = 'A')
print(W85_hist_plot)

#biomarker diff vs AFR
check_diff_AFR <- function(df,pheno,lab_pheno="pheno_name"){
  diff_AFR <-
    ggplot(df, aes(x = pop2, y = pheno)) + 
    geom_point()+
    labs(x = "AFR ancestry proportion", y = lab_pheno)+
    theme(text = element_text(size = 16),  # Adjusts overall text size
          axis.title = element_text(size = 14),  # Specific size for axis titles
          axis.text = element_text(size = 12))  # Specific size for axis text
  return(diff_AFR)}
pTau181_diff_AFR <- check_diff_AFR(W5_W8_diff,W5_W8_diff$pTau181_diff,lab_pheno="\u0394(p-Tau181)")
Abeta40_dif_AFR <- check_diff_AFR(W5_W8_diff,W5_W8_diff$Abeta40_diff,lab_pheno="Abeta40")
Abeta42_dif_AFR <- check_diff_AFR(W5_W8_diff,W5_W8_diff$Abeta42_diff,lab_pheno="Abeta42")
GFAP_dif_AFR <- check_diff_AFR(W5_W8_diff,W5_W8_diff$GFAP_diff,lab_pheno="\u0394(GFAP)")
NFlight_dif_AFR <- check_diff_AFR(W5_W8_diff,W5_W8_diff$NFlight_diff,lab_pheno="\u0394(NFL)")


AFR_vs_biomarker_diff <- pTau181_diff_AFR + GFAP_dif_AFR + NFlight_dif_AFR +
  plot_layout(guides = 'collect') + # Optional: collect all legends into one
  plot_annotation(title = "AFR Ancestry Proportion vs Difference in Biomarker", tag_levels = 'A',
                  theme = theme(plot.title = element_text(size = 20)))  # Increase and bold the title siz)

print(AFR_vs_biomarker_diff)







# model1: biomarker ~ ancestry% + sex + age
Model1 <- function(df, values, fileinput){
  header <- c("sampleSize","phenotype","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","pop_rsq","sex_rsq","age_rsq")
  write.table(t(as.data.frame(header)),file=fileinput,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  for (i in values){
    sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
    fitModel <- glm(unlist(df[,i])~pop2+gsex+AgePCSCw8,data=df)
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

Model1binary <- function(df, values, fileinput){
  header <- c("sampleSize","phenotype","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","pop_rsq","sex_rsq","age_rsq")
  write.table(t(as.data.frame(header)),file=fileinput,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  for (i in values){
    sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
    fitModel <- glm(unlist(df[,i])~pop2+gsex+AgePCSCw8,data=df, family = "binomial")
    ml_summary <- summary(fitModel)
    pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
    sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
    age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
    result <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(df[i]),pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue))))
    write.table(result,file=fileinput,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
    pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
  }
}
Model1(W5_W8_diff, c(24,27,28), "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W8-5unrel_diff_regression.txt")
Model1binary(W5_W8_diff,c(29,30), "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W8-5unrel_diff_regression.txt")


# #Try AFR>0.9
# AA0.9_W5_W8_diff <-W5_W8_diff[W5_W8_diff$pop2 >= 0.1, ]
# # model1: biomarker ~ ancestry% + sex + age
# header <- c("sampleSize","phenotype","population","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","pop_rsq","sex_rsq","age_rsq")
# write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W8-5unrel_diff_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
# for (i in 19:23){
#   sampleSize <- sum(!is.na(AA0.9_W5_W8_diff[,i]))
#   fitModel <- glm(unlist(AA0.9_W5_W8_diff[,i])~pop2+gsex+AgePCSCw8,data=AA0.9_W5_W8_diff)#,family = "binomial")
#   ml_summary <- summary(fitModel)
#   pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
#   sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
#   age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
#   rsqValue <- rsq.partial(fitModel, adj = TRUE)
#   pop_rsq <- rsqValue$partial.rsq[1]
#   sex_rsq <- rsqValue$partial.rsq[2]
#   age_rsq <- rsqValue$partial.rsq[3]
#   result <- t(as.data.frame(c(sampleSize,colnames(AA0.9_W5_W8_diff[i]),"AFR0.9",pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,pop_rsq,sex_rsq,age_rsq)))
#   write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W8-5unrel_diff_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
#   pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
# }


#Remove missing columns

W5_W8_diff_clean <- W5_W8_diff %>% filter(!is.na(GFAP_diff))
##add race and count individuals
W5_W8_diff_clean <- merge(W5_W8_diff_clean, sexinfo[c("PID", "Race")], by.x = "PID", by.y = "PID")
table(W5_W8_diff_clean$Race)



# model1: biomarker ~ ancestry% + sex + age
Model1(W5_W8_diff, c(39,42,43), "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W8-5unrel_diff_regression.txt")
Model1binary(W5_W8_diff,c(44,45), "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W8-5unrel_diff_regression.txt")



##check ptau increase(1) and no increase(0)
W5_W8_diff$pTau181_diff_cat <- ifelse(W5_W8_diff$pTau181_diff < 0.00001, "0", "1")
W5_W8_diff$pTau181_diff_cat <- factor(W5_W8_diff$pTau181_diff_cat,ordered = F,levels=c("0","1"))

Model1binary(W5_W8_diff,46, "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W8-5unrel_ptau_0vs1_regression.txt")




##plot boxplot based on quintiles of ptau ancestry
df_quintile <- data.frame(W5_W8_diff$pop2,W5_W8_diff$pTau181_diff)
colnames(df_quintile) <- c("AFR","ptau_diff")
df_quintile <- df_quintile %>% filter(!is.na(ptau_diff))
df_quintile <- df_quintile %>%
  mutate(quintile = ntile(AFR, 5))  # Create a new column for quintiles

# Display the plot
quintile_cutoffs <- quantile(df_quintile$AFR, probs = seq(0, 1, by = 0.2))
print(quintile_cutoffs)

quintile_labels <- sapply(2:length(quintile_cutoffs), function(i) {
  paste0("Q", i-1, ": [", round(quintile_cutoffs[i-1], 2), "-", round(quintile_cutoffs[i], 2), "]")
})

#Plot bar plot of means and standard error
df_quintile <- df_quintile %>%
  mutate(quintile_label = quintile_labels[quintile])
# Plot boxplots with new quintile range labels
quintile_plot <- ggplot(df_quintile, aes(x = quintile_label, y = ptau_diff)) +
  geom_boxplot() +
  labs(x = "Quintile of AFR Proportion", y = "pTau_diff") +
  theme_minimal()

# Calculate mean and standard error
stats <- df_quintile %>%
  group_by(quintile_label) %>%
  summarise(
    mean = mean(ptau_diff),
    se = sd(ptau_diff) / sqrt(n()),
    .groups = 'drop'
  )
stats$quintile_label <- factor(as.character(stats$quintile_label),ordered = F,levels = c("Q1: [0.01-0.75]","Q2: [0.75-0.81]","Q3: [0.81-0.85]","Q4: [0.85-0.89]","Q5: [0.89-1]"))

bar_colors <- as.vector(c("#0073C2FF","#0073C2FF","#EFC000FF","#EFC000FF","#EFC000FF")) #see pal_jco

quintile_plot1 <- ggplot(stats, aes(x = quintile_label, y = mean, fill = quintile_label)) +
  geom_col(width = 0.7) +  # This creates the bar plot
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1,linewidth = 0.2) +  # Adds error bars
  labs(x = "Quintiles of AFR Proportion", y = "Average change in pTau-181") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)      # Change axis text size
  )+
  scale_fill_manual(values = as.vector(bar_colors))
#scale_fill_brewer(palette = "Paired", guide = FALSE)  # No need for a legend here
quintile_plot1 <- quintile_plot1 + theme(legend.position = "none")
# Display the plot
print(quintile_plot1)



## Plot ancestry proportion in categories (ptau increase/no increase)
plot_ptau <- data.frame(W5_W8_diff$IID,W5_W8_diff$pop2,W5_W8_diff$pTau181_diff)
colnames(plot_ptau) <- c("IID","pop2","pTau181_diff")
plot_ptau <- plot_ptau %>%
  filter(!is.na(pTau181_diff))

plot_ptau$ptau_status <- ifelse(plot_ptau$pTau181_diff > 0.00001, "Increase", "No change/decrease")
plot_ptau$ptau_status <- factor(plot_ptau$ptau_status,ordered = F,levels=c("Increase","No change/decrease"))
# Calculate the counts for each group
counts <- plot_ptau %>%
  group_by(ptau_status) %>%
  summarize(n = n())
p <- ggplot(plot_ptau, aes(x = factor(ptau_status), y = pop2, fill = ptau_status)) +
  geom_boxplot(na.rm = TRUE,width = 0.6) +
  scale_fill_jco(name = "Change in pTau-181") +
  labs(x = "Change in pTau-181", y = "African Ancestry Proportion") +#, title = "Change in pTau vs. AFR ancestry proportion") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)
  ) +
  # Add text annotations for sample size
  geom_text(
    data = counts,
    aes(x = ptau_status, y = max(plot_ptau$pop2, na.rm = TRUE) + 0.02, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    size = 4
  )
print(p)

#perform wilcox test
wilcox_increase_decrease <- wilcox.test(pop2 ~ ptau_status, data = plot_ptau)
print(wilcox_increase_decrease)







# ## Plot ancestry proportion in quintiles and make proportion histogram
plot_ptau1 <- W5_W8_diff[complete.cases(W5_W8_diff["pTau181_diff"]), ]

plot_ptau1 <- plot_ptau1 %>%
  mutate(AFR_quintile = ntile(pop2, 5))

# Calculate quantile ranges separately
quantile_ranges <- data.frame(
  AFR_quintile = 1:5,
  lower = round(quantile(plot_ptau1$pop2, probs = seq(0, 0.8, by = 0.2)), 2),
  upper = round(quantile(plot_ptau1$pop2, probs = seq(0.2, 1, by = 0.2)), 2)
)



# Create labels based on quantile ranges
quantile_ranges$label <- paste("Q", quantile_ranges$AFR_quintile, ": [",
                               quantile_ranges$lower, "-", quantile_ranges$upper, "]", sep = "")

# Merge quantile ranges back to the dataframe
plot_ptau1 <- plot_ptau1 %>%
  left_join(quantile_ranges, by = "AFR_quintile") %>%
  mutate(AFR_quintile = label) %>%
  select(-label, -lower, -upper)  # Remove extra columns

# Prepare data for plotting by categorizing 'diff' and summarizing
plot_ptau1 <- plot_ptau1 %>%
  mutate(diff_cat = ifelse(pTau181_diff > 0, "Increase", "No change/decrease")) %>%
  group_by(AFR_quintile, diff_cat) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(AFR_quintile) %>%
  mutate(prop = count / sum(count))

# Plotting
library("ggsci")


p2 <- ggplot(plot_ptau1, aes(x = AFR_quintile, y = prop, fill = diff_cat)) +
  geom_bar(stat = "identity", position = position_stack(),width = 0.7) +
  labs(x = "Quintiles of AFR Proportion", y = "Individual Proportion") +#, title = "Proportions of diff Categories by AFR Quintile") +
  scale_fill_jco(name = "Change in pTau-181") + #previously have palette = "Darjeeling1", but didn't work anymore
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 14),  # Change legend title size
        legend.text = element_text(size = 12),    # Change legend text size
        axis.title = element_text(size = 15),    # Change axis title size
        axis.text = element_text(size = 12)      # Change axis text size
  )
print(p2)

# Extract the legend
legend_p2 <- ggplotGrob(p2)$grobs[[which(sapply(ggplotGrob(p2)$grobs, function(x) x$name) == "guide-box")]]
# Remove the original legend from p2
p2 <- p2 + theme(legend.position = "none")


library(gridExtra)
library(ggpubr)
grid.arrange(
  arrangeGrob(legend_p2, vp = grid::viewport(height = unit(2, "lines"))),  # The legend on top
  ggarrange(quintile_plot1, p2, p, ncol = 3, labels = c("A", "B", "C")),   # Arrange plots
  heights = c(2, 9),  # Adjust the height ratio to give less space to legend and more to plots
  layout_matrix = rbind(c(1, 1, 1),
                        c(2, 2, 2))
)

#save to pdf
pdf("C:/Users/luwan/Desktop/UGA/KY_lab/manu/alz_submission1.6/FIGURE2.pdf", width = 10, height = 5.5)  # Adjust width and height as needed

# Create the plot
grid.arrange(
  arrangeGrob(legend_p2, vp = grid::viewport(height = unit(2, "lines"))),  # The legend on top
  ggarrange(quintile_plot1, p2, p, ncol = 3, labels = c("A", "B", "C")),   # Arrange plots
  heights = c(2, 9),  # Adjust the height ratio to give less space to legend and more to plots
  layout_matrix = rbind(c(1, 1, 1),
                        c(2, 2, 2))
)

# Close the PDF device
dev.off()


# Redo with new clean dataset
Model1(W5_W8_diff_clean, c(19,22,23), "C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/RegressionResult/W8-5unrel_diff_regression.txt")
Model1(W5_W8_diff_clean, c(34,37,38), "C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/RegressionResult/W8-5unrel_diff_regression.txt")



#Sensitivity analysis
#Filter out AFR<0.1
W5_W8_diff_clean_sens  <- W5_W8_diff_clean[W5_W8_diff_clean$pop2 >= 0.3, ]
i=23
df=W5_W8_diff_clean_sens
sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
fitModel <- glm(unlist(df[,i])~pop2+gsex+AgePCSCw8,data=df)
ml_summary <- summary(fitModel)
pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
rsqValue <- rsq.partial(fitModel, adj = TRUE)
pop_rsq <- rsqValue$partial.rsq[1]
sex_rsq <- rsqValue$partial.rsq[2]
age_rsq <- rsqValue$partial.rsq[3]
result <- as.data.frame(t(as.data.frame(c(sampleSize,colnames(df[i]),pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,pop_rsq,sex_rsq,age_rsq))))



### Try PC1
# read in PC 
popinfo <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Admixture/20130606_g1k_3202_samples_ped_population.txt",header=TRUE)
PCA_1kGP_ref <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/PCA_1kGP_Gibbons/1kGP_PCA_results_ref.eigenvec")
PCA_1kGP_eigenval <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/PCA_1kGP_Gibbons/1kGP_PCA_results_ref.eigenval")
PCA_gibbons_proj <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/PCA_1kGP_Gibbons/Gibbons_PCA_1kGP_projection.sscore")
colnames(PCA_1kGP_ref) <- c("FID","SampleID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA_gibbons_proj <- PCA_gibbons_proj[, !(names(PCA_gibbons_proj) %in% c("V3", "V4"))]
colnames(PCA_gibbons_proj) <- c("FID","SampleID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA_gibbons_proj$Population <- "Cross-sectional"
PCA_gibbons_proj$Superpopulation <- "Cross-sectional"
# Keep only Combined dataset
unrel_list <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/all_keep_indiv_unrel.txt",header=F, sep=" ")
PCA_gibbons_proj <- PCA_gibbons_proj[PCA_gibbons_proj$SampleID %in% unrel_list$V2, ]

#scale PC projection
PCA_gibbons_proj$PC1=PCA_gibbons_proj$PC1/((-PCA_1kGP_eigenval[1,1]**0.5)/2)
PCA_gibbons_proj$PC2=PCA_gibbons_proj$PC2/((-PCA_1kGP_eigenval[2,1]**0.5)/2)
PCA_gibbons_proj$PC3=PCA_gibbons_proj$PC3/((-PCA_1kGP_eigenval[3,1]**0.5)/2)
PCA_gibbons_proj$PC4=PCA_gibbons_proj$PC4/((-PCA_1kGP_eigenval[4,1]**0.5)/2)
PCA_gibbons_proj$PC5=PCA_gibbons_proj$PC5/((-PCA_1kGP_eigenval[5,1]**0.5)/2)
PCA_gibbons_proj$PC6=PCA_gibbons_proj$PC6/((-PCA_1kGP_eigenval[6,1]**0.5)/2)
PCA_gibbons_proj$PC7=PCA_gibbons_proj$PC7/((-PCA_1kGP_eigenval[7,1]**0.5)/2)
PCA_gibbons_proj$PC8=PCA_gibbons_proj$PC8/((-PCA_1kGP_eigenval[8,1]**0.5)/2)
PCA_gibbons_proj$PC9=PCA_gibbons_proj$PC9/((-PCA_1kGP_eigenval[9,1]**0.5)/2)
PCA_gibbons_proj$PC10=PCA_gibbons_proj$PC10/((-PCA_1kGP_eigenval[10,1]**0.5)/2)

#attach to biomarkers
colnames(PCA_gibbons_proj)[colnames(PCA_gibbons_proj) == "SampleID"] <- "IID"
W5_W8_diff_pca <- merge(W5_W8_diff, PCA_gibbons_proj[c("IID", "PC1")], by.x = "IID", by.y = "IID")
# model1: biomarker ~ ancestry% + sex + age
Model_PCA <- function(df, values, fileinput){
  header <- c("sampleSize","phenotype","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","pop_rsq","sex_rsq","age_rsq")
  write.table(t(as.data.frame(header)),file=fileinput,col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  for (i in values){
    sampleSize <- sum(!is.na(df[,i]) & !is.nan(df[,i]))
    fitModel <- glm(unlist(df[,i])~PC1+gsex+AgePCSCw8,data=df)
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
Model_PCA(W5_W8_diff_pca, c(19,22,23), "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/W8-5unrel_diff_regression_PCA.txt")




### Test APOE4 vs change
APOE4_data <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_status_pheno_data.txt",header=T,sep="\t")
APOE4_biomarker_diff <- merge(W5_W8_diff, APOE4_data[c("IID", "APOE4")], by.x = "IID", by.y = "IID")

#Read in PCA data
# PCA_ID <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/GT_for_PCA/W5W8_indiv_unrel_AFR0.9.fam",header=F,sep=" ")
PCA_data <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/PCA/W5W8_PCA.eigenvec",header=F,sep=" ")
colnames(PCA_data) <- c("FID", "IID", "sex","age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
APOE4_biomarker_diff <- merge(APOE4_biomarker_diff, PCA_data[c("IID", "PC1","PC2","PC3","PC4")], by.x = "IID", by.y = "IID")

i=23
# header <- c("sampleSize","phenotype","APOE_coef","APOE_se","APOE_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","APOE_rsq","sex_rsq","age_rsq")
# write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_regression_clean.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
# for (i in 15:17){
sampleSize <- sum(!is.na(APOE4_biomarker_diff[,i]) & !is.nan(APOE4_biomarker_diff[,i]))
fitModel <- glm(unlist(APOE4_biomarker_diff[,i])~APOE4+gsex+AgePCSCw8+PC1+PC2+PC3+PC4,data=APOE4_biomarker_diff)
ml_summary <- summary(fitModel)
APOE_coef <- ml_summary$coefficients[2,1]; APOE_se <- ml_summary$coefficients[2,2]; APOE_pvalue <- ml_summary$coefficients[2,4]
sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
rsqValue <- rsq.partial(fitModel, adj = TRUE)
APOE_rsq <- rsqValue$partial.rsq[1]
sex_rsq <- rsqValue$partial.rsq[2]
age_rsq <- rsqValue$partial.rsq[3]
result <- (as.data.frame(c(sampleSize,colnames(APOE4_biomarker_diff[i]),APOE_coef,APOE_se,APOE_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,APOE_rsq,sex_rsq,age_rsq)))
write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/APOE4/APOE4_regression_clean.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
APOE_coef="NA"; APOE_se="NA";APOE_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";APOE_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
# }





##############-------------------------------------
####Try divide APOE4 carrier (>0) and non-carrier and redo regression
APOE4_biomarker_diff_0 <- APOE4_biomarker_diff %>% filter(APOE4 == 0)
APOE4_biomarker_diff_carrier <- APOE4_biomarker_diff %>% filter(APOE4 > 0)
Model1(APOE4_biomarker_diff_0, c(19,22,23), "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/APOE4carrier_regression.txt")
Model1(APOE4_biomarker_diff_carrier, c(19,22,23), "C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/APOE4carrier_regression.txt")

#Interaction term
##Interaction term
header <- c("sampleSize","phenotype","population","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","APOE4_coef","APOE4_se","APOE4_pvalue","AFRxAPOE4_coef","AFRxAPOE4_se","AFRxAPOE4_pvalue","pop_rsq","sex_rsq","age_rsq","APOE4_rsq","AFRxAPOE4_rsq")
write.table(t(as.data.frame(header)),file="C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/APOE4interaction_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
numbers <- c(19, 22, 23)
for (i in numbers){
  sampleSize <- sum(!is.na(APOE4_biomarker_diff[,i]))
  fitModel <- glm(unlist(APOE4_biomarker_diff[,i])~pop2+gsex+AgePCSCw8+APOE4+pop2*APOE4,data=APOE4_biomarker_diff)#,family = "binomial")
  ml_summary <- summary(fitModel)
  pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
  sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
  age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
  APOE4_coef <- ml_summary$coefficients[5,1]; APOE4_se <- ml_summary$coefficients[5,2]; APOE4_pvalue <- ml_summary$coefficients[5,4]
  AFRxAPOE4_coef <- ml_summary$coefficients[6,1]; AFRxAPOE4_se <- ml_summary$coefficients[6,2]; AFRxAPOE4_pvalue <- ml_summary$coefficients[6,4]
  rsqValue <- rsq.partial(fitModel, adj = TRUE)
  pop_rsq <- rsqValue$partial.rsq[1]
  sex_rsq <- rsqValue$partial.rsq[2]
  age_rsq <- rsqValue$partial.rsq[3]
  APOE4_rsq <- rsqValue$partial.rsq[4]
  AFRxAPOE4_rsq <- rsqValue$partial.rsq[5]
  result <- t(as.data.frame(c(sampleSize,colnames(APOE4_biomarker_diff[i]),"AFR",pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,APOE4_coef,APOE4_se,APOE4_pvalue,AFRxAPOE4_coef,AFRxAPOE4_se,AFRxAPOE4_pvalue,pop_rsq,sex_rsq,age_rsq,APOE4_rsq,AFRxAPOE4_rsq)))
  write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/Regression/results/APOE4interaction_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
  pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
}

##compare ANOVA interaction term
APOE4_ptau_diff_anova <- lm(pTau181_diff ~ pop2 + APOE4 + gsex + AgePCSCw8,data = APOE4_biomarker_diff)
APOE4_ptau_diff_ixn_anova <- lm(pTau181_diff ~ pop2 + APOE4 + gsex + AgePCSCw8 + pop2*APOE4,data = APOE4_biomarker_diff)
summary(APOE4_ptau_ixn_anova)
anova(APOE4_ptau_diff_anova,APOE4_ptau_diff_ixn_anova)






###-------------------------------####
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
diff_pgs <- merge(W5_W8_diff, pTauPGS[c("IID", "pTau_PGS")], by.x = "IID", by.y = "IID")
diff_pgs <- merge(diff_pgs, abeta40PGS[c("IID", "abeta40_PGS")], by.x = "IID", by.y = "IID")
diff_pgs <- merge(diff_pgs, abeta42PGS[c("IID", "abeta42_PGS")], by.x = "IID", by.y = "IID")
diff_pgs <- merge(diff_pgs, GFAPPGS[c("IID", "GFAP_PGS")], by.x = "IID", by.y = "IID")
diff_pgs <- merge(diff_pgs, NFLPGS[c("IID", "NFL_PGS")], by.x = "IID", by.y = "IID")


sampleSize <- sum(!is.na(diff_pgs[,19]))
fitModel <- glm(unlist(diff_pgs[,19])~pTau_PGS+gsex+AgePCSCw8,data=diff_pgs)#,family = "binomial")
ml_summary <- summary(fitModel)
PGS_coef <- ml_summary$coefficients[2,1]; PGS_se <- ml_summary$coefficients[2,2]; PGS_pvalue <- ml_summary$coefficients[2,4]
sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
rsqValue <- rsq.partial(fitModel, adj = TRUE)
PGS_rsq <- rsqValue$partial.rsq[1]
sex_rsq <- rsqValue$partial.rsq[2]
age_rsq <- rsqValue$partial.rsq[3]
result <- t(as.data.frame(c(sampleSize,colnames(diff_pgs[19]),PGS_coef,PGS_se,PGS_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,PGS_rsq,sex_rsq,age_rsq)))
write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/PGS/PGS_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
PGS_coef="NA"; PGS_se="NA";PGS_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";PGS_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"

#GFAP
sampleSize <- sum(!is.na(diff_pgs[,22]))
fitModel <- glm(unlist(diff_pgs[,22])~GFAP_PGS+gsex+AgePCSCw8,data=diff_pgs)#,family = "binomial")
ml_summary <- summary(fitModel)
PGS_coef <- ml_summary$coefficients[2,1]; PGS_se <- ml_summary$coefficients[2,2]; PGS_pvalue <- ml_summary$coefficients[2,4]
sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
rsqValue <- rsq.partial(fitModel, adj = TRUE)
PGS_rsq <- rsqValue$partial.rsq[1]
sex_rsq <- rsqValue$partial.rsq[2]
age_rsq <- rsqValue$partial.rsq[3]
result <- t(as.data.frame(c(sampleSize,colnames(diff_pgs[22]),PGS_coef,PGS_se,PGS_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,PGS_rsq,sex_rsq,age_rsq)))
write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/PGS/PGS_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
PGS_coef="NA"; PGS_se="NA";PGS_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";PGS_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"

#NFL
sampleSize <- sum(!is.na(diff_pgs[,23]))
fitModel <- glm(unlist(diff_pgs[,23])~NFL_PGS+gsex+AgePCSCw8,data=diff_pgs)#,family = "binomial")
ml_summary <- summary(fitModel)
PGS_coef <- ml_summary$coefficients[2,1]; PGS_se <- ml_summary$coefficients[2,2]; PGS_pvalue <- ml_summary$coefficients[2,4]
sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
rsqValue <- rsq.partial(fitModel, adj = TRUE)
PGS_rsq <- rsqValue$partial.rsq[1]
sex_rsq <- rsqValue$partial.rsq[2]
age_rsq <- rsqValue$partial.rsq[3]
result <- t(as.data.frame(c(sampleSize,colnames(diff_pgs[16]),PGS_coef,PGS_se,PGS_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,PGS_rsq,sex_rsq,age_rsq)))
write.table(result,file="C:/Users/luwan/Desktop/UGA/KY_lab/PGS/PGS_regression.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
PGS_coef="NA"; PGS_se="NA";PGS_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";PGS_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"








########_____-----------------
##save diff biomarker for SNP association and APOE4 status analysis
#pheno
pheno_data_diff <- data.frame(0,W5_W8_diff$IID,W5_W8_diff$pTau181_diff,W5_W8_diff$Abeta40_diff,W5_W8_diff$Abeta42_diff,W5_W8_diff$GFAP_diff,W5_W8_diff$NFlight_diff)
colnames(pheno_data_diff) <- c("FID","IID", "ptau_diff","ab40_diff","ab42_diff","gfap_diff","nfl_diff")
# #covar
# PC1_20<-read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/PCA/W5W8_PCA.eigenvec",header=FALSE,sep=" ")
# covar_data <- merge(W5_W8_diff[c("IID", "Age", "gsex")],PC1_20, by.x = "IID", by.y = "V2")
# covar_data <- covar_data[, c("V1", setdiff(names(covar_data), "V1"))]
# colnames(covar_data) <- c("FID","IID", "age","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")

write.table(pheno_data_diff,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/biomarker_diff_for_SNPassociation.txt",col.names = TRUE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
# write.table(covar_data,file="C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/covar_PC_for_SNPassociation.txt",col.names = TRUE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')







#################################
## Add covariates
more_cov <- read.csv("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/more_variables/more_variables/FACHS_W8_covariates_03-12-2025.csv", header = TRUE)
W5_W8_diff_morecov <- merge(W5_W8_diff, more_cov, by.x = "PID", by.y = "PID")

W5_W8_diff_morecov$G8A5023_factor <- ifelse(W5_W8_diff_morecov$G8A5023 <= 12, 0, 1)
W5_W8_diff_morecov$G8A5023_factor <- as.factor(W5_W8_diff_morecov$G8A5023_factor)



#Regression one by one
# fitModel <- glm(pTau181_diff~pop2+gsex+AgePCSCw8+W5_W8_diff_morecov[,i],data=W5_W8_diff_morecov)#,family = "binomial")
# ml_summary <- summary(fitModel)
# diff_ptau_pvalue <- ml_summary$coefficients[2,4]
for (i in 25:31) {
  W5_W8_diff_morecov[, i] <- as.factor(W5_W8_diff_morecov[, i])
}
for (i in 47:56) {
  W5_W8_diff_morecov[, i] <- as.factor(W5_W8_diff_morecov[, i])
}

W5_W8_diff_morecov$Systolic <- (W5_W8_diff_morecov$LM1Systolictop + W5_W8_diff_morecov$LM2Systolictop)/2
W5_W8_diff_morecov$Diastolic <- (W5_W8_diff_morecov$LM1Diastolicbottom + W5_W8_diff_morecov$LM2Diastolicbottom)/2

p_values_df <- data.frame(Variable = character(), p_value = numeric(), stringsAsFactors = FALSE)
for (i in 25:60) {
  sample_size <- sum(!is.na(W5_W8_diff_morecov[, i]))
  fitModel <- glm(pTau181_diff ~ pop2 + gsex + AgePCSCw8 + W5_W8_diff_morecov[, i], 
                  data = W5_W8_diff_morecov)
  ml_summary <- summary(fitModel)
  diff_pop2_pvalue <- ml_summary$coefficients[2, 4]
  diff_pop2_beta <- ml_summary$coefficients[2, 1]
  
  p_values_df <- rbind(p_values_df, data.frame(
    Variable = colnames(W5_W8_diff_morecov)[i], 
    SampleSize = sample_size,
    EffectSize = diff_pop2_beta,
    p_value = diff_pop2_pvalue
  ))
}

# Print the final dataframe
print(p_values_df)




##add Social eco status and income
income_data <- read.csv("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/more_variables/more_variables/Income_Ws5&8_03-18-2025.csv", header = TRUE)
ses_data <- read.csv("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/more_variables/more_variables/SES_Ws5&8_03-18-2025.csv", header = TRUE)
w5_cov <- read.csv("C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/more_variables/more_variables/FACHS_W5_covariates_03-12-2025.csv", header = TRUE)
W5_W8_diff_income_ses <- merge(W5_W8_diff, ses_data, by.x = "PID", by.y = "PID")
W5_W8_diff_income_ses <- merge(W5_W8_diff_income_ses, w5_cov[c("PID", "BMI5")], by.x = "PID", by.y = "PID", all.x = TRUE)
W5_W8_diff_income_ses <- merge(W5_W8_diff_income_ses, income_data, by.x = "PID", by.y = "PID", all.x = TRUE)
for (i in 24:25) {
  W5_W8_diff_income_ses[, i] <- as.numeric(W5_W8_diff_income_ses[, i])
}
for (i in 27:28) {
  W5_W8_diff_income_ses[, i] <- as.numeric(W5_W8_diff_income_ses[, i])
}

p_values_more <- data.frame(Variable = character(), p_value = numeric(), stringsAsFactors = FALSE)
for (i in 24:28) {
  sample_size <- sum(!is.na(W5_W8_diff_income_ses[, i]))
  fitModel <- glm(pTau181_diff ~ pop2 + gsex + AgePCSCw8 + W5_W8_diff_income_ses[, i], 
                  data = W5_W8_diff_income_ses)
  ml_summary <- summary(fitModel)
  diff_pop2_pvalue <- ml_summary$coefficients[2, 4]
  diff_pop2_beta <- ml_summary$coefficients[2, 1]
  
  p_values_more <- rbind(p_values_more, data.frame(
    Variable = colnames(W5_W8_diff_income_ses)[i], 
    SampleSize = sample_size,
    EffectSize = diff_pop2_beta,
    p_value = diff_pop2_pvalue
  ))
}

# Print the final dataframe
print(p_values_more)
