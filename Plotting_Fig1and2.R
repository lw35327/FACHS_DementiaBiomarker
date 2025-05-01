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

PCA_1kGP_pop <- merge(PCA_1kGP_ref, popinfo[c("SampleID", "Population","Superpopulation")], by.x = "SampleID", by.y = "SampleID")


PCA_1kGP_Gibbons <- rbind(PCA_1kGP_pop, PCA_gibbons_proj)

library(RColorBrewer)
mypalette <- brewer.pal(5,"Set2")
plinkpca_plot <- function(f_pca){
  require(ggplot2)
  
  df_plot <- f_pca[,c("PC1","PC2","Superpopulation")]
  df_plot$PC1 <- as.numeric(as.character(df_plot$PC1))
  df_plot$PC2 <- as.numeric(as.character(df_plot$PC2))
  # df_plot$Superpopulation <- as.factor(as.character(df_plot$Superpopulation))
  df_plot$Superpopulation <- factor(as.character(df_plot$Superpopulation),
                                    levels = c("AFR","AMR","EAS","EUR","SAS","Cross-sectional"))
  
  p_pca <- ggplot(df_plot,aes(x=PC1,y=PC2,color=Superpopulation)) + geom_point() +
    # geom_point(size = 0.5) +
    scale_x_continuous(name = "PC1") + 
    scale_y_continuous(name = "PC2") + 
    scale_color_manual(values=c(mypalette,"#626262"), name="Groups") +
    
    theme(panel.grid.major =element_line(colour = '#B5B5B5', linetype = 'dashed'),
          panel.grid.minor = element_line(colour = '#B5B5B5', linetype = 'dashed'),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none"
          # legend.title = element_text(size = 14),
          # legend.text = element_text(size = 12)
    )
  
  return(p_pca)
}

# PCA12 = read.table(filename, header=T, as.is=T)
# png(file="population_stratify_pc12.png", width=500, height=500)
# plinkpca_plot(PCA_1kGP_pop)
# dev.off()


plinkpca_plot(PCA_1kGP_Gibbons)













#PC34
plinkpca34_plot <- function(f_pca){
  require(ggplot2)
  
  df_plot <- f_pca[,c("PC3","PC4","Superpopulation")]
  df_plot$PC3 <- as.numeric(as.character(df_plot$PC3))
  df_plot$PC4 <- as.numeric(as.character(df_plot$PC4))
  # df_plot$Superpopulation <- as.factor(as.character(df_plot$Superpopulation))
  df_plot$Superpopulation <- factor(as.character(df_plot$Superpopulation),
                                    levels = c("AFR","AMR","EAS","EUR","SAS","Cross-sectional"))
  
  p_pca <- ggplot(df_plot,aes(x=PC3,y=PC4,color=Superpopulation)) + geom_point() +
    # geom_point(size = 0.5) +
    scale_x_continuous(name = "PC3") + 
    scale_y_continuous(name = "PC4") + 
    scale_color_manual(values=c(mypalette,"#626262"), name="Groups") +
    
    theme(panel.grid.major =element_line(colour = '#B5B5B5', linetype = 'dashed'),
          panel.grid.minor = element_line(colour = '#B5B5B5', linetype = 'dashed'),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)
    )
  
  return(p_pca)
}

plinkpca34_plot(PCA_1kGP_Gibbons)


png(file="C:/Users/luwan/Desktop/UGA/KY_lab/PCA_1kGP_Gibbons/PC12_1kGP_proj_gibbons.png",width=500, height=500)
plinkpca_plot(PCA_1kGP_Gibbons)
dev.off()

png(file="C:/Users/luwan/Desktop/UGA/KY_lab/PCA_1kGP_Gibbons/PC34_1kGP_proj_gibbons.png", width=500, height=500)
plinkpca34_plot(PCA_1kGP_Gibbons)
dev.off()



# library(gridExtra)
p1 <- plinkpca_plot(PCA_1kGP_Gibbons)+ 
  theme(
    # legend.title = element_text(size = 14),  # Change legend title size
    # legend.text = element_text(size = 12),    # Change legend text size
    axis.title = element_text(size = 15),    # Change axis title size
    axis.text = element_text(size = 12)      # Change axis text size
  )
p2 <- plinkpca34_plot(PCA_1kGP_Gibbons)+ 
  theme(
    legend.title = element_text(size = 14),  # Change legend title size
    legend.text = element_text(size = 12),    # Change legend text size
    axis.title = element_text(size = 15),    # Change axis title size
    axis.text = element_text(size = 12)      # Change axis text size
  )
# # Combine the plots side by side
# combined_plot <- grid.arrange(p1, p2, ncol = 2)

# combined_plot <- p1 + p2 +
#   plot_layout(guides = 'collect') + # Optional: collect all legends into one
#   plot_annotation(title = "PCA Projection to 1kGP", tag_levels = 'A')
library(ggpubr)
p12 <- ggarrange(p1,p2,ncol = 2,common.legend = T,legend = "right")

# annotate_figure(p12, top = text_grob("PCA Projection", 
#                                       face = "bold", size = 14))

# Save the combined plot to a PDF file
ggsave("C:/Users/luwan/Desktop/UGA/KY_lab/PCA_1kGP_Gibbons/PC_1kGP_proj_gibbons.pdf", plot = combined_plot, width = 10, height = 5)




#Admixture
##Try project k=4
##Read in files
popinfo <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Admixture/Individual_subpop_list.txt",header=T,sep="\t")
gibbonsID <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Admixture/Gibbons_flipped_admix_allchr.fam",header=F,sep="\t")

admixresult <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Admixture/1kGP_high_coverage_Illumina.indpdt_overlap4_qc_hw_ld_allchr.4.Q",header=F,sep=" ")
gibbonsresult <- read.table("C:/Users/luwan/Desktop/UGA/KY_lab/Admixture/Gibbons_flipped_admix_allchr.4.Q",header=F,sep=" ")

##Try concat Q file and population info file
admixresult$Individual <- popinfo$Individual
gibbonsresult$Individual <- gibbonsID$V2
all.equal(admixresult$Individual,popinfo$Individual)
all.equal(gibbonsresult$Individual,gibbonsID$V2)
# Keep only Combined dataset
gibbonsresult <- gibbonsresult[gibbonsresult$Individual %in% unrel_list$V2, ]


refadmix_merge <- merge(admixresult,popinfo,by.y = "Individual") 
gibbonsresult$Population <- "Cross-sectional"
gibbonsresult$SubPopulation <- "Cross-sectional"

library(dplyr)

tbl <- rbind(refadmix_merge,gibbonsresult)

str(tbl)


tbl_order <- tbl[rev(order(tbl$V2,tbl$V3,tbl$V4,tbl$V1)),] #change this

t.tbl_order.anc <- tbl_order[which(tbl_order$Population=="AFR"),]
ID.AFR_ord <- t.tbl_order.anc$Individual[order(t.tbl_order.anc$V2,t.tbl_order.anc$V1,t.tbl_order.anc$V4,t.tbl_order.anc$V3)] #,decreasing = F)]
t.tbl_order.anc <- tbl_order[which(tbl_order$Population=="AMR"),]
ID.AMR_ord <- t.tbl_order.anc$Individual[order(t.tbl_order.anc$V4,t.tbl_order.anc$V3,t.tbl_order.anc$V1,t.tbl_order.anc$V2)]
t.tbl_order.anc <- tbl_order[which(tbl_order$Population=="SAS"),]
ID.SAS_ord <- t.tbl_order.anc$Individual[order(t.tbl_order.anc$V3,t.tbl_order.anc$V4,t.tbl_order.anc$V1,t.tbl_order.anc$V2)]
t.tbl_order.anc <- tbl_order[which(tbl_order$Population=="EAS"),]
ID.EAS_ord <- t.tbl_order.anc$Individual[order(t.tbl_order.anc$V2,t.tbl_order.anc$V1,t.tbl_order.anc$V4,t.tbl_order.anc$V3)]
t.tbl_order.anc <- tbl_order[which(tbl_order$Population=="EUR"),]
ID.EUR_ord <- t.tbl_order.anc$Individual[order(t.tbl_order.anc$V2,t.tbl_order.anc$V1,t.tbl_order.anc$V4,t.tbl_order.anc$V3)]
t.tbl_order.anc <- tbl_order[which(tbl_order$Population=="Cross-sectional"),]
ID.COMB_ord <- t.tbl_order.anc$Individual[order(t.tbl_order.anc$V2,t.tbl_order.anc$V1,t.tbl_order.anc$V4,t.tbl_order.anc$V3)]


# tbl_order <- tbl_order[order(tbl_order$Population),]




library(reshape2)
cols <- colnames(tbl_order)
tbl_order <- tbl_order[,c("Individual","V1","V2","V3","V4")] 
tbl_order_melt <- melt(data=tbl_order[,1:5],id.vars="Individual") 
colnames(tbl_order_melt) <- c("ID","variable","value")

tbl_order_melt$Population <- tbl$Population[match(tbl_order_melt$ID,tbl$Individual)]
tbl_order_melt$SubPopulation <- tbl$SubPopulation[match(tbl_order_melt$ID,tbl$Individual)]

tbl_order_melt$ID <- factor(tbl_order_melt$ID,ordered = F,levels = c(ID.AFR_ord,ID.AMR_ord,ID.SAS_ord,ID.EAS_ord,ID.EUR_ord,ID.COMB_ord))
tbl_order_melt$variable <- factor(as.character(tbl_order_melt$variable),ordered = F,levels = rev(c("V2","V1","V3","V4")))  #order to color?

# assign color
#plate <- brewer.pal(10,"Paired")
plate <- brewer.pal(9,"Set1")

# id_levels <- tbl_order_melt$ID[1:(length(tbl_order_melt$ID)/4)] #change this

tbl_order_melt$Population <- factor(tbl_order_melt$Population,ordered = F,levels=c("AFR","AMR","SAS","EAS","EUR","Cross-sectional")) #order of 5 plots
str(tbl_order_melt)

#pdf("afr_UKB.k2.supervised.structure.v1.pdf",width = 8,height = 6)
gibbons_k4_projection <- ggplot()+
  geom_bar(data=tbl_order_melt, 
           aes(x=ID,
               y=value,
               fill=variable,
               color=NULL),
           #group=ID),
           stat="identity",
           width = 1)+
  scale_y_continuous(limits = c(0,1.001),breaks = seq(0,1.001,0.25))+
  # ylim(0,1.001)+
  theme_bw()+
  theme(panel.background = element_blank(),panel.border=element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        axis.text.x =element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 14),
        legend.background=element_rect(fill = 'white', colour = 'black'),
        legend.title=element_text(size=14),
        legend.text=element_text(size=15),
        legend.position= "none",#"right",#c(0.11,0.14),
        legend.box = "horizontal"
  )+
  xlab("")+
  ylab("Ancestry Proportion")+
  ggtitle("k=4")+
  labs(fill="Ancestry")+#,colour="Ancestry")+
  scale_fill_manual(values =c(plate[2],plate[4],plate[5],plate[3]),
                    breaks=c("V1","V2","V3","V4"))+
  facet_wrap(~Population,
             scales ="free_x",
             nrow = 1,strip.position = "bottom"
  )

print(gibbons_k4_projection)

library(patchwork)

library(gridExtra)
grid.arrange(ggarrange(p1,labels = c("A")),
             ggarrange(p2,labels = c("B")),
             ggarrange(gibbons_k4_projection,labels = c("C")),
             layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),c(3,3,3,3,3,3,3,3,3,3,3,3,3)))


#save to pdf
pdf("C:/Users/luwan/Desktop/UGA/KY_lab/manu/alz_submission1.6/FIGURE1.pdf", width = 10, height = 8)  # Adjust width and height as needed

# Create the plot
grid.arrange(
  ggarrange(p1, labels = c("A")),
  ggarrange(p2, labels = c("B")),
  ggarrange(gibbons_k4_projection, labels = c("C")),
  layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2),
    c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
  )
)

# Close the PDF device
dev.off()
