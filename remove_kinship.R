library(bigsnpr)
library(ggplot2)

bedfile <- "C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/GT_for_PCA/W5W8_indiv.bed"

plink2 <- download_plink2("data")

# Relatedness
rel <- snp_plinkKINGQC(
  plink2.path = plink2,
  bedfile.in = bedfile,
  thr.king = 2^-4.5,
  make.bed = FALSE,
  ncores = nb_cores()
)
str(rel)


#remove least related indiv
least_rm_samID_relative <- function(X,ID1=NA,ID2=NA,seed=NA){
  
  pair.now <- X
  rmID_multi.now <- c()
  print(paste0("Total number of pairs: ",nrow(pair.now)))
  
  set.seed(seed)
  
  for(i in 1:nrow(X)){
    
    if(i%%5000==0){print(paste0("Loop number: ",i))}
    
    ## Count the number of repetitions
    total_ID.now <- c(as.character(pair.now[,ID1]),as.character(pair.now[,ID2]))
    table_ID.now <- table(total_ID.now)
    
    #identify which ID(s) has/have the most repetitions.
    max_count.now <- max(table_ID.now)
    if(max_count.now==1){break}
    if(max_count.now>1){
      max_ID.now <- names(table_ID.now)[which(table_ID.now==max_count.now)]  ##sample ID which have most close relatives
      max1_ID.now <- max_ID.now[sample(x=c(1:length(max_ID.now)),size = 1,replace = T)]
      
      rmID_multi.now <- c(rmID_multi.now,max1_ID.now)
      pair.now <- pair.now[-which(pair.now[,ID1]%in%rmID_multi.now | pair.now[,ID2]%in%rmID_multi.now),]
    }
  }
  rmID_multi <- rmID_multi.now
  pair1to1 <- X[-which(X[,ID1]%in%rmID_multi | X[,ID2]%in%rmID_multi),]
  n_pair1to1 <- nrow(pair1to1)
  
  rd.col <- sample(x=c(1,2),size = n_pair1to1,replace = T)
  rmID_rd1in2 <- pair1to1[as.matrix(data.frame(c(1:n_pair1to1),rd.col))]
  
  rmID_related.least <- c(rmID_multi,rmID_rd1in2)
  
  print(paste0("Total number of removed samples: ",length(rmID_related.least)))
  return(rmID_related.least)
}

data <- data.frame(
  ID1 = rel$IID1,
  ID2 = rel$IID2
)
rmID <- least_rm_samID_relative(data,ID1 = "ID1",ID2 = "ID2",seed = 12345)

rmID_list <- data.frame(
  ID1 = "0",
  ID2 = rmID
)
write.table(rmID_list,file = "C:/Users/luwan/Desktop/UGA/KY_lab/Pheno_W5W8_combined/W5W8_min_related_list_to_remove.txt",row.names = F,col.names = F,quote = F,sep = "\t")



# kin <- read.table("C:/Users/luwan/Desktop/UGA/Rotation3/2ndtry_stuff/kinship_3rd_beach_gibbons.txt",header=T,sep="\t")
# kin_filtered <- kin[kin$POP1 != "Gibbons", ]
# kin_filtered <- kin_filtered[kin_filtered$POP2 != "Gibbons", ]
# kin_filtered <- kin_filtered[kin_filtered$POP1 != "Beach", ]
# kin_filtered <- kin_filtered[kin_filtered$POP2 != "Beach", ]