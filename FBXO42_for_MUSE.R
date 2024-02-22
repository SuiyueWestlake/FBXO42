#MUSE
###
Cao=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/fbxo42.csv")
Cao=na.omit(Cao)

NOTCH_analysis=read.csv("/Users/tangsuiyue/Documents/Documents/Project/MUSE/data/input_datasets/NOTCH_pathway_analisis.csv")
NOTCH_Controls=NOTCH_analysis[21032:72348,1:5]
colnames(Cao)=colnames(NOTCH_Controls)
Cao_data_for_MUSE_2=rbind(Cao,NOTCH_Controls)
#write.csv(Cao_data_for_MUSE_2,"/Users/tangsuiyue/Documents/Documents/Project/Cao_QI/Cao_data_for_MUSE_2.csv",row.names = F)
#write.csv(Cao_data_for_MUSE_2,"/Users/tangsuiyue/Documents/Documents/Project/MUSE/data/input_datasets/Cao_data_for_MUSE_2.csv",row.names = F)


###
tranform_log <- NULL
Cao_modify_prey <- as.character(Cao_data_for_MUSE_2[,4])
preys_vect <- c(paste0("KRT", 1:9), paste0("HSPA", 1:9),
                paste0("HIST", 1:4, "H"), paste0("TUBA", 1:9),
                paste0("TUBB", 1:9), "HNRNPA", "HSPB", "EEF1", "HSP90A",
                "ACTA", "KPNA", "HLA-A", "HLA-B", "MRPS", "MRPL")

## combine the homology proteins: KRT*
for(i in 1:length(preys_vect)){
  prey_name <- preys_vect[i]
  index2 <- grep(pattern = prey_name , Cao[,4])
  print(c(i, length(index2)))
  if(length(index2) > 0){
    new_prey_name <- as.character(Cao_data_for_MUSE_2[index2[1] ,4])
    Cao_modify_prey[index2] <- new_prey_name
    tranform_log <- rbind(tranform_log,
                          c(paste0(prey_name, "*"), new_prey_name))
    print(c(paste0(prey_name, "*"), new_prey_name))
  }
}
Cao_modify <- Cao_data_for_MUSE_2
Cao_modify[, 4] <- Cao_modify_prey
#write.csv(tranform_log,
#          file=paste0(pathstr_save, "notch_combprey.log"),
#          row.names=FALSE)
### delete krt
head(Cao_modify)
index3 <- grep(pattern = "KRT", Cao_modify[, 4])
length(index3)
Cao_modified <- Cao_modify[- index3, ]
write.csv(Cao_modified,
          file="/Users/tangsuiyue/Documents/Documents/Project/FBXO42/fbxo42_for_MUSE.csv",
          row.names=FALSE)


