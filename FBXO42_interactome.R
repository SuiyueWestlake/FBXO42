###FBXO42 interactome
MUSE_Result=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/Interactome/FBXO42_MUSE_Score.csv")
MUSE_Result$Total_PSM=apply(MUSE_Result[,7:8],1,sum)
CRAPome_cut_off=411*0.5
Contaminants_list=MUSE_Result[which(MUSE_Result[,9]>=CRAPome_cut_off),2]
Delete_index=which(MUSE_Result$Prey %in% Contaminants_list)
Sub_network=MUSE_Result[-Delete_index,]
Sub_network=Sub_network[order(Sub_network$Total_PSM,decreasing = T),]
#Sub_network$Total_PSM=scale(Sub_network$Total_PSM)
#Sub_network$weighted_score=Sub_network$Score*Sub_network$Total_PSM
Sub_network=Sub_network[-3,]
Sub_network=Sub_network[1:100,]
write.csv(Sub_network,"/Users/tangsuiyue/Documents/Documents/Project/FBXO42/Interactome/FBXO42_interactome.csv",quote = F,row.names = F)

