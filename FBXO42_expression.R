###FBXO42 expression box plot
library(ggplot2)
FBXO42_expression=read.csv("/Users/tangsuiyue/Documents/Documents/Project/FBXO42/FBXO42_mRNA.csv")
colnames(FBXO42_expression)=c("Patients","Cancer_type","Group","Expression")
unique(FBXO42_expression$Cancer_type)
FBXO42_expression_HNSC=FBXO42_expression[which(FBXO42_expression$Cancer_type == "HNSC"),]
FBXO42_expression_BRCA=FBXO42_expression[which(FBXO42_expression$Cancer_type == "BRCA"),]
FBXO42_expression_LAML=FBXO42_expression[which(FBXO42_expression$Cancer_type == "LAML"),]
FBXO42_expression_DLBC=FBXO42_expression[which(FBXO42_expression$Cancer_type == "DLBC"),]

ggplot(data = FBXO42_expression) + 
  geom_boxplot(aes(x = Cancer_type, y = Expression,
                   fill = factor(Group))) +
  theme_classic()

ggplot(data = FBXO42_expression_HNSC) + 
  geom_boxplot(aes(x = Group, y = Expression)) +
  scale_fill_brewer(palette = "Set3")
ggplot(data = FBXO42_expression_BRCA) + 
  geom_boxplot(aes(x = Group, y = Expression))
ggplot(data = FBXO42_expression_LAML) + 
  geom_boxplot(aes(x = Group, y = Expression))
ggplot(data = FBXO42_expression_DLBC) + 
  geom_boxplot(aes(x = Group, y = Expression))



p2 <- ggplot(FBXO42_expression, aes
             (x=V2,Y=V3)) +
  geom_hline(yintercept = seq(25, 100, 25), color = 'gray') +
  geom_box(color = "black", width = .7, position = 'stack') +
  labs( y = 'Relative abundance (%)') +
  scale_fill_brewer(palette = "Set3")+
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2
