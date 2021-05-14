
library(reshape2)
library(ggplot2)
library(ggpubr)


## Load data--------------------------------------------------------------------------------------------------
genes1<- read.delim(file="D:/Gene_Data/New_data/Nkx2-1_Sftpa1.txt",
                    header = TRUE, sep = "\t")

genes2<- read.delim(file="D:/Gene_Data/New_data/Nkx2-1_Sftpb.txt",
                    header = TRUE, sep = "\t")

genes3<- read.delim(file="D:/Gene_Data/New_data/Nkx2-1_Sftpc.txt",
                    header = TRUE, sep = "\t")

genes4<- read.delim(file="D:/Gene_Data/gene_data_old/copula_proximal_ASE_vgn028_1012h_atac_k27.tab",
                    header = TRUE, sep = "\t")



genes5<- read.delim(file="D:/Gene_Data/gene_data_old/copula_proximal_ASE_vgn028_1012h_atac_rna.tab",
                    header = TRUE, sep = "\t")

genes6<- read.delim(file="D:/Gene_Data/gene_data_old/copula_proximal_ASE_vgn028_1012h_k4_rna.tab",
                    header = TRUE, sep = "\t")


genes4 <- data.frame(genes4[,2],genes4[,4])
colnames(genes4) <- c("ATAC","k27")
genes5 <- data.frame(genes5[,2],genes5[,4])
colnames(genes5) <- c("ATAC","RNA")
genes6 <- data.frame(genes6[,2],genes6[,4])
colnames(genes6) <- c("k4","RNA")




## Histograms--------------------------------------------------------------------------------------------



data1<- melt(genes1)
data2<- melt(genes2)
data3<- melt(genes3)
data4<- melt(genes4)
data5<- melt(genes5)
data6<- melt(genes6)


hist1 <- ggplot(data1,aes(x=value, fill=variable)) + geom_histogram(position = "dodge",alpha=0.55,bins=40)+ scale_fill_brewer(palette = "Set1")+ theme_light()
hist2 <- ggplot(data2,aes(x=value, fill=variable)) + geom_histogram(position = "dodge",alpha=0.55,bins=40)+ scale_fill_brewer(palette = "Set1")+ theme_light()
hist3 <- ggplot(data3,aes(x=value, fill=variable)) + geom_histogram(position = "dodge",alpha=0.55,bins=40)+ scale_fill_brewer(palette = "Set1")+ theme_light()
hist4 <- ggplot(data4,aes(x=value, fill=variable)) + geom_histogram(position = "dodge",alpha=0.55,bins=40)+ scale_fill_brewer(palette = "Set1")+ theme_light()
hist5 <- ggplot(data5,aes(x=value, fill=variable)) + geom_histogram(position = "dodge",alpha=0.55,bins=40)+ scale_fill_brewer(palette = "Set1")+ theme_light()
hist6 <- ggplot(data6,aes(x=value, fill=variable)) + geom_histogram(position = "dodge",alpha=0.55,bins=40)+ scale_fill_brewer(palette = "Set1")+ theme_light()

ggarrange(hist1, hist2,hist3,ncol = 1, nrow = 3)
ggsave("histsA.png", width = 25, height = 25, units = "cm")
ggarrange(hist4,hist5,hist6,ncol = 3, nrow = 2)
ggsave("histsB.png", width = 25, height = 25, units = "cm")







#### Boxplots Single-Cell---------------------------------------------------------------------------------------

old <- data.frame(genes1$Nkx2.1, genes1$Sftpa1, genes2$Sftpb, genes3$Sftpc)
colnames(old) <- c("Nkx2.1","Sftpa1","Sftpb","Sftpc")
dat <- melt(old)
colnames(dat) <- c("genes","value")




dat_mean <- dat %>% 
  group_by(genes) %>% 
  summarize(average = mean(value)) %>%
  ungroup()



dat %>% 
  ggplot(aes(x=genes, y=value, color = genes)) +
  geom_boxplot() +geom_jitter(width=0.15, alpha=0.3)+
  theme_light() + geom_point(data = dat_mean, 
                             mapping = aes(x = genes, y = average),
                             color="black") +
  geom_line(data = dat_mean, 
            mapping = aes(x = genes, y = average, group=1),color = "black")

ggsave("boxplotA.png", width = 25, height = 25, units = "cm")




#### Boxplots Epigenome---------------------------------------------------------------------------------------

new <- data.frame(genes4, genes5, genes6)
colnames(old) <- c("Nkx2.1","Sftpa1","Sftpb","Sftpc")
dat4 <- melt(genes4)
dat5 <- melt(genes5)
dat6 <- melt(genes6)
colnames(dat4) <- c("genes","value")
colnames(dat5) <- c("genes","value")
colnames(dat6) <- c("genes","value")




dat_mean4 <- dat4 %>% 
  group_by(genes) %>% 
  summarize(average = mean(value)) %>%
  ungroup()
dat_mean5 <- dat5 %>% 
  group_by(genes) %>% 
  summarize(average = mean(value)) %>%
  ungroup()
dat_mean6 <- dat6 %>% 
  group_by(genes) %>% 
  summarize(average = mean(value)) %>%
  ungroup()



dat4 %>% 
  ggplot(aes(x=genes, y=value, color = genes)) +
  geom_boxplot() +geom_jitter(width=0.15, alpha=0.3)+
  theme_light() + geom_point(data = dat_mean4, 
                             mapping = aes(x = genes, y = average),
                             color="black") +
  geom_line(data = dat_mean4, 
            mapping = aes(x = genes, y = average, group=1),color = "black")

ggsave("boxplotB.png", width = 25, height = 25, units = "cm")


dat5 %>% 
  ggplot(aes(x=genes, y=value, color = genes)) +
  geom_boxplot() +geom_jitter(width=0.15, alpha=0.3)+
  theme_light() + geom_point(data = dat_mean5, 
                             mapping = aes(x = genes, y = average),
                             color="black") +
  geom_line(data = dat_mean5, 
            mapping = aes(x = genes, y = average, group=1),color = "black")

ggsave("boxplotC.png", width = 25, height = 25, units = "cm")

dat6 %>% 
  ggplot(aes(x=genes, y=value, color = genes)) +
  geom_boxplot() +geom_jitter(width=0.15, alpha=0.3)+
  theme_light() + geom_point(data = dat_mean6, 
                             mapping = aes(x = genes, y = average),
                             color="black") +
  geom_line(data = dat_mean6, 
            mapping = aes(x = genes, y = average, group=1),color = "black")

ggsave("boxplotD.png", width = 25, height = 25, units = "cm")


