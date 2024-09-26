p1 <- c("tidyverse", "vegan", "BiocManager")
p2 <- c("phyloseq", "ANCOMBC", "DESeq2", "ComplexHeatmap")
load_package <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    ifelse(p %in% p1, 
           install.packages(p, repos = "http://cran.us.r-project.org/"), 
           BiocManager::install(p))
  }
  library(p, character.only = TRUE, quietly = TRUE)
}
invisible(lapply(c(p1,p2), load_package))
install.packages("ggpubr", repos = "https://cloud.r-project.org/", dependencies = TRUE)
warnings()

#Building a phyloseq object

otu <- read.table(file = "feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")
taxonomy <- read.table(file = "taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

# clean the taxonomy, Greengenes format
tax <- taxonomy %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)
#The stringsAsFactors = FALSE argument means that all character columns in the 
#tax.clean data frame will remain as character data types instead of being converted to factors. 

tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

metadata <- read.table(file = "metadata.tsv", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)
TREE = read_tree("tree.nwk")
# merge the data
ps <- phyloseq(OTU, TAX, SAMPLE,TREE)

######Alpha Diversity#########
##Shannon plot##
library(qiime2R)
#metadata1<-read_q2metadata("metadata.tsv")
#view(metadata1)

#shannon<-read_qza("shannon_vector.qza")
#shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged
#view(shannon)
#gplots::venn(list(metadata=metadata1$SampleID, shannon=shannon$SampleID)) #from diagram, we see that all samples in metadata have assigned Shannon diversity value

#metadata2<- metadata1 %>% left_join(shannon) 
#head(metadata2)
#write.table(metadata2, file = "shannonmetadata.csv", sep = ",")

##Shannon boxplots (line?)
plot_richness(ps, x="line", measures=c("Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

##Shannon metadata
library(phyloseq)
library(qiime2R)
shannon <- rownames_to_column(estimate_richness(ps, split = TRUE, measures = "Shannon"), var = "SampleID")
metadata1<-read_q2metadata("metadata.tsv")
gplots::venn(list(metadata=metadata1$SampleID, shannon=shannon$SampleID))
metadata_shannon<- metadata1 %>% left_join(shannon)

head(metadata_shannon)

write.table(metadata_shannon, file = "shannonmetadata.csv", sep = ",")

#Plot for D0 and D5**
  
met_D5 <- metadata_shannon[c(1:33,64:95),c(1:17)]
view(met_D5)

met_D5$diet_f = factor(met_D5$diet, levels=c('D0','CF','OKA','PKM','RIB','SBM', 'M', 'P', 'MHP', 'MLP', 'NUS'))
str(met_D5) #check col types
legend_title <- "Line" #rename the legend

plot1<-met_D5 %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x= line, y=Shannon, fill=line)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.1, height=0) +
  coord_cartesian(ylim=c(0,7)) + # adjust y-axis
  facet_grid(~diet_f) + # create a panel for each body site
  xlab("Diet") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(legend_title, values=c("black","white")) + #use renamed legend and specify custom colors
  theme(legend.position="right") + theme(axis.text.x=element_blank())+ theme(text=element_text(size = 15)) 
plot1
**Plot for D0 and D10**
  
  ```
met3 <- metadata1[c(1:3,34:66,96:125),c(1:17)]
view(met3)

met3$diet_f = factor(met3$diet, levels=c('D0','CF','OKA','PKM','RIB','SBM', 'M', 'P', 'MHP', 'MLP', 'NUS'))
legend_title <- "Line" #rename the legend

plot2<-met3 %>%
  filter(!is.na(shannon_entropy)) %>%
  ggplot(aes(x= line, y=shannon_entropy, fill=line)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.1, height=0) +
  coord_cartesian(ylim=c(0,7)) + # adjust y-axis
  facet_grid(~diet_f) + # create a panel for each body site
  xlab("Diet") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(legend_title, values=c("black","white")) + #use renamed legend and specify custom colors
  theme(legend.position="right") + theme(axis.text.x=element_blank())+ theme(text=element_text(size = 15)) 



plot <- ggarrange(plot1 + rremove("ylab") + rremove("xlab"), plot2 + rremove("ylab") + rremove("xlab"), # remove axis labels from plots
                  ncol = 1, nrow = 2,
                  common.legend = TRUE, legend = "bottom")



annotate_figure(plot, left = textGrob("Shannon Diversity", rot = 90, vjust = 0.5, gp = gpar(cex = 1.3)),
                bottom = textGrob("Diet", gp = gpar(cex = 1.3)))

```

# significance of variables

```
hist(metadata1$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)

aov.shannon.line = aov(shannon_entropy ~ line, data=metadata1)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.shannon.line)

aov.shannon.diet = aov(shannon_entropy ~ diet, data=metadata1)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.shannon.diet)


aov.shannon.age = aov(shannon_entropy ~ age, data=metadata1)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.shannon.age)


##For Shannon diversity,  pairwise test with Wilcoxon rank-sum test, corrected by FDR method:
metadata_clean <- metadata %>%
  mutate(diet = paste0(toupper(substr(metadata$line, 1, 2)), "_", diet)) 
SAMPLE1 = sample_data(metadata_clean)
wilcox_ps <- phyloseq(OTU, TAX, SAMPLE1,TREE)
rich = estimate_richness(ps, measures = c("Observed", "Shannon"))
wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(wilcox_ps)$diet, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

tab.shannon
# Exporting to CSV
write.csv(tab.shannon, "tab_shannon.csv", row.names = FALSE)

##Beta Diversity



# PCOA plots using the same phyloseq object
library(ggplot2)
library(ggpubr)
#PCOA: for each diet based on experimental days**
###CF**###

ps.sub.CTRL <- subset_samples(ps, diet %in% c("D0", "CF"))

dist = phyloseq::distance(ps.sub.CTRL, method="bray")
ordination = ordinate(ps.sub.CTRL, method="PCoA", distance=dist)

p1<-plot_ordination(ps.sub.CTRL, ordination, color="age", shape = "line") + 
  theme_bw()+ geom_point(size=3)+ ggtitle("CF") + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p1)

###NUS**###
ps.sub.NUS <- subset_samples(ps, diet %in% c("D0", "NUS"))
#str(ps.sub.NUS)
#sample_data(ps.sub.NUS)$age_factor <- as.factor(sample_data(ps.sub.NUS)$age)
#sample_data(ps.sub.NUS)$line_factor <- as.factor(sample_data(ps.sub.NUS)$line)
dist = phyloseq::distance(ps.sub.NUS, method="bray")
ordination = ordinate(ps.sub.NUS, method="PCoA", distance=dist)

p2<- plot_ordination(ps.sub.NUS, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("NUS") + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank()) + labs(shape="Line", colour="Experimental Day")

print(p2)


####MHP**###
ps.sub.MHP <- subset_samples(ps, diet %in% c("D0", "MHP"))

dist = phyloseq::distance(ps.sub.MHP, method="bray")
ordination = ordinate(ps.sub.MHP, method="PCoA", distance=dist)

p3<- plot_ordination(ps.sub.MHP, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("MHP")  + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p3)

###MLP**
ps.sub.MLP <- subset_samples(ps, diet %in% c("D0", "MLP"))

dist = phyloseq::distance(ps.sub.MLP, method="bray")
ordination = ordinate(ps.sub.MLP, method="PCoA", distance=dist)

p4 <- plot_ordination(ps.sub.MLP, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("MLP")  + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p4)

###**M**
ps.sub.M <- subset_samples(ps, diet %in% c("D0", "M"))

dist = phyloseq::distance(ps.sub.M, method="bray")
ordination = ordinate(ps.sub.M, method="PCoA", distance=dist)

p5 <- plot_ordination(ps.sub.M, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("M")  + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p5)

###P**
ps.sub.P <- subset_samples(ps, diet %in% c("D0", "P"))

dist = phyloseq::distance(ps.sub.P, method="bray")
ordination = ordinate(ps.sub.P, method="PCoA", distance=dist)

p6 <- plot_ordination(ps.sub.P, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("P")  + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p6)

###OKA**
ps.sub.OKA <- subset_samples(ps, diet %in% c("D0", "OKA"))

dist = phyloseq::distance(ps.sub.OKA, method="bray")
ordination = ordinate(ps.sub.OKA, method="PCoA", distance=dist)

p7 <- plot_ordination(ps.sub.OKA, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("OKA")  + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p7)

###PKM**
ps.sub.PKM <- subset_samples(ps, diet %in% c("D0", "PKM"))

dist = phyloseq::distance(ps.sub.PKM, method="bray")
ordination = ordinate(ps.sub.PKM, method="PCoA", distance=dist)

p8 <- plot_ordination(ps.sub.PKM, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("PKM")  + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p8)

###RIB**
ps.sub.RIB <- subset_samples(ps, diet %in% c("D0", "RIB"))

dist = phyloseq::distance(ps.sub.RIB, method="bray")
ordination = ordinate(ps.sub.RIB, method="PCoA", distance=dist)

p9 <- plot_ordination(ps.sub.RIB, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("RIB")  + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p9)

###SBM**
ps.sub.SBM <- subset_samples(ps, diet %in% c("D0", "SBM"))

dist = phyloseq::distance(ps.sub.SBM, method="bray")
ordination = ordinate(ps.sub.SBM, method="PCoA", distance=dist)

p10 <- plot_ordination(ps.sub.SBM, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("SBM")  + stat_ellipse(aes(group = age), linetype = 2) + 
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(p10)


ggarrange(p1, p2, p8, p10, p5, p3, p9, p6, p7, p4, nrow=5, ncol=2, common.legend = TRUE, legend="bottom")
