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

##Alpha Diversity
#Shannon plot##
library(qiime2R)
metadata1<-read_q2metadata("metadata.tsv")
view(metadata1)

shannon<-read_qza("shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged
view(shannon)
gplots::venn(list(metadata=metadata1$SampleID, shannon=shannon$SampleID)) #from diagram, we see that all samples in metadata have assigned Shannon diversity value

metadata2<- metadata1 %>% left_join(shannon) 
head(metadata2)
write.table(metadata2, file = "shannonmetadata.csv", sep = ",")

#For Shannon diversity,  pairwise test with Wilcoxon rank-sum test, corrected by FDR method:
rich = estimate_richness(ps, measures = c("Observed", "Shannon"))
wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(ps)$diet, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
tab.shannon

##Beta Diversity



# PCOA plots using the same phyloseq object
#PCOA: for each diet based on experimental days**
###CF**###

ps.sub.CTRL <- subset_samples(ps, diet %in% c("D0", "CF"))

dist = phyloseq::distance(ps.sub.CTRL, method="bray")
ordination = ordinate(ps.sub.CTRL, method="PCoA", distance=dist)

p1<-plot_ordination(ps.sub.CTRL, ordination, color="age", shape = "line") + 
  theme_bw()+ geom_point(size=3)+ ggtitle("CF") + stat_ellipse(aes(group = age), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ labs(shape="Line", colour="Experimental Day")+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank())

print(P1)

###NUS**###
  
ps.sub.NUS <- subset_samples(ps, diet %in% c("D0", "NUS"))
#str(ps.sub.NUS)
#sample_data(ps.sub.NUS)$age_factor <- as.factor(sample_data(ps.sub.NUS)$age)
#sample_data(ps.sub.NUS)$line_factor <- as.factor(sample_data(ps.sub.NUS)$line)
dist = phyloseq::distance(ps.sub.NUS, method="bray")
ordination = ordinate(ps.sub.NUS, method="PCoA", distance=dist)

p2<- plot_ordination(ps.sub.NUS, ordination, color="age", shape = "line") + geom_point(size=3) +
  theme_bw() + ggtitle("NUS") + stat_ellipse(aes(group = age_factor), linetype = 2) +
  scale_color_manual(values=c("orange", "deeppink", "brown4"),limits = c("D0", "D5", "D10"))+ scale_shape_manual(values=c(16, 1))+
  theme(strip.background = element_blank()) + labs(shape="Line", colour="Experimental Day")

print(p2)

############
