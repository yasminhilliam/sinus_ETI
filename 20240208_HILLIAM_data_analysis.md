---
title: "Cystic fibrosis pathogens persist in the upper respiratory tract following initiation of elexacaftor/tezacaftor/ivacaftor therapy"
subtitle: "Analysis coding notebook"
author: "Yasmin Hilliam, PhD"
output: html_notebook
---
```{r setup, include = FALSE}
knitr::opts_chunk$set(echo = TRUE,
                      collapse = TRUE)
```
```{r, include = FALSE}
library(tidyverse)
library(lubridate)
library(phyloseq)
library(vegan)
library(ape)
library(decontam)
library(ggtext)
library(lme4)
library(ggh4x)
library(Maaslin2)
library(cowplot)
library(ggrepel)

theme_set(theme_bw()) # set theme
theme_replace(axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20,angle = 90),
              axis.text.x = element_text(size = 18),
              axis.text.y = element_text(size = 18),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 18),
              strip.text = element_text(size = 16, margin = margin(10,10,10,10))) # set individual theme elements for plotting
```
##### Data tables
```{r, message = FALSE}
read.table(file.path(wd,
                     "primers.tsv"))
read.table(file.path(wd,
                     "sample_data",
                     "adult_metadata.tsv"),
           header = TRUE)
```
##### Import primer-sorted reads as *QIIME2* artifacts into separate folders
```{bash, message = FALSE, results = FALSE}
for i in $(cut -f1 primers.tsv); \
do qiime tools import \
--type EMPPairedEndSequences \
--input-path reads/primer_sorted/"$i" \
--output-path qiime/"$i"/paired-end-sequences.qza; \
done
```
##### Demultiplex primer-sorted sequences
```{bash, message = FALSE, results = FALSE}
for i in $(cut -f1 primers.tsv); \
do qiime demux emp-paired \
--i-seqs qiime/"$i"/paired-end-sequences.qza \
--m-barcodes-file adult_metadata.tsv \
--m-barcodes-column barcode \
--p-no-golay-error-correction \
--o-per-sample-sequences qiime/"$i"/demultiplexed-seqs.qza \
--o-error-correction-details /home/yah71/adult_CF_CRS/final_data/qiime/"$i"/error-correction.qza; \
done
	
for i in $(cut -f1 primers.tsv); \
do qiime demux summarize \
--i-data qiime/"$i"/demultiplexed-seqs.qza \
--o-visualization qiime/"$i"/demux-viz.qzv; \
done
```
##### Denoise 16S rRNA amplicons with *DADA2* plug-in
```{bash, message = FALSE, results = FALSE}
qiime dada2 denoise-paired \
--i-demultiplexed-seqs qiime/V4/demultiplexed-seqs.qza \
--p-trim-left-f 13 \
--p-trim-left-r 13 \
--p-trunc-len-f 250 \
--p-trunc-len-r 250 \
--o-table qiime/V4/denoise-table.qza \
--o-representative-sequences qiime/V4/rep-seqs.qza \
--o-denoising-stats qiime/V4/denoising-stats.qza
	
qiime feature-table summarize \
--i-table qiime/V4/denoise-table.qza \
--o-visualization qiime/V4/feature-table.qzv \
--m-sample-metadata-file adult_metadata.tsv
	
qiime feature-table tabulate-seqs \
--i-data qiime/V4/rep-seqs.qza \
--o-visualization qiime/V4/rep-seqs.qzv
```
##### Identify and remove chimeric sequences using *VSEARCH* plug-in
```{bash, message = FALSE, results = FALSE}
qiime vsearch uchime-denovo \
--i-table qiime/V4/denoise-table.qza  \
--i-sequences qiime/V4/rep-seqs.qza  \
--output-dir qiime/V4/uchime-out

qiime metadata tabulate \
--m-input-file qiime/V4/uchime-out/stats.qza \
--o-visualization qiime/V4/uchime-out/stats.qzv

qiime feature-table filter-features \
--i-table qiime/V4/denoise-table.qza \
--m-metadata-file qiime/V4/uchime-out/nonchimeras.qza \
--o-filtered-table qiime/V4/uchime-out/table-nonchimeric-wo-borderline.qza

qiime feature-table filter-seqs \
--i-data qiime/V4/rep-seqs.qza \
--m-metadata-file qiime/V4/uchime-out/nonchimeras.qza \
--o-filtered-data qiime/V4/uchime-out/rep-seqs-nonchimeric-wo-borderline.qza

qiime feature-table summarize \
--i-table qiime/V4/uchime-out/table-nonchimeric-wo-borderline.qza \
--o-visualization qiime/V4/uchime-out/table-nonchimeric-wo-borderline.qzv
```
##### Generate phylogenetic tree
```{bash}
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences qiime/V4/uchime-out/rep-seqs-nonchimeric-wo-borderline.qza \
--o-alignment qiime/V4/aligned-rep-seqs-no-chi.qza \
--o-masked-alignment qiime/V4/masked-aligned-rep-seqs-no-chi.qza \
--o-tree qiime/V4/unrooted-tree-no-chi.qza \
--o-rooted-tree qiime/V4/rooted-tree-no-chi.qza
```
##### Perform alpha rarefaction to identify optimal sampling depth
```{bash, message = FALSE, results = FALSE}
qiime diversity alpha-rarefaction \
--i-table qiime/V4/uchime-out/table-nonchimeric-wo-borderline.qza  \
--i-phylogeny qiime/V4/rooted-tree-no-chi.qza \
--p-max-depth 909638 \
--p-min-depth 11730 \
--p-steps 20 \
--m-metadata-file adult_metadata.tsv \
--o-visualization qiime/V4/alpha-rarefaction-no-chi.qzv
```
##### Plot sampling depth
```{r, message = FALSE, include = FALSE}
read.csv(file.path(wd,
                   "qiime",
                   "observed_features.csv")) -> oFeatures

oFeatures %>%
  pivot_longer(2:201,
               names_to = "depth",
               values_to = "features") -> oFeatures

str_remove_all(oFeatures$depth,"depth.") -> oFeatures$depth

str_remove_all(oFeatures$depth,"_iter.\\d+") -> oFeatures$depth

oFeatures %>%
  group_by(sample.id,depth) %>%
  mutate(mn_features = mean(features,
                            na.rm = FALSE)) %>%
  ungroup() %>%
  select(sample.id,barcode,patientid,sampleloc,date,depth,mn_features) %>%
  distinct() -> oFeatures

as.numeric(oFeatures$depth) -> oFeatures$depth

ggplot(oFeatures,
       aes(x = depth,
           y = mn_features,
           color = sampleloc,
           group = sample.id)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  scale_color_discrete(name = "Sample location",
                       labels = c("Environmental control",
                                  "H2O control",
                                  "Sinus",
                                  "Sputum",
                                  "Throat")) +
  scale_y_continuous(name = "Mean observed features") +
  scale_x_continuous(name = "Sampling depth")
```
##### Rarefy samples to specified depth with *QIIME2*
```{bashmessage = FALSE, results = FALSE}
qiime feature-table rarefy \
--i-table qiime/V4/uchime-out/table-nonchimeric-wo-borderline.qza \
--p-sampling-depth 11500 \
--o-rarefied-table qiime/V4/rarefied-denoise-table-no-chi.qza
```
##### Download SILVA database and build training set
```{bash, message = FALSE, results = FALSE}
wget -O "/home/yah71/build/silva-138-99-515-806-nb-classifier.qza" "https://data.qiime2.org/2021.4/common/silva-138-99-515-806-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
--i-classifier build/silva-138-99-515-806-nb-classifier.qza \
--i-reads qiime/V4/uchime-out/rep-seqs-nonchimeric-wo-borderline.qza \
--o-classification qiime/V4/taxonomy-no-chi.qza
```
##### Visualize taxonomy assignments with *QIIME2*
```{bash, message = FALSE, results = FALSE}
qiime metadata tabulate \
--m-input-file qiime/V4/taxonomy-no-chi.qza \
--o-visualization qiime/V4/taxonomy-no-chi.qzv

qiime taxa barplot \
--i-table qiime/V4/rarefied-denoise-table-no-chi.qza \
--i-taxonomy qiime/V4/taxonomy-no-chi.qza \
--m-metadata-file adult_metadata.tsv \
--o-visualization qiime/V4/rarefied-taxa-barplots-no-chi.qzv
```
##### Set seed
```{r, message = FALSE}
set.seed(9792)
```
```{r, include = FALSE}
wd <- c("/Users/f006n1k/Documents/adult_CF_CRS/") # set working directory variable
setwd(wd) # set working directory
```
##### Import and process *QIIME2* export data
```{r, message = FALSE}
asvTable <- read.table(file.path(wd,
                                 "16S",
                                 "feature-table-no-chi.tsv"),
                  sep = "\t",
                  row.names = 1,
                  header = TRUE) # import feature table

colnames(asvTable) <- gsub("(\\d{1,2}_\\d{1,2}_)(\\d{2})(\\d{2})", 
                           "\\1\\3", 
                           colnames(asvTable)) # remove 4 digit date format from first row of asv table

colnames(asvTable) <- gsub("GESCF13_sinus_3_8_17", # sample date mis-recorded - no samples in 2017
                           "GESCF13_sinus_3_8_18", # replace incorrect index with correct index
                           colnames(asvTable))

asvTable[] <- sapply(asvTable, as.numeric) # convert asv values to numeric
```
##### Import and process metadata file
```{r, message = FALSE}
metadata <- read.csv(file.path(wd,
                               "sample_data",
                               "2023-07-31_metadata.csv")) # read in metadata file

setdiff(metadata$index,colnames(asvTable)) # 391 observations in asvTable, 392 observations in metadata file - identify missing index
# "EVNCTRL14_11_12_21" - environmental control that did not amplify during sequencing

metadata <- metadata[metadata$index!="EVNCTRL14_11_12_21",] # remove missing sample

# package decontam does not accept 0 values for DNA concentration - 0s were included in dna_conc
# column where samples fell below LOD for Qubit HS dsDNA kit
# quant column sets 0 values to LOD for kit (0.005 ng/uL DNA)

metadata <- metadata %>%
  mutate(quant = if_else(dna_conc == 0, 
                         0.005,
                         dna_conc))

print(metadata[metadata$index == "GESCF13_sinus_3_8_17",])

metadata$index <- gsub("GESCF13_sinus_3_8_17", # sample date mis-recorded - no samples in 2017
                       "GESCF13_sinus_3_8_18", # replace incorrect index with correct index
                       metadata$index)

print(metadata[metadata$date == "2017-03-08",]) # only one sample with this incorrect date so can be replaced
metadata$date <- gsub("2017-03-08",
                      "2018-03-08",
                      metadata$date)
```
##### Import and process ASV taxonomy file
```{r, message = FALSE}
taxonomy <- read.table(file.path(wd,
                                 "16S",
                                 "taxonomy-no-chi.tsv"),
           sep = "\t",
           row.names = 1,
           header = TRUE) # import taxonomy table

taxonomy <- taxonomy %>% # split Taxon column in to multiple columns
  mutate(Taxon = str_replace_all(Taxon,".\\_","")) %>%
  mutate(Taxon = str_replace_all(Taxon, "\\.","")) %>%
  mutate(Taxon = str_replace_all(Taxon, "\\;","")) %>%
  mutate(Taxon = str_replace_all(Taxon, "\\ ","")) %>%
  mutate(Taxon = str_replace_all(Taxon, "\\[","")) %>%
  mutate(Taxon = str_sub(Taxon,2)) %>%
  separate(Taxon,
           into=c("kingdom","phylum","class","order","family","genus"),
           sep="_") %>%
  mutate(phylum=str_replace_all(phylum,"^$","Unassigned")) %>%
  mutate(class=str_replace_all(class,"^$","Unassigned")) %>%
  mutate(order=str_replace_all(order,"^$","Unassigned")) %>%
  mutate(family=str_replace_all(family,"^$","Unassigned")) %>%
  mutate(genus=str_replace_all(genus,"^$","Unassigned"))

taxonomy[is.na(taxonomy)] <- "Unassigned" # replace NA values with Unassigned

taxonomy <- as.matrix(taxonomy) # set as matrix
```
##### Set up data for *phyloseq*
```{r, message = FALSE}
samp <- metadata # copy metadata object to new samp object

rownames(samp) <- metadata$index # set index column as rownames 

otu <- otu_table(asvTable, taxa_are_rows = TRUE) # generate otu table object
samp <- sample_data(samp) # generate sample data object
tax <- tax_table(taxonomy) # generate taxonomy table object

ps <- merge_phyloseq(otu,samp,tax) # merge into phyloseq object

sample_data(ps)
ps
```
##### *phyloseq* analysis
```{r, message = FALSE}
# following vignette for decontam package
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

sd <- as.data.frame(sample_data(ps)) # put sample_data into ggplot friendly data frame

sd$librarysize <- sample_sums(ps) #generate librarysize column

plot(sd$librarysize) # plot library size
# data were rarefied in QIIME2 analysis so librarysize plot is straight line
```
##### Decontaminate taxa using *decontam*
```{r, message = FALSE}
# package 'decontam' identified contaminants using frequency and prevalence 
# and are binned as contaminants if either method identifies it as such

sample_data(ps)$is.neg <- sample_data(ps)$sample_or_control == "control" # utilizes sample_or_control column from metadata where samples are classified as either patient sample or control

decon <- isContaminant(ps,
                        method = "either", # specify method
                        neg = "is.neg", # identify negative samples - logical
                        conc = "quant", # column specifying input DNA concentration
                        threshold = 0.5, # probability threshold
                        normalize = FALSE) # no normalization - samples are rarefied

table(decon$contaminant) #370 contaminants

contams <- decon[decon$contaminant == TRUE,] # filter decon object to only contaminants

contams$asv <- rownames(contams) # generate asv column from rownames

taxonomyTable <- as.data.frame(taxonomy) # convert taxonomy matrix to dataframe

taxonomyTable$asv <- rownames(taxonomyTable) # generate asv column from rownames

contams <- left_join(contams,
                     taxonomyTable,
                     by = "asv") # join taxonomyTable to contams as record of contaminants

noncontams <- decon[decon$contaminant == FALSE,] # filter decon object to only non-contaminants

noncontams$asv <- rownames(noncontams) # generate asv column from rownames

noncontams <- left_join(noncontams,
                        taxonomyTable,
                        by = "asv") # join taxonomyTable to noncontams as record of non-contaminants

filter <- decon$contaminant == TRUE | rowSums(asvTable) <= 10 # create logical filter that filters out contaminants and low abundance taxa 

filter <- as.logical(filter) # convert filter to logical

summary(filter) # summarise filter values

psFiltered <- prune_taxa(!filter, 
                          ps) # prune contaminants and low abundance taxa from phyloseq object

psFiltered
```
##### Filter low DNA input samples in *phyloseq*
```{r, message = FALSE}
sampFilter <- (metadata$sample_or_control == "sample" & metadata$dna_conc < 0.01) # create logical to filter samples (but not controls) that have low input DNA measurements

psSampFiltered <- prune_samples(!sampFilter,
                                psFiltered) # filter low DNA samples

psSampFiltered
```
##### Generate tree with *ape*
```{r}
randomTree <- rtree(ntaxa(psSampFiltered),
                    rooted = TRUE,
                    tip.label = taxa_names(psSampFiltered)) # generate random tree

plot_tree(randomTree,
          label.tips = NULL) # plot tree

psSampFiltered <- merge_phyloseq(psSampFiltered,
                                 randomTree)

psSampFiltered
```
##### Filter out non-Bacterial taxa and controls in *phyloseq*
```{r}
psBact <- subset_taxa(psSampFiltered, 
                      kingdom == "Bacteria") # remove non bacterial taxa

psBact # check dimensions

psBact <- subset_samples(psBact,
                         sampleloc != "control") # remove controls from dataset

psBact # check dimensions
```
##### Generate phylogenetic tree with *ape*
```{r, message = FALSE}
bactTree <- rtree(ntaxa(psBact),
                  rooted = TRUE,
                  tip.label = taxa_names(psBact)) # generate random tree

plot_tree(bactTree,
          label.tips = NULL) # plot tree

psBact <- merge_phyloseq(psBact,
                         bactTree)

psBact
```
##### Merge duplicate genera together in *phyloseq*
```{r}
psBactgl <- tax_glom(psBact,
                     "genus",
                     NArm = FALSE)

psBactgl
```
##### Generate dataframe from *phyloseq* object
```{r, message = FALSE}
allDf <- psmelt(psBactgl) # use psmelt to generate dataframe

allDf <- allDf %>%
  dplyr::select(-c(Sample,
                   is.neg)) # remove columns generated in phyloseq filtering 

write.csv(allDf,
          file.path(wd,
                    "16S",
                    "all_sample_dataframe.csv"),
          row.names = FALSE) # write to file 
```
##### Tidy genera with multiple hyphenated names
```{r, message = FALSE}
allDf$genus <-  sapply(strsplit(allDf$genus,"-"), `[`, 1)
```
##### Extract unassigned genera to *BLAST* 16S amplicons for identification
```{r, message = FALSE}
taxaUnkn <- allDf %>% 
  filter_all(any_vars(. == "Unassigned")) # filter to any taxa with Unassigned values

taxaUnkn <- taxaUnkn %>%
  select(OTU,kingdom,phylum,class,order,family,genus)

taxaUnkn <- unique(taxaUnkn)

write.table(x = taxaUnkn$OTU, 
            file = file.path(wd,
                             "16S",
                             "unassigned_taxa.txt"), 
            sep = "\n", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

taxaUnkn
```
##### Export representative sequences to FASTA file using *QIIME2*
```{bash, message = FALSE, results = FALSE}
qiime tools export \
--input-path qiime/V4/uchime-out/rep-seqs-nonchimeric-wo-borderline.qza \
--output-path qiime/V4/sequences-no-chi.fasta
```
##### Use *grep* to pull out ASVs from unassigned taxa for *BLAST*
```{bash}
grep -A 1 -f adult_CF_CRS/16S/unassigned_taxa.txt adult_CF_CRS/16S/sequences-no-chi.fasta > adult_CF_CRS/16S/unassigned_taxa_BLAST.txt

# results saved adult_CF_CRS/16S/C2DTRGSE013-Alignment.txt
```
##### Manually replace genera of unknown ASVs
```{r, message = FALSE}
allDf <- allDf %>%
  mutate(genus = case_when(OTU == "5f29d043bc53ad90a02d787d8798a4d2" ~ "Achromobacter",
                           OTU == "ac50f0d54ca5ff006c20f3052c5c4a10" ~ "Haemophilus",
                           OTU == "48c5e222f357b39de396cda86d2abb7b" ~ "Bradyrhizobium",
                           OTU == "fd661f32753c511bb0bb7160c9ca786b" ~ "Eikenella",
                           OTU == "67499d31f8b9c0f847e46c4ccd9d47df" ~ "Acidovorax",
                           OTU == "db22fe162ae9252b4e298875316a7748" ~ "Granulicatella",
                           OTU == "c7d69f1bab9c945505272848b78a0eaf" ~ "Selemonas",
                           OTU == "3d0cca66b082ac0333b14cee32616a72" ~ "Klebsiella",
                           OTU == "0a8045bfff248a724207b9459fecfee7" ~ "Tessaracoccus",
                           OTU == "257d01490aaff13da6c7375a67cc2d5a" ~ "Pandoraea",
                           OTU == "15c824f83de849f7a9f4cd155e104bee" ~ "Arthrobacter",
                           TRUE ~ genus))

allDf <- allDf %>%
  filter(OTU != "56a649f237e8ef55982e6b5662561545") # BLAST identified as human mitochondrial DNA - removed from dataset 

allDf <- allDf %>%
  filter(phylum != "Cyanobacteria") # remove all ASVs identified as cyanobacteria chloroplasts
```
##### Generate relative abundance column
```{r, message = FALSE}
allDf <- allDf %>%
  group_by(index) %>%
  mutate(rel_ab = (Abundance/sum(Abundance))*100) %>%
  ungroup()
```
##### Write to file
```{r, message = FALSE}
write.csv(allDf,
          file.path(wd,
                    "16S",
                    "all_sample_dataframe_tidy.csv"),
          row.names = FALSE)
```
##### Convert tidied dataframe back to *phyloseq* obejcts for alpha diversity calculations
```{r, message = FALSE}
taxonomyTidy <- allDf %>%
  dplyr::select(OTU,kingdom,phylum,class,order,family,genus) # select columns for tidy taxonomy table

taxonomyTidy <- distinct(taxonomyTidy) # remove duplicate rows

taxonomyTidy <- as.data.frame(taxonomyTidy) # convert to dataframe to allow rownames to be set

rownames(taxonomyTidy) <- taxonomyTidy[,1] # set ASVs as rownames

taxonomyTidy <- taxonomyTidy[,-1] # remove ASV column

taxonomyTidy <- as.matrix(taxonomyTidy)

taxTidy <- tax_table(taxonomyTidy) # generate taxonomy table object

asvTableTidy <- allDf %>%
  select(OTU,index,Abundance) # select colummns for tidy ASV table

asvTableTidy <- pivot_wider(asvTableTidy, 
                            names_from = "index", 
                            values_from = "Abundance") # pivot table to wide format

asvTableTidy <- as.data.frame(asvTableTidy) # convert to dataframe to allow rownames to be set

rownames(asvTableTidy) <- asvTableTidy[,1] # set ASVs as rownames

asvTableTidy <- asvTableTidy[,-1] # remove ASV column

otuTidy <- otu_table(asvTableTidy, taxa_are_rows = TRUE) # generate otu table object

sampTidy <- allDf %>%
  select(index,dna_conc,patientid,sampleloc,date,virus,virus_2,hemt,sample_or_control,virus_yn,quant) # select columns for tidy sample data table

sampTidy <- distinct(sampTidy) # remove duplicate rows

sampTidy <- as.data.frame(sampTidy) # set as dataframe to allow rownames to be set

rownames(sampTidy) <- sampTidy[,1] # set index column as rownames

sampTidyps <- sample_data(sampTidy) # generate sample data object

psTidy <- merge_phyloseq(otuTidy,taxTidy,sampTidyps) # merge tidied phyloseq objects into phyloseq object

sample_data(psTidy)
psTidy
```
##### Estimate richness with *phyloseq*
```{r, message = FALSE}
allAlpha <- estimate_richness(psTidy, 
                              split = TRUE, 
                              c("Observed","Chao1","Shannon","Simpson"))

allAlpha$index <- rownames(allAlpha)
```
##### Write to file
```{r, message = FALSE}
write.csv(allAlpha,
          file.path(wd,
                    "16S",
                    "all_sample_alpha_diversity_tidy.csv"),
          row.names = FALSE)
```
##### Import data
<h6>Data acquired from viewing "demux-viz.qzv" file in *QIIME* web viewer</h6>
```{r, message = FALSE}
exos <- read_tsv(file.path(wd,
                           "primer_data",
                           "exoS_per_sample_counts.tsv"))
exou <- read_tsv(file.path(wd,
                           "primer_data",
                           "exoU_per_sample_counts.tsv"))
gsea <- read_tsv(file.path(wd,
                           "primer_data",
                           "gseA_per_sample_counts.tsv"))
ldh  <- read_tsv(file.path(wd,
                           "primer_data",
                           "ldh1_per_sample_counts.tsv"))
meca <- read_tsv(file.path(wd,
                           "primer_data",
                           "mecA_per_sample_counts.tsv"))
pela <- read_tsv(file.path(wd,
                           "primer_data",
                           "pelA_per_sample_counts.tsv"))
pf   <- read_tsv(file.path(wd,
                           "primer_data",
                           "Pf5_per_sample_counts.tsv"))
psla <- read_tsv(file.path(wd,
                           "primer_data",
                           "pslA_per_sample_counts.tsv"))
tset <- read_tsv(file.path(wd,
                           "primer_data",
                           "tseT_per_sample_counts.tsv"))
tsit <- read_tsv(file.path(wd,
                           "primer_data",
                           "tsiT_per_sample_counts.tsv"))
```
##### Create preproc function to process individual primer files
```{r, message = FALSE}
preproc <- function(primer) {
  primer$'sample ID' <- str_replace_all(primer$'sample ID',c("\\_2018" = "\\_18","\\_2019" = "\\_19",
                                                   "\\_2020" = "\\_20","\\_2021" = "\\_21","\\_2022" = "\\_22"))
  primer <- primer %>% select('sample ID','forward sequence count')
  primer <- primer %>% rename(index = sample,
                              count = 'forward sequence count')
}
```
```{r, message = FALSE}
preproc <- function(primer) {
  primer$index <- str_replace_all(primer$'sample ID',
                                c("\\_2018" = "\\_18",
                                  "\\_2019" = "\\_19",
                                  "\\_2020" = "\\_20",
                                  "\\_2021" = "\\_21",
                                  "\\_2022" = "\\_22"))
  
  primer <- primer %>%
    dplyr::select(c(2,4))
  
  primer <- primer %>%
    rename(index = 2,
           count = 1)
  
  primer <- primer %>%
    relocate(index,
             .before = count)
}
```
##### Use preproc function to process imported data
```{r, message = FALSE}
exos <- preproc(exos)
exou <- preproc(exou)
gsea <- preproc(gsea)
ldh  <- preproc(ldh)
meca <- preproc(meca)
pela <- preproc(pela)
pf   <- preproc(pf)
psla <- preproc(psla)
tset <- preproc(tset)
tsit <- preproc(tsit)
```
##### Join primer dataframes together
```{r, message = FALSE}
primerCounts <- full_join(exos,exou,by="index") %>%
  rename(exos = 2,exou = 3) %>%
  full_join(.,gsea,by="index") %>%
  rename(gsea = 4) %>%
  full_join(.,ldh,by="index") %>%
  rename(ldh = 5) %>%
  full_join(.,meca,by="index") %>%
  rename(meca = 6) %>%
  full_join(.,pela,by="index") %>%
  rename(pela = 7) %>%
  full_join(.,pf,by="index") %>%
  rename(pf = 8) %>%
  full_join(.,psla,by="index") %>%
  rename(psla = 9) %>%
  full_join(.,tset,by="index") %>%
  rename(tset = 10) %>%
  full_join(.,tsit,by="index") %>%
  rename(tsit = 11)
```
##### Fix mislabeled sample index in primerCounts
```{r, message = FALSE}
gsub("GESCF13_sinus_3_8_17",
     "GESCF13_sinus_3_8_18",
     primerCounts$index) -> primerCounts$index
```
##### Create binary dataframe of positive/negative results from primerCounts dataframe based on 1% maximum read count cut-off for positivity
```{r, message = FALSE}
primerBinary <- primerCounts %>%
  mutate(exos = case_when(exos >= (max(exos,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(exos) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(exou = case_when(exou >= (max(exou,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(exou) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(gsea = case_when(gsea >= (max(gsea,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(gsea) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(ldh  = case_when(ldh >= (max(ldh,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(ldh) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(meca = case_when(meca >= (max(meca,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(meca) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(pela = case_when(pela >= (max(pela,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(pela) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(pf   = case_when(pf >= (max(pf,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(pf) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(psla = case_when(psla >= (max(psla,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(psla) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(tset = case_when(tset >= (max(tset,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(tset) ~ "NA",
                          TRUE ~ "0")) %>%
  mutate(tsit = case_when(tsit >= (max(tsit,
                                       na.rm = TRUE)*0.01) ~ "1",
                          is.na(tsit) ~ "NA",
                          TRUE ~ "0"))

primerBinary[primerBinary == "NA"] <- NA
```
##### Write data to file
```{r, message = FALSE}
write.csv(primerCounts,
          file.path(wd,
                    "primer_data",
                    "all_primer_counts.csv"),
          row.names = FALSE)

write.csv(primerBinary,
          file.path(wd,
                    "primer_data",
                    "all_primer_binary.csv"),
          row.names = FALSE)
```
##### Generate standard curve for 16S 515f - 806r universal primers for qPCR copy number calculation
```{r, message = FALSE}
cq <- read.csv(file.path(wd,
                         "qPCR",
                         "16S_515f806r_mod_standard_curve_2023-03-29_Cq.csv")) #read in Cq data for standard curve

wells <- c("A01","B01","C01","D01","E01","F01","G01","H01", #create vector of wells used in plate
           "A02","B02","C02","D02","E02","F02","G02","H02",
           "A03","B03","C03","D03","E03","F03","G03","H03",
           "A04","B04","C04","D04")

samples <- c(1.17e+11,1.17e+11,1.17e+10,1.17e+10,1.17e+9,1.17e+9,1.17e+8,1.17e+8, #create vector of samples in wells 
             1.17e+7,1.17e+7,1.17e+6,1.17e+6,1.17e+5,1.17e+5,1.17e+4,1.17e+4,
             1.17e+3,1.17e+3,1.17e+2,1.17e+2,1.17e+1,1.17e+1,1.17,1.17,
             1.17e-1,1.17e-1,0,0)    

plate <- data.frame(wells,samples) #join wells and samples to create plate df

cq <- cq %>%
  rename(well = Well) #rename Well column

plate <- plate %>%
  rename(well = wells) #rename wells column

cq <- left_join(plate,
                cq,
                by = "well") #join plate and cq to keep only data from wells with samples

cq <- cq %>%
  select(well,samples,Cq) #select only desired columns

cq[cq=="NaN"] <- NA #change NaN generated by qPCR machine to NA
cq[cq==0] <- NA #change 0s to NA

cq <- cq %>%
  mutate(log10_cn = log10(samples)) #create log10 dilution column for plotting - must plot log10 serial dilution in order to generate straight line plot

cq_nolow <- cq %>%
  filter(!is.na(samples) & !is.na(Cq)) %>%
  filter(samples != 1.17e-01) #create df with no NA samples and without dilution with fewer than 1 copy of target gene

lm(cq_nolow$log10_cn ~ cq_nolow$Cq) #generate regression

ggplot(cq,
       aes(x = Cq,
           y = log10_cn,
           groups = log10_cn)) + #plot cq against log10 copy number to generate standard curve
  geom_point() + #plot points
  geom_abline(intercept=9.6124,
              slope=-0.2853,
              color="red") #plot regression line
```
##### Calculate 16S copy number per sample from regression equation
```{r, message = FALSE}
locs <- read.csv(file.path(wd,
                           "qPCR",
                           "dna_dilutions_v2.csv")) #import dna location info

locs <- locs[1:391,] #extra rows added in excel - remove

locs$well1 <- str_replace_all(locs$well1,"^1$","01") %>% #replace single digit number with two digit format
  str_replace_all(.,"^3$","03") %>%
  str_replace_all(.,"^5$","05") %>%
  str_replace_all(.,"^7$","07") %>%
  str_replace_all(.,"^9$","09") %>%
  str_replace_all(.,"^2$","02") %>% 
  str_replace_all(.,"^4$","04") %>%
  str_replace_all(.,"^6$","06") %>%
  str_replace_all(.,"^8$","08")

locs$well2 <- str_replace_all(locs$well2,"^2$","02") %>% #replace single digit number with two digit format
  str_replace_all(.,"^4$","04") %>%
  str_replace_all(.,"^6$","06") %>%
  str_replace_all(.,"^8$","08") %>%
  str_replace_all(.,"^1$","01") %>% 
  str_replace_all(.,"^3$","03") %>%
  str_replace_all(.,"^5$","05") %>%
  str_replace_all(.,"^7$","07") %>%
  str_replace_all(.,"^9$","09")


locs <- locs %>%
  mutate(df = if_else(dna_input == 10,1,(10/dna_input))) #add dilution factor column

locs <- locs %>%
  mutate(well1 = paste(row,well1,sep="")) %>% #paste together row numbers and letters to form well IDs
  mutate(well2 = paste(row,well2,sep="")) %>%#paste together row numbers and letters to form well IDs
  pivot_longer(cols=c(well1,well2),names_to="well_no",values_to="well") %>% #pivot to create one column for well IDs instead of two
  select(sample_number,box,cell,sampleloc,patientid,date,df,plate,well) #select to remove extra well info columns

#formatting in excel means some repeated wells are labeled incorrectly - need to replace
locs[101,8] <- 9
locs[101,9] <- "D09"
locs[142,9] <- "E01"
locs[184,8] <- 9
locs[184,9] <- "E02"
locs[295,8] <- 9
locs[295,9] <- "E11"
locs[300,9] <- "F01"
locs[696,9] <- "G01"

locs$date <- as.Date(locs$date,format="%m/%d/%y") #convert dates to date format

#import plate data
p1_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_one_2023-04-04_Cq.csv"))
p2_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_two_2023-04-05_Cq.csv"))
p3_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_three_2023-04-05_Cq.csv"))
p4_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_four_2023-04-06_Cq.csv"))
p5_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_five_2023-04-06_Cq.csv"))
p6_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_six_2023-04-06_Cq.csv"))
p7_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_seven_2023-04-07_Cq.csv"))
p8_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_eight_2023-04-07_Cq.csv"))
p9_cq <- read.csv(file.path(wd,
                            "qPCR",
                            "16S_515f806r_mod_GILEAD_plate_nine_2023-04-10_Cq.csv"))

#add plate info column
p1_cq <- p1_cq %>%
  mutate(plate = 1)
p2_cq <- p2_cq %>%
  mutate(plate = 2)
p3_cq <- p3_cq %>%
  mutate(plate = 3)
p4_cq <- p4_cq %>%
  mutate(plate = 4)
p5_cq <- p5_cq %>%
  mutate(plate = 5)
p6_cq <- p6_cq %>%
  mutate(plate = 6)
p7_cq <- p7_cq %>%
  mutate(plate = 7)
p8_cq <- p8_cq %>%
  mutate(plate = 8)
p9_cq <- p9_cq %>%
  mutate(plate = 9)


cq <- rbind(p1_cq,p2_cq,p3_cq,p4_cq,p5_cq,p6_cq,p7_cq,p8_cq,p9_cq) %>% #join plates together
  select(Well,Cq,plate) #select only relevant columns

cq[cq=="NaN"] <- NA #change NaN to NA

cq <- cq %>%
  rename("cq"="Cq",
         "well"="Well") #rename columns to allow joining

na_cq <- cq %>%
  filter(is.na(cq)) %>% #transformations for cq will not work on NAs so remove to separate df for now
  mutate(log10_cn = 0) %>%
  mutate(cn = 0)

nona_cq <- cq %>%
  filter(!is.na(cq)) #no NA df

nona_cq <- nona_cq %>%
  mutate(log10_cn = (-0.2853*cq)+9.6124) %>% #calculate log10 copy no from standard curve slope and intercept
  mutate(cn = 10^log10_cn) #antilog to calculate copy no

cq <- rbind(nona_cq,na_cq) #join NA and non NA dfs back together

cq_locs <- left_join(cq,
                     locs,
                     by = c("plate","well"))

cq_locs <- cq_locs %>%
  filter(!is.na(sample_number)) #filter out NAs 

cq_locs <- cq_locs %>%
  mutate(cn_df = cn*df) #multiply copy number by dilution factor

cq_locs$date <- as.Date(cq_locs$date,
                        format = "%m/%d/%y") #format as date for join

cq_locs$sampleloc <- str_replace_all(cq_locs$sampleloc,"Sinus","sinus") %>% #remove capital letters from sampleloc for join
  str_replace_all(.,"Sputum","sputum") %>%
  str_replace_all(.,"Throat","throat")

allDfTidy %>%
  select(patientid,sampleloc,date,hemt) %>%
  distinct() -> meta

meta$date <- as.Date(meta$date)

cq_locs_meta <- left_join(cq_locs,
                          meta,
                          by = c("patientid","sampleloc","date")) #join with metadata

cq_locs_meta$hemt <- factor(cq_locs_meta$hemt,
                            levels =c ("pre","post","not prescribed")) #set hemt column as factor

cq_locs %>% group_by(sample_number) %>%
  mutate(mean_cn = mean(cn_df)) %>%
  ungroup() %>%
  select(patientid,sampleloc,date,mean_cn) -> cn

distinct(cn) -> cn
```
##### Write data to file
```{r, message = FALSE}
write.csv(cn,
          file.path(wd,
                    "qPCR",
                    "all_sample_16S_copy_number.csv"),
                    row.names = FALSE)
```
#### Sinus data analysis
##### Import data
```{r, messsage = FALSE}
allDf <- read.csv(file.path(wd,
                            "16S",
                            "all_sample_dataframe_tidy.csv"))
                  
allAlpha <- read.csv(file.path(wd,
                               "16S",
                               "all_sample_alpha_diversity_tidy.csv"))

allCn <- read.csv(file.path(wd,
                            "qPCR",
                            "all_sample_16S_copy_number.csv"))
```
##### Join meta data tables
```{r, message = FALSE}
allDfTidy$date <- as.Date(allDfTidy$date)

allDfTidy %>%
  select(index,patientid,sampleloc,date,hemt) %>%
  distinct() %>%
  left_join(.,
            meta,
            by = c("patientid","sampleloc","date","hemt")) -> meta

left_join(allAlpha,
          meta,
          by = "index") -> allAlpha
```
##### Filter data to only sinus samples
```{r, message = FALSE}
sinusDf <- allDf %>%
  filter(sampleloc == "sinus")

sinusAlpha <- allAlpha %>%
  filter(grepl("sinus",index))

sinusCn <- allCn %>%
  filter(sampleloc == "sinus")
```
##### Ensure data columns are correctly formatted
```{r, message = FALSE}
as.Date(sinusDf$date) -> sinusDf$date

as.Date(sinusAlpha$date) -> sinusAlpha$date

as.Date(sinusCn$date) -> sinusCn$date
```
##### Join dataframes together
```{r, message = FALSE}
sinusDfCn <- left_join(sinusDf,
                       sinusCn,
                       by = c("patientid","sampleloc","date"))

sinusDfCnAlpha <- left_join(sinusDfCn,
                            sinusAlpha,
                            by = c("index","patientid","sampleloc","date","hemt"))
```
##### Filter out samples with no copies of 16S by qPCR
```{r, message = FALSE}
sinusDfCnAlpha <- sinusDfCnAlpha %>%
  filter(mean_cn != 0)
```
##### Format data for hierarchical clustering
```{r, message = FALSE}
sinusHemtPID <- c("GESCF04","GESCF05","GESCF08","GESCF09","GESCF12","GESCF14","GESCF16","GESCF24",
                  "GESCF26","GESCF28","GESCF33","GESCF36","GESCF37","GESCF39") # generate vector of patient IDs of interest for ETI comparison

etiClust <- sinusDfCnAlpha %>%
  filter(patientid %in% sinusHemtPID) %>%
  dplyr::select(index,genus,rel_ab) # select relevant columns for clustering

etiMeta <- sinusDfCnAlpha %>%
  filter(patientid %in% sinusHemtPID) %>%
  dplyr::select(index,patientid,hemt) %>%
  distinct() # generate metadata df

etiClustWide <- pivot_wider(etiClust,
                             names_from = genus,
                             values_from = rel_ab,
                             values_fn = mean) # pivot wider for clustering

etiClustWide <- etiClustWide %>%
  column_to_rownames(var = "index") # remove index column and assign to rownames

etiClustWide[is.na(etiClustWide)] <- 0 # replace NAs with 0s
```
##### Generate Morisita-Horn distance matrix using *vegan*
```{r, message = FALSE}
etiMh <- vegdist(etiClustWide,
                 method = "horn")
```
##### Perform hierarchical clustering
```{r, message = FALSE}
etiDend <- etiClustWide %>%
  vegdist(method = "horn") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram()

plot(etiDend)

etiHc <- etiClustWide %>%
  vegdist(method = "horn") %>%
  hclust(method = "ward.D2")
```
##### Calculate centroid distances per patient
```{r, message = FALSE}
etiDat <- sinusDfCnAlpha %>%
  filter(patientid %in% sinusHemtPID) %>%
  dplyr::select(index,patientid,hemt) %>%
  distinct()
  
etiCent <- betadisper(etiMh,
                      group = as.factor(etiDat$patientid),
                      type = "centroid")

plot(etiCent)
boxplot(etiCent)

etiCentPt <- as.data.frame(etiCent$distances)

etiCentPt$index <- rownames(etiCentPt)

etiDfCnAlpha <- sinusDfCnAlpha %>%
  right_join(.,
             etiCentPt,
             by = "index")

etiDfCnAlpha <- etiDfCnAlpha %>%
  rename("centroid_distance" = `etiCent$distances`)
```
##### Cut dendrogram and assign cluster numbers
```{r, message = FALSE}
etik3 <- cutree(etiHc,
                k = 3)

etik3 <- as.data.frame(etik3)

etik3$index <- rownames(etik3)

etiDfCnAlpha <- left_join(etiDfCnAlpha,
                          etik3,
                          by = "index")

etiDfCnAlpha <- etiDfCnAlpha %>%
  rename("k3" = etik3)

etiDfCnAlpha$k3 <- factor(etiDfCnAlpha$k3)
```
##### Generate column indicating whether patient samples change cluster over time
```{r, message = FALSE}
etiChange <- etiDfCnAlpha %>%
  count(patientid,k3) %>%
  count(patientid) %>%
  mutate(change = if_else(n == 1,
                          "no change",
                          "change")) %>%
  select(patientid,change)

left_join(etiDfCnAlpha,
          etiChange,
          by = "patientid") -> etiDfCnAlpha
```
##### Perform linear regression to determine if increased distance to centroid per sample is associated with patients with samples from multiple clusters
```{r, message = FALSE}
etiDfCnAlpha %>%
  select(index,patientid,centroid_distance,change) %>%
  distinct() -> etiClusChange

etiCClm <- lmer(data = etiClusChange, centroid_distance ~ change + (1 | patientid))

summary(etiCClm)
```
##### Perform linear regression to determine if Shannon diversity index is significantly altered post-ETI
```{r, message = FALSE}
sinusAlphaNnp <- sinusAlpha %>%
  filter(patientid %in% sinusHemtPID)

sinusAlphaLm <- lmer(data = sinusAlphaNnp,
                Shannon ~ hemt + (1 | patientid))

summary(sinusAlphaLm)
```
##### Perform linear regression to determine if total bacterial abundance is significantly altered post-ETI
```{r, message = FALSE}
sinusCnDis <- sinusDfCn %>%
  dplyr::select(index,patientid,hemt,mean_cn) %>%
  filter(patientid %in% sinusHemtPID) %>%
  distinct()

sinusCnLm <- lmer(data = sinusCnDis,
                  mean_cn ~ hemt + (1 | patientid))

summary(sinusCnLm)
```
##### Prepare data for MaAsLin 2 analysis
```{r, message = FALSE}
sinusDfCnAlpha %>% 
  filter(patientid %in% sinusHemtPID) %>%
  select(index,genus,rel_ab) %>%
  pivot_wider(names_from = genus, values_from = rel_ab, values_fn = mean) -> counts

column_to_rownames(counts, var = "index") -> counts

sinusDfCnAlpha %>%
  filter(patientid %in% sinusHemtPID) %>%
  select(index,patientid,hemt) %>%
  distinct() -> meta

column_to_rownames(meta, var = "index") -> meta
```
##### Perform MaAsLin 2 analysis
```{r, collapse = TRUE, eval = FALSE}
Maaslin2(input_data = counts,
         input_metadata = meta,
         min_abundance = 0.0001,
         min_prevalence = 0.01,
         analysis_method = "LM",
         normalization = "NONE",
         output = wd,
         fixed_effects = "hemt",
         random_effects = "patientid")
```

***

#### Sinus and sputum comparison analysis  
##### Generate list of patients and dates with paired samples
```{r, message = FALSE}
paired <- allDf %>%
  filter(hemt != "post") %>%
  filter(sampleloc != "throat") %>%
  dplyr::select(index,patientid,sampleloc,date) %>%
  distinct() %>%
  count(patientid,date) %>% # count how many occurrences of each patient and sample date
  filter(n == 2) # filter to only patientids and sample dates with 2 samples (indicates sinus and sputum samples)

paired <- paired %>%
  mutate(id = paste(patientid,date,
                    sep = "_"),
         .before=patientid)

paired <- paired$id

pairedDf <- allDf %>%
  mutate(id = paste(patientid,date,
                    sep = "_"),
         .after = index)

pairedDf <- pairedDf %>%
  filter(hemt != "post") %>%
  filter(sampleloc != "throat") %>%
  filter(id %in% paired)

pairedData <- pairedDf %>%
  select(index,id,patientid,sampleloc) %>%
  distinct()

pairedIndex <- pairedData$index

pairedAlpha <- allAlpha %>%
  filter(index %in% pairedIndex)

pairedAlpha <- left_join(pairedAlpha,
                         pairedData,
                         by = "index")

allCn$date <- as.Date(allCn$date)

pairedDf$date <- as.Date(pairedDf$date)

pairedDfCn <- left_join(pairedDf,
                        allCn,
                        by = c("patientid","date","sampleloc"))
```
##### Calculate concordance metrics
```{r, message = FALSE}
# sensitivity = true_positive / (true_positive + false_negative) [AKA recall]
# specificity = true_negative / (true_negative + false_positive)
# precision   = true_positive / (true_positive + false_positive)
# accuracy    = (true_positive + true_negative) / total
# F measure   = (2 * sensitivity * precision) / (sensitivity + precision)
# true positive: non-zero abundance in both sinus and sputum
# true negative: 0 abundance in both sinus and sputum
# false positive: non-zero abundance in sinus, 0 abundance in sputum
# false negative: 0 abundance in sinus, non-zero abundance in sputum

pairedTNTPFNFP <- pairedDfCn %>%
  dplyr::select(id,patientid,sampleloc,date,OTU,genus,rel_ab) %>%
  pivot_wider(names_from = sampleloc, values_from = rel_ab) %>%
  mutate(true_pos = if_else(sinus > 0 & sputum > 0,1,0)) %>%
  mutate(true_neg = if_else(sinus == 0 & sputum == 0,1,0)) %>%
  mutate(false_pos = if_else(sinus > 0 & sputum == 0,1,0)) %>%
  mutate(false_neg = if_else(sinus == 0 & sputum > 0,1,0)) %>%
  group_by(OTU) %>%
  mutate(true_posSum = sum(true_pos,na.rm = TRUE)) %>%
  mutate(true_negSum = sum(true_neg,na.rm = TRUE)) %>%
  mutate(false_posSum = sum(false_pos,na.rm = TRUE)) %>%
  mutate(false_negSum = sum(false_neg,na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(genus,true_posSum,true_negSum,false_posSum,false_negSum) %>%
  distinct() 
  

pairedSensSpec <- pairedTNTPFNFP %>%
  group_by(genus) %>%
  summarise(true_posSum = mean(true_posSum),
            true_negSum = mean(true_negSum),
            false_posSum = mean(false_posSum),
            false_negSum = mean(false_negSum)) %>%
  ungroup() %>%
  mutate(sensitivity = true_posSum / (true_posSum + false_negSum)) %>%
  mutate(specificity = true_negSum / (true_negSum + false_posSum)) %>%
  mutate(precision = true_posSum / (true_posSum + false_posSum)) %>%
  mutate(accuracy = (true_posSum + true_negSum) / (true_posSum + true_negSum + false_posSum + false_negSum)) %>%
  mutate(F_measure = (2 * sensitivity * precision) / (sensitivity + precision)) %>%
  select(genus,sensitivity,specificity,precision,accuracy,F_measure)

pairedProps <- pairedTNTPFNFP %>%
  mutate(tpProp = true_posSum/53) %>%
  mutate(tnProp = true_negSum/53) %>%
  mutate(fpProp = false_posSum/53) %>%
  mutate(fnProp = false_negSum/53) %>%
  pivot_longer(cols = c(tpProp,tnProp,fpProp,fnProp),
               names_to = "meas", 
               values_to = "prop") %>%
  dplyr::select(genus,meas,prop)

pairedProps$meas <- factor(pairedProps$meas,
                           levels = c("tpProp","tnProp","fpProp","fnProp"))
```
##### Print session info
```{r, echo = TRUE, results = TRUE}
sessionInfo()
```
