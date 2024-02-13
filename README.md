# Cystic fibrosis pathogens persist in the upper respiratory tract following initiation of elexacaftor/tezacaftor/ivacaftor therapy

This repository contains code and subset data files generated and used in analysis in the publication of this manuscript. Please contact Jennifer Bomberger (corresponding author; jbomb@dartmouth.edu) for other data availability.

### File list:
* <a href="https://github.com/yasminhilliam/sinus_ETI/blob/main/20240208_HILLIAM_data_analysis.md">20240208_HILLIAM_data_analysis.md</a>: Bash and R code used in the preparation of this manuscript.<br>
* <a href="https://github.com/yasminhilliam/sinus_ETI/blob/main/20240208_HILLIAM_sinus_data.csv">20240208_HILLIAM_sinus_data.csv</a>: OTU and abundance data table from subjects with pre- and post-ETI samples.<sup><span>&#8224;</span></sup><br>
* <a href="https://github.com/yasminhilliam/sinus_ETI/blob/main/20240208_HILLIAM_paired_data.csv">20240208_HILLIAM_paired_data.csv</a>:  OTU and abundance data table from subjects with paired sinus and sputum samples.<sup><span>&#8224;</span></sup>
<br>
<br>
<span>&#8224;</span> Data in these tables have been processed through QIIME2, phyloseq, and vegan and then filtered to only include samples from relevant subjects.
<br>
<hr>
<br>
### Data file headers:<br>
OTU: Amplicon sequence variant (ASV) assigned through QIIME2 and Silva database<br>
Abundance: Count per ASV<br>
index: Unique sample identifier (SUBJECTID_SAMPLELOCATION_SAMPLEDATE)<br>
id: Identified for paired samples (SUBJECTID_SAMPLEDATE) <span>&#42;</span><br>
sampleloc: Sampling location (sinus or sputum) <span>&#42;</span><br>
patientid: Anonymised subject identifier <span>&#167;</span><br>
date: Sample collection date (MM/DD/YY)<br>
hemt: Elexacaftor/tezacaftor/ivacaftor (ETI) status of subject at sample date <span>&#167;</span><br>
kingdom: Assigned kingdom<br>
phylum: Assigned phylum<br>
class: Assigned class<br>
order: Assigned order<br>
family: Assigned family<br>
genus: Assigned genus<br>
rel_ab: Calculated relative abundance (%)<br>
mean_cn: Mean 16S copy number from duplicate qPCR results<br>
Observed: Calculated Observed diversity <span>&#167;</span><br>
Chao1: Calculated CHAO1 diversity <span>&#167;</span><br>
Shannon: Calculated Shannon diversity <span>&#167;</span><br>
Simpson: Calculated Simpson diversity <span>&#167;</span><br>
<br>
<br>
<span>&#42;</span> 20240208_HILLIAM_paired_data.csv only<br>
<span>&#167;</span> 20240208_HILLIAM_paired_data.csv only

