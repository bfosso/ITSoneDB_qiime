# ITSoneDB_qiime
[**ITSoneDB**](http://itsonedb.cloud.ba.infn.it) is a curated database [1] collecting the **Internal transcribed spacer 1** regions of the *Eukaryotic* ribosomal gene clusters.  
This repository contains the resources for the implemetation of ITSoneDB in the [**QIIME2**](https://qiime2.org) pipeline [2]. 

# Sequence clustering
ITS1 sequences available in the ITSoneDB were collected in the `ITSoneDB_total_fasta_rel138.fa` fasta file and clustered at 99% of similarity by using **VSEARCH** [3].
```
vsearch --cluster_fast ITSoneDB_total_fasta_rel138.fa -uc vsearch_cluster_99 \
--centroids 99_ref_ITSoneDB.fa --id 0.99
```
# Taxonomic annotation
![tree_figure](./tree.png) 

***Figure1:*** *A representation of a list of nodes represented belonging to a cluster in a reference taxonomy.* 

The obtained cluters represenative sequences were taxonomically labelled by using two approaches:  
* **Lowest common ancestor (LCA)**: this algorithm assigne a taxonomica label by considerign the lowest taxonomic rank common to all the sequences belonging to a cluster.   
* **TANGO** algorithm: this algorithm finds the most likely node by inferrig F-measures [4].  

# Data
Data available in this repository:  
1. [`ITSoneDB_total_fasta_rel138.fa.gz`](https://www.dropbox.com/s/0mozfmhgamq7ems/ITSoneDB_total_fasta_rel138.fa.gz?dl=0): compressed fasta file containing the ITSoneDB ITS1 sequences (use the link to download the file);  
2. `99_ref_ITSoneDB.fa.gz`: ITS1 clusters representative sequences;  
3. `tax_99_ref_ITSoneDB.txt.gz`: NCBI taxonomic path associated to the clusters representative sequences inferred by using the LCA algorithm;  
4. `tax_99_ref_ITSoneDB_tango_based.txt.gz`: NCBI taxonomic path associated to the clusters representative sequences inferred by using the TANGO algorithm;  

# QIIME2 import commands
Following are listed the QIIME2 command to import ITSoneDB data as `qza`artifacts:  
```
qiime tools import \
  --input-path 99_ref_ITSoneDB.fasta \
  --output-path 99_ref_ITSoneDB.qza \
  --type 'FeatureData[Sequence]'

qiime tools import \
 --type FeatureData[Taxonomy] \
 --input-path tax_99_ref_ITSoneDB.txt \
 --input-format HeaderlessTSVTaxonomyFormat \
 --output-path ITSoneDB_taxonomy.qza
```

# References
1. Santamaria M, Fosso B, Licciulli F, Balech B, Larini I, Grillo G, et al. ITSoneDB: a comprehensive collection of eukaryotic ribosomal RNA Internal Transcribed Spacer 1 (ITS1) sequences. Nucleic Acids Res. 2018;46: D127–D132. doi:10.1093/nar/gkx855
2. Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet C, Al-Ghalith GA, et al. QIIME 2: Reproducible, interactive, scalable, and extensible microbiome data science. PeerJ Inc.; 2018 Dec. Report No.: e27295v2. doi:10.7287/peerj.preprints.27295v2
3. Rognes T, Flouri T, Nichols B, Quince C, Mahe F. VSEARCH: a versatile open source tool for metagenomics. PeerJ. 2016;4: e2584. doi:10.7717/peerj.2584
4. Alonso-Alemany D, Barre A, Beretta S, Bonizzoni P, Nikolski M, Valiente G. Further Steps in TANGO: improved taxonomic assignment in metagenomics. Bioinformatics. 2014;30: 17–23. doi:10.1093/bioinformatics/btt256
