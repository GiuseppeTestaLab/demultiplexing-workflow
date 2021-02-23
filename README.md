# demultiplexing-workflow

This Snakefile pipeline performs demultiplexing starting from Bamfile and VCF file. It will output demultiplexing results for:

- Demuxlet v1
- Demuxlet v2
- SoupOrCell
- Vireo
- SCanSNP

# Content overview

- Snakefile Main workflow for demultiplexing 
- scripts/Benchmarks_index3.0.0Mod_Plus_Synth.Rmd contains the benchmark script
- Sample* contains demultiplexing results for Demuxlet, Demuxletv2, Soup or Cell, Vireo and SCanSNP
- Sample*/Demultiplexing.report.html contains the result of scripts/Benchmarks_index3.0.0Mod_Plus_Synth.Rmd (sadly no preview is available because of size)

### SCanSNP output folder
- Diagnostic plot on DBLs call and LowQuality
- doubletsMarked.tsv with final ID assignment in the format 

| Barcode  | ID_1 | ID_2 | ID_N | DBL_FirstID_Score | DBL_SecondID_Score | FirstID | SecondID | DropletType | ID      | Qual        |
|----------|------|------|------|-------------------|--------------------|---------|----------|-------------|---------|-------------|
| AAATGC-1 | 262  | 2    | 10   | 1077.41           | 39.29              | ID1     | IDN      | Singlet     | ID1     | GoodQuality |
| AAAGGC-1 | 0    | 450  | 302  | 550.15            | 475.20             | ID2     | IDN      | Doublet     | Doublet | Doublet     |
| ACGTCC-1 | 0    | 150  | 1    | 70.00             | 50.00              | ID2     | IDN      | Singlet     | ID2     | LowQuality  |

In a nutshell: 

ID_*: (corrected) reads overlapping ID-private SNPs

ID: ID and doublet information 

Qual: Present only if  --flagLowQual set to True and contains Quality information (use carefully and looking at mixture fitting QC plot)


### Vireo output folder

- Column "donor_id" in file vireoOut/donor_ids.tsv  contains informations about type, and assignment


### SoC output folder

- Column "status" in file clusters.tsv  contains informations about type
- Column "assignment" in file clusters.tsv  contains informations about assignment


### Demuxlet V1

- Columns "BEST" in file *.best ations about type and assignment e.g: DBL-ID1
- Columns "SNG.1ST"in file *.best ontains informations about assignment assuming singlet type

### Demuxlet V2

- Column "DROPLET.TYPE" contains information about type
- Column "BEST.GUESS" contains information about most probable assignment (ordered pair) e.g. IDX,IDY,0.5 (use the first reported id)

### R packages



