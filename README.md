# HTG_Analysis
HTG analysis of Sarcoidosis, Cocci and TB

To be able to run this program you need to clone this repository into a linux machine.

You need to install the following software
ruby 2.6.6p146

Salmon version 1.3.0 (from https://github.com/COMBINE-lab/salmon/releases)



R version 3.6.3

BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so


attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] pheatmap_1.0.12    RColorBrewer_1.1-2 circlize_0.4.10    ggrepel_0.8.2      ggplot2_3.3.2      r2excel_1.0.0
[7] xlsx_0.6.3         edgeR_3.28.1       limma_3.42.2       hash_2.2.6.1       dplyr_1.0.0
[12] r2excel_1.0.0    optparse_1.6.6      tximport_1.14.2

loaded via a namespace (and not attached):
[1] Rcpp_1.0.5          pillar_1.4.6        compiler_3.6.3      tools_3.6.3         digest_0.6.25       evaluate_0.14
[7] lifecycle_0.2.0     tibble_3.0.3        gtable_0.3.0        lattice_0.20-41     pkgconfig_2.0.3     rlang_0.4.7
[13] rstudioapi_0.11     yaml_2.2.1          xfun_0.15           rJava_0.9-13        withr_2.2.0         knitr_1.29
[19] GlobalOptions_0.1.2 generics_0.0.2      vctrs_0.3.2         xlsxjars_0.6.1      locfit_1.5-9.4      tidyselect_1.1.0
[25] glue_1.4.1          R6_2.4.1            rmarkdown_2.3       farver_2.0.3        purrr_0.3.4         magrittr_1.5
[31] scales_1.1.1        ellipsis_0.3.1      htmltools_0.5.0     shape_1.4.4         colorspace_1.4-1    labeling_0.3
[37] munsell_0.5.0       crayon_1.3.4      jsonlite_1.7.0       hms_0.5.3    grid_3.6.3    getopt_1.20.3   readr_1.3.1   magrittr_1.5  crayon_1.3.4


After installing and testing that all the packages are working properly. You need to request the HTG files to ncasanova@arizona.edu. The package with the fastq files is fastqFiles.tar.bz2 with a md5sum of 1cfc2d7943ea1327d69bba84d8194c55. You will need to unpack the package by running tar -xvf fastqFiles.tar.bz2 into a directory. The same directory should contain the scripts directory from this repository. The scripts directory should contain the following files:

e5118d0352675df1fd9eb969213c9810  analyzeCounts.R

cb3140240c74a9e0a4baba56524c4fa3  generateRFiles.sh

24bc6b02c6f4fc1e735b42f49ceafa0c  generateTemplate.sh

8022b5f7b75eea6df9c71b0c488037a6  htgDataProcessing.R

c7680eeaddcfdbeddb8bbdbb56f8233b  mainScript.sh

2b97090783ea69f342168daafdbb30d4  minRNASeqFuncLib.R

92a98050e45ff7e7ac76593b01eff2f0  reNumberFile.rb

The main directory should contain the scripts and fastqFiles folders. Go to the main directory and run the following command:
scripts/mainScript.sh "<full path to the salmon binary> $( pwd -LP . ) analysisDir

If everything was installed correctly you should see the progress message from salmon and eventually the number of differential expressed genes. The analysisDir will contain volcano and heatmap plots and excel files with the top differential expressed genes.

Unfortunatelly, I won't be able to answer any questions about the software, since I am not working at UA .  Please contact Nancy Casanova at ncasanova@arizona.edu if you need any help.


Good luck!


Manuel L. Gonzalez-Garay Ph.D.
