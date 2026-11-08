# Endothelial cells as sources of immune checkpoint ligands and candidate regulators of PD-1 signalling in non-small cell lung cancer 

## Repository layout

```
config/paths.R            input/output locations (set with environment variables)
R/utils.R                 shared helpers (Seurat workflow, QC, inferCNV, marker plots)
analysis/
  01_scrna_processing/    QC and initial processing of public and internal sc/snRNA-seq datasets
  02_endothelial_cells/   EC identification, filtering, batch correction, and subtyping
  03_t_cells/             T-cell lineages, CD4 subtypes, CD8 trajectories
  04_icl_expression/      Statistical comparisons of ICL expression
  05_cellchat/            Scripts for running CellChat and generating figures
  06_spatial/             CosMx, Xenium and CODEX processing and proximity statistics
figures.R    code for every computational figure panel
```
## Data

| Dataset | Source |
|---|---|
| Internal sc/snRNA-seq (6 fresh-frozen, 28 FFPE LUSC) | See manuscript for data accession numbers |
| NSCLC and normal lung scRNA-seq atlas | Salcher et al., 2022, CZ CELLxGENE |
| Cell-type annotation reference | Guimarães et al., 2024, CZ CELLxGENE |
| Human Primary Cell Atlas (SingleR) | `celldex::HumanPrimaryCellAtlasData()` |
| CosMx SMI NSCLC FFPE dataset | Bruker Spatial Biology |
| CODEX mIF, 35 NSCLC patients | Zheng et al., 2025, Dryad |
| Xenium Prime 5K human LUAD | 10x Genomics |

Data files are not stored in this repository.

## Software
```
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] speckle_0.0.3                           
 [2] edgeR_3.40.2                            
 [3] limma_3.54.2                            
 [4] emmeans_1.11.2                          
 [5] fs_1.6.2                                
 [6] cytomapper_1.10.1                       
 [7] EBImage_4.40.1                          
 [8] imcRtools_1.4.2                         
 [9] SpatialExperiment_1.8.1                 
[10] patchwork_1.3.2.9000                    
[11] ggforce_0.4.2                           
[12] RANN_2.6.1                              
[13] sfheaders_0.4.4                         
[14] sf_1.0-21                               
[15] STopover_1.1.2                          
[16] zellkonverter_1.8.0                     
[17] caret_6.0-94                            
[18] lattice_0.20-45                         
[19] qs2_0.1.7                               
[20] httpgd_2.0.4                            
[21] glmmTMB_1.1.11                          
[22] slingshot_2.6.0                         
[23] TrajectoryUtils_1.6.0                   
[24] princurve_2.1.6                         
[25] viridis_0.6.3                           
[26] viridisLite_0.4.2                       
[27] CellChat_2.1.2                          
[28] bigmemory_4.6.4                         
[29] igraph_2.0.3                            
[30] SeuratExtend_1.1.0                      
[31] SeuratExtendData_0.2.1                  
[32] SeuratWrappers_0.3.2                    
[33] ggrepel_0.9.3                           
[34] RColorBrewer_1.1-3                      
[35] ggsci_3.0.0                             
[36] ggpubr_0.6.0                            
[37] harmony_1.2.4                           
[38] Rcpp_1.0.11                             
[39] infercnv_1.3.3                          
[40] scDblFinder_1.12.0                      
[41] BSgenome.Hsapiens.UCSC.hg38_1.4.5       
[42] BSgenome_1.66.3                         
[43] rtracklayer_1.58.0                      
[44] Biostrings_2.66.0                       
[45] XVector_0.38.0                          
[46] TFBSTools_1.36.0                        
[47] JASPAR2020_0.99.10                      
[48] motifmatchr_1.20.0                      
[49] biomaRt_2.54.1                          
[50] TxDb.Hsapiens.UCSC.hg38.knownGene_3.16.0
[51] pheatmap_1.0.12                         
[52] EnsDb.Hsapiens.v86_2.99.0               
[53] ensembldb_2.22.0                        
[54] AnnotationFilter_1.22.0                 
[55] GenomicFeatures_1.50.4                  
[56] AnnotationDbi_1.60.2                    
[57] scran_1.26.2                            
[58] scRNAseq_2.12.0                         
[59] scater_1.26.1                           
[60] scuttle_1.8.4                           
[61] SingleCellExperiment_1.20.1             
[62] SummarizedExperiment_1.28.0             
[63] Biobase_2.58.0                          
[64] GenomicRanges_1.50.2                    
[65] GenomeInfoDb_1.34.9                     
[66] IRanges_2.32.0                          
[67] S4Vectors_0.36.2                        
[68] BiocGenerics_0.44.0                     
[69] MatrixGenerics_1.10.0                   
[70] matrixStats_1.0.0                       
[71] Signac_1.14.0                           
[72] Seurat_5.2.1                            
[73] SeuratObject_5.0.2                      
[74] sp_2.0-0                                
[75] Matrix_1.6-4                            
[76] data.table_1.14.8                       
[77] lubridate_1.9.2                         
[78] forcats_1.0.0                           
[79] stringr_1.5.0                           
[80] dplyr_1.1.2                             
[81] purrr_1.0.1                             
[82] readr_2.1.4                             
[83] tidyr_1.3.0                             
[84] tibble_3.2.1                            
[85] ggplot2_4.0.2                           
[86] tidyverse_2.0.0                         

loaded via a namespace (and not attached):
  [1] KEGGREST_1.38.0               spam_2.9-1                   
  [3] circlize_0.4.15               coin_1.4-2                   
  [5] locfit_1.5-9.8                remotes_2.4.2                
  [7] BiocVersion_3.16.0            distances_0.1.13             
  [9] spatstat.utils_3.1-3          vctrs_0.6.2                  
 [11] utf8_1.2.3                    blob_1.2.4                   
 [13] svgPanZoom_0.3.4              R.oo_1.25.0                  
 [15] withr_2.5.0                   ggnetwork_0.5.13             
 [17] S7_0.2.0                      registry_0.5-1               
 [19] uuid_1.1-0                    lifecycle_1.0.3              
 [21] ScaledMatrix_1.6.0            ggalluvial_0.12.5            
 [23] codetools_0.2-18              lmtest_0.9-40                
 [25] org.Hs.eg.db_3.16.0           annotate_1.76.0              
 [27] magick_2.9.0                  parallelly_1.36.0            
 [29] AnnotationHub_3.6.0           fastmatch_1.1-3              
 [31] basilisk_1.10.2               metapod_1.6.0                
 [33] Rtsne_0.16                    stringi_1.7.12               
 [35] RcppRoll_0.3.0                sctransform_0.4.1            
 [37] grid_4.2.2                    polyclip_1.10-4              
 [39] rhdf5filters_1.10.1           sandwich_3.0-2               
 [41] goftest_1.2-3                 cluster_2.1.4                
 [43] ggraph_2.2.1                  ape_5.7-1                    
 [45] TFMPvalue_0.0.9               pkgconfig_2.0.3              
 [47] prettyunits_1.1.1             sparseMatrixStats_1.10.0     
 [49] ggridges_0.5.6                timechange_0.2.0             
 [51] estimability_1.5.1            nnls_1.6                     
 [53] org.Mm.eg.db_3.16.0           httr_1.4.6                   
 [55] fftwtools_0.9-11              DropletUtils_1.18.1          
 [57] progress_1.2.2                fastcluster_1.2.3            
 [59] GetoptLong_1.0.5              modeltools_0.2-23            
 [61] terra_1.7-39                  beachmat_2.14.2              
 [63] graphlayouts_1.1.0            basilisk.utils_1.10.0        
 [65] htmltools_0.5.5               miniUI_0.1.1.1               
 [67] NMF_0.26                      yaml_2.3.7                   
 [69] prodlim_2023.03.31            pillar_1.9.0                 
 [71] jquerylib_0.1.4               later_1.3.1                  
 [73] fitdistrplus_1.1-11           glue_1.6.2                   
 [75] DBI_1.1.3                     BiocParallel_1.32.6          
 [77] plyr_1.8.8                    foreach_1.5.2                
 [79] ProtGenerics_1.30.0           dotCall64_1.0-2              
 [81] ExperimentHub_2.6.0           gtable_0.3.6                 
 [83] rsvd_1.0.5                    raster_3.6-26                
 [85] caTools_1.18.2                GlobalOptions_0.1.2          
 [87] fastmap_1.1.1                 broom_1.0.5                  
 [89] promises_1.2.0.1              FNN_1.1.3.2                  
 [91] hms_1.1.3                     png_0.1-8                    
 [93] clue_0.3-64                   spatstat.explore_3.1-0       
 [95] lazyeval_0.2.2                crayon_1.5.2                 
 [97] gridBase_0.4-7                svglite_2.1.1                
 [99] RTriangle_1.6-0.15            boot_1.3-28                  
[101] tidyselect_1.2.0              BiocSingular_1.14.0          
[103] unigd_0.1.3                   splines_4.2.2                
[105] survival_3.4-0                rappdirs_0.3.3               
[107] xgboost_1.7.7.1               bit64_4.0.5                  
[109] rngtools_1.5.2                CNEr_1.34.0                  
[111] lambda.r_1.2.4                jpeg_0.1-10                  
[113] futile.logger_1.4.3           ggsignif_0.6.4               
[115] R.methodsS3_1.8.2             stringfish_0.18.0            
[117] xtable_1.8-4                  DT_0.28                      
[119] cachem_1.0.8                  DelayedArray_0.24.0          
[121] ipred_0.9-14                  abind_1.4-5                  
[123] systemfonts_1.0.4             mime_0.12                    
[125] rjson_0.2.21                  rstatix_0.7.2                
[127] processx_3.8.2                spatstat.sparse_3.0-1        
[129] numDeriv_2016.8-1.1           tools_4.2.2                  
[131] cli_3.6.1                     Rdpack_2.6                   
[133] magrittr_2.0.3                proxy_0.4-27                 
[135] dichromat_2.0-0.1             future.apply_1.11.0          
[137] bigmemory.sri_0.1.8           libcoin_1.0-9                
[139] DelayedMatrixStats_1.20.0     ggbeeswarm_0.7.2             
[141] sna_2.7-2                     HDF5Array_1.26.0             
[143] ica_1.0-3                     pbapply_1.7-2                
[145] argparse_2.2.2                BiocFileCache_2.6.1          
[147] mgcv_1.8-41                   plotly_4.10.2                
[149] R.utils_2.12.2                tweenr_2.0.2                 
[151] multcomp_1.4-25               zlibbioc_1.44.0              
[153] restfulr_0.0.15               tzdb_0.4.0                   
[155] ps_1.7.5                      fansi_1.0.4                  
[157] tidygraph_1.3.1               TH.data_1.1-2                
[159] tensor_1.5                    ROCR_1.0-11                  
[161] KernSmooth_2.23-20            backports_1.4.1              
[163] farver_2.1.1                  bit_4.0.5                    
[165] Rsamtools_2.14.0              gplots_3.1.3                 
[167] rjags_4-14                    shiny_1.7.4.1                
[169] BiocIO_1.8.0                  sess_3.2.1                   
[171] scattermore_1.2               rbibutils_2.2.16             
[173] hardhat_1.3.0                 sass_0.4.7                   
[175] RcppAnnoy_0.0.21              futile.options_1.0.1         
[177] pROC_1.18.4                   rstudioapi_0.15.0            
[179] iterators_1.0.14              minqa_1.2.5                  
[181] Rhdf5lib_1.20.0               spatstat.geom_3.1-0          
[183] RcppHNSW_0.4.1                nlme_3.1-160                 
[185] DirichletMultinomial_1.40.0   shape_1.4.6                  
[187] gtools_3.9.5                  beeswarm_0.4.0               
[189] network_1.18.2                bslib_0.5.0                  
[191] TMB_1.9.17                    listenv_0.9.0                
[193] reshape2_1.4.4                rhdf5_2.42.1                 
[195] generics_0.1.3                colorspace_2.1-0             
[197] XML_3.99-0.13                 ModelMetrics_1.2.2.2         
[199] e1071_1.7-13                  spatstat.data_3.0-1          
[201] dbplyr_2.3.3                  dqrng_0.3.0                  
[203] GenomeInfoDbData_1.2.9        timeDate_4022.108            
[205] progressr_0.18.0              memoise_2.0.1                
[207] coda_0.19-4                   ComplexHeatmap_2.15.1        
[209] doParallel_1.0.17             vipor_0.4.5                  
[211] httpuv_1.6.11                 class_7.3-20                 
[213] irlba_2.3.5.1                 BiocManager_1.30.27          
[215] formatR_1.14                  classInt_0.4-9               
[217] seqLogo_1.64.0                vroom_1.6.3                  
[219] jsonlite_1.8.7                RSpectra_0.16-1              
[221] dir.expiry_1.6.0              digest_0.6.33                
[223] poweRlaw_0.70.6               cowplot_1.2.0.9000           
[225] bitops_1.0-7                  RSQLite_2.3.1                
[227] globals_0.16.2                compiler_4.2.2               
[229] nnet_7.3-18                   reticulate_1.44.1            
[231] statmod_1.5.0                 zoo_1.8-12                   
[233] carData_3.0-5                 pracma_2.4.2                 
[235] interactiveDisplayBase_1.36.0 rlang_1.1.0                  
[237] nloptr_2.0.3                  uwot_0.1.16                  
[239] lava_1.7.2.1                  recipes_1.0.6                
[241] future_1.33.0                 mvtnorm_1.2-2                
[243] htmlwidgets_1.6.2             curl_5.0.1                   
[245] concaveman_1.2.0              parallel_4.2.2               
[247] fastDummies_1.7.3             BiocNeighbors_1.16.0         
[249] filelock_1.0.2                scales_1.4.0                 
[251] RcppParallel_5.1.7            lme4_1.1-34                  
[253] deldir_1.0-9                  gridExtra_2.3                
[255] tiff_0.1-12                   bluster_1.8.0                
[257] RCurl_1.98-1.12               car_3.1-2                    
[259] GO.db_3.16.0                  MASS_7.3-58.1                
[261] ellipsis_0.3.2                spatstat.random_3.1-4        
[263] xml2_1.3.5                    shinydashboard_0.7.3         
[265] gower_1.0.1                   reshape_0.8.9                
[267] reformulas_0.4.0              rpart_4.1.19                 
[269] R6_2.5.1                      statnet.common_4.9.0         
[271] GenomicAlignments_1.34.1      units_0.8-3      
```