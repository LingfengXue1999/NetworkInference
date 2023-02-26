# Codes to investigate network inference in real datasets

1. Download datasets, reads mapping using CellRanger, calculate spliced/unspliced matrix using Velocyto
2. Read spliced/unspliced matrix, QC by Seurat (read_velocyto_output_hBCell.R)
<br> The input spliced/unspliced matrices have been submitted to figshare (https://figshare.com/articles/dataset/Datasets/21342117)
4. Network inference by GENIE3 (genie3_for_exon_data_hBCell.R, genie3_for_inex_data_hBCell.R)
5. Evaluate network inference resuls (Evaluate_network_inference_scRNAseq.R)
6. Summary for evalulation and plot (Summary_evaluation_30_samples.R)
7. Factor analysis for results (art_factor_dependency_analysis.R, analysis_var_of_TFs.R)


