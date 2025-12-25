# This script calculates FDR and processes Phase 1 results.

# 获取命令行参数 (每次 run 时传入变量)
args <- commandArgs(trailingOnly = TRUE)
merged_file  <- args[1]
sig_out_file <- args[2]
exp_name     <- args[3]
#exp_name     <- 'AA'
exp_path <- paste0('exposure/', exp_name,'/')
#merged_file  <-'CD4_Dynamic_demo_NECROSIS_PERIAPICAL/Phase1_All_Results.txt'


suppressPackageStartupMessages({
library(qs)
library(dplyr)
library(data.table)
})

# A. Load Phase 1 Results
tryCatch({
  res <- read.csv(merged_file,sep='\t')
}, error = function(e) { stop('Error reading merged file') })

if (nrow(res) > 0) {
  # B. Calculate FDR
  target_methods <- c('Inverse variance weighted', 'Wald ratio', 'IVW')
  ivw_res <- subset(res, method %in% target_methods)
  
  if (nrow(ivw_res) > 0) {
    ivw_res$fdr <- p.adjust(ivw_res$p, method = 'fdr')
    
    # C. Identify Significant Genes (FDR < 0.05)
    sig_genes <- unique(ivw_res$exposure[ivw_res$fdr < 0.05])
    
    cat(paste0('Found ', length(sig_genes), ' significant genes (FDR < 0.05).\n'))
    
    if (length(sig_genes) > 0) {
      # 保存显著基因列表
      write.table(sig_genes, sig_out_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      # D. Prepare Phase 2 Data (R2=0.01)
      r2_01_file <- paste0(exp_path, exp_name,'_r2_0.01.qs')
      
      if (file.exists(r2_01_file)) {
        cat('Loading R2=0.01 data and subsetting...\n')
        exp_01 <- qread(r2_01_file)
        
        exp_01_sub <- subset(exp_01, exposure %in% sig_genes)
        
        # 保存为新的临时文件
        out_file <- paste0(exp_path, exp_name, '_r2_0.01_phase2.qs')
        qsave(exp_01_sub, out_file)
        cat('Phase 2 input data prepared:', out_file, '\n')
      } else {
        cat('Error: R2=0.01 source file not found.\n')
      }
    } else {
      cat('No significant genes found. Phase 2 will be skipped.\n')
    }
  } else {
    cat('No IVW/Wald ratio results found.\n')
  }
} else {
  cat('Merged file is empty.\n')
}
