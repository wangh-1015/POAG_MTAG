##### TGFM
library(data.table)
library(dplyr)
library(tidyr)

### 0. QTL data preprocess
# upload
setwd("/home/wangh/DATA/GTEx/GTEx_Analysis_v8_eQTL/sig_s2g_pair")
MTAG <- fread("/home/wangh/project1_glaucoma/analysis/TGFM/MTAG_hg38.txt", sep = "\t", header = TRUE)
MTAG <- MTAG[,c(2,1,3)]
file_list <- list.files(pattern = "\\.v8\\.signif_variant_gene_pairs\\.txt$")
sample_size <- fread("/home/wangh/DATA/GTEx/GTEx_Analysis_v8_eQTL/sample_size.txt", sep = "\t", header = TRUE)

# munge
lapply(file_list, function(file_path) {
  
  # 
  file_name <- sub("\\.v8\\.signif_variant_gene_pairs\\.txt$", "", basename(file_path))
  file_name <- tolower(file_name)
  print(paste("处理文件:", file_name))
  
  # 
  df <- fread(file_path)
  
  # add N
  if (file_name %in% sample_size$V1) {
    df$N <- sample_size$V2[sample_size$V1 == file_name]
  } else {
    df$N <- 141
  }
  print("add N!")
  
  #
  library(tidyr)
  library(dplyr)
  
  df <- df %>% 
    separate(variant_id, into = c("CHR", "SNP_BP", "A1", "A2", "build"), sep = "_") %>% 
    mutate(SNP_BP = as.numeric(SNP_BP), GENE_COORD = SNP_BP - tss_distance) %>% 
    setnames(old = c("gene_id", "slope", "slope_se"), new = c("GENE", "BETA", "BETA_VAR")) %>% 
    select(GENE, CHR, GENE_COORD, SNP_BP, A1, A2, N, BETA, BETA_VAR)
  print("pro1!")
  
  df$CHR <- gsub("chr", "", df$CHR)
  df <- df %>%
    mutate(CHR = as.numeric(CHR)) %>% 
    left_join(MTAG, by = c("CHR" = "CHR", "SNP_BP" = "BP")) %>% 
    filter(!is.na(SNP)) %>% 
    select(GENE, SNP, CHR, GENE_COORD, SNP_BP, A1, A2, N, BETA, BETA_VAR)
  print("pro2!")
  
  # 
  output_file <- paste0("/home/wangh/DATA/GTEx/GTEx_Analysis_v8_eQTL/eQTL_summary/", file_name, ".txt")
  write.table(df, 
             file = output_file, 
             sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  #
  rm(df)
  
  return(NULL)
})


##########################################################
###！                      TGFM！                      ###
##########################################################
### 1. Estimating causal eQTL effect size distributions### 
nohup bash -c 'source /home/wangh/miniconda3/bin/activate && \
cd /home/wangh/project1_glaucoma/analysis2/TGFM/qtl_dist_small

conda activate TGFM
export R_LIBS_USER=/home/wangh/miniconda3/envs/tgfm_r/lib/R/library
export LD_LIBRARY_PATH=/home/wangh/miniconda3/envs/tgfm_r/lib:$LD_LIBRARY_PATH

for chrom in {1..22}
do
  tissue_out_dir="/home/wangh/project1_glaucoma/analysis2/TGFM/qtl_dist_small/$chrom/qtl"
  dir="/home/wangh/project1_glaucoma/analysis2/TGFM/qtl_dist_small/$chrom"
  mkdir -p "$dir"
  
    python /home/wangh/software/TGFM/susie_eqtl_fine_mapping_for_tgfm.py \
        --eqtl-data-type SumStat \
        --chrom $chrom \
        --genotype-stem /home/wangh/DATA/ref/LDREF_hg38/small_ref/1000G.EUR. \
        --eqtl-sumstat /home/wangh/project1_glaucoma/analysis2/TGFM/retina_sig_eqtl.txt \
        --gwas-sumstat /home/wangh/project1_glaucoma/analysis/TGFM/MTAG_hg38.txt \
        --filter-strand-ambiguous \
        --min-cis-snps-per-gene 5 \
        --out "$tissue_out_dir"
done

' > /home/wangh/project1_glaucoma/analysis2/TGFM/qtl_dist_small.log 2>&1 &

 # tissue chr double loop
nohup bash -c 'source /home/wangh/miniconda3/bin/activate && \

conda activate TGFM && \
export R_HOME=${CONDA_PREFIX}/lib/R
export R_ENVIRON=${CONDA_PREFIX}/lib/R/etc/Renviron
export R_PROFILE=${CONDA_PREFIX}/lib/R/etc/Rprofile.site

QTL_DIR="/home/wangh/DATA/GTEx/GTEx_Analysis_v8_eQTL/eQTL_summary/TGFM"
RESULT_BASE_DIR="/home/wangh/project1_glaucoma/analysis/TGFM/result/QTL_dis1"

tissue_files=$(ls ${QTL_DIR}/*.txt)

for tissue_file in $tissue_files; do
  tissue_name=$(basename "$tissue_file" .txt)
  echo "Processing tissue: $tissue_name"
  
  tissue_out_dir="${RESULT_BASE_DIR}/${tissue_name}"
  mkdir -p "$tissue_out_dir"
  
  for chrom in {1..22}; do
    echo "  Processing chromosome: $chrom"
    
    chrom_out_dir="${tissue_out_dir}/${chrom}"
    mkdir -p "$chrom_out_dir"
    
    python /home/wangh/software/TGFM/susie_eqtl_fine_mapping_for_tgfm.py \
      --eqtl-data-type SumStat \
      --chrom $chrom \
      --genotype-stem /home/wangh/DATA/ref/LDREF_hg38/small_ref/1000G.EUR. \
      --eqtl-sumstat "$tissue_file" \
      --gwas-sumstat /home/wangh/project1_glaucoma/analysis/TGFM/MTAG_hg38.txt \
      --filter-strand-ambiguous \
      --out "$chrom_out_dir/"
  done
done
 
' > /home/wangh/project1_glaucoma/analysis/TGFM/result/QTL_dis1/scDRS.log 2>&1 &

library(dplyr)
library(readr)

input_dir <- "/home/wangh/project1_glaucoma/analysis/TGFM/result/QTL_dis1"
output_dir <- file.path(input_dir, "A_susie_eqtl_output_stem")
dir.create(output_dir)

tissue_folders <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)

for (tissue_path in tissue_folders) {
  tissue_name <- basename(tissue_path)
  
  cat("processing：", tissue_name, "\n")
  
  combined_list <- list()
  
  for (chr_num in 1:22) {
    chr_path <- file.path(tissue_path, as.character(chr_num))
    
    summary_file <- list.files(chr_path, pattern = "_chr.*_gene_summary.txt$", full.names = TRUE)
    
    if (length(summary_file) == 1) {
      
      df <- tryCatch({
  tmp <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  tmp[] <- lapply(tmp, as.character)
  tmp
}, error = function(e) {
  cat("fail：", summary_file, "\n")
  return(NULL)
})
      tmp$CHR <- as.numeric(tmp$CHR)
      tmp$GENE_COORD <- as.numeric(tmp$GENE_COORD)
      
      if (!is.null(df)) {
        combined_list[[length(combined_list) + 1]] <- df
      }
    }
  }
  
  if (length(combined_list) > 0) {
    combined_df <- bind_rows(combined_list)
    
    output_file <- file.path(output_dir, paste0(tissue_name, ".txt"))
    write.table(combined_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    cat("跳过组织", tissue_name, "no_data \n")
  }
}

cat("所有组织文件处理完毕。\n")

### 
### 2. 

## （1） GWAS beta se

MTAG_filter$BETA <- approx(MTAG_filter$BP, MTAG_filter$BETA, xout = MTAG_filter$BP, method = "linear")$y


MTAG_filter$BETA_VAR <- approx(MTAG_filter$BP, MTAG_filter$BETA_VAR, xout = MTAG_filter$BP, method = "linear")$y


## （2） window_df
generate_windows <- function(bim_data) {
  windows <- data.frame(
    window_name = character(),
    chr = numeric(),
    window_start = numeric(),
    window_end = numeric(),
    stringsAsFactors = FALSE
  )
  
  chr_list <- unique(bim_data$CHR)
  
  for (chr_num in chr_list) {
    chr_data <- bim_data[bim_data$CHR == chr_num,]
    
    min_pos <- min(chr_data$BP)
    max_pos <- max(chr_data$BP)
    
    start_pos <- max(1, min_pos - 1000)  # 使用1kb (1000bp)而不是1MB
    
    current_start <- start_pos
    while (current_start <= max_pos) {
      window_start <- current_start
      window_end <- current_start + 3000000  # 3MB窗口
      
      # 创建窗口名称，格式为 chr:start:end
      window_name <- paste(chr_num, window_start, window_end, sep = ":")
      
      # 添加窗口到我们的数据框
      window_row <- data.frame(
        window_name = window_name,
        chr = chr_num,
        window_start = window_start,
        window_end = window_end,
        stringsAsFactors = FALSE
      )
      windows <- rbind(windows, window_row)
      
      # 移动到下一个窗口，步长为1MB
      current_start <- current_start + 1000000
    }
    
    # 确保最后一个窗口延伸到最后一个变异位点之后1MB
    if (max(windows$window_end[windows$chr == chr_num]) < max_pos + 1000000) {
      window_start <- max_pos - 2000000  # 确保是3MB窗口
      window_end <- max_pos + 1000000
      window_name <- paste(chr_num, window_start, window_end, sep = ":")
      
      window_row <- data.frame(
        window_name = window_name,
        chr = chr_num,
        window_start = window_start,
        window_end = window_end,
        stringsAsFactors = FALSE
      )
      windows <- rbind(windows, window_row)
    }
  }
  
  return(windows)
}

# 生成窗口
window_df <- generate_windows(bim)

## (3) 生成每个窗口的bim文件 共2752个 variant_info_file

for (i in 1:nrow(window_df)) {
  # Get the current window information
  window_start <- window_df$window_start[i]
  window_end <- window_df$window_end[i]
  chr <- window_df$chr[i]
  window_name <- window_df$window_name[i]
  
  # Filter the bim data to find variants within the current window
  variants_in_window <- bim[bim$CHR == chr & bim$BP >= window_start & bim$BP <= window_end, ]
  
  # If there are variants in this window, create a file
  if (nrow(variants_in_window) > 0) {
    # Create file name based on the format specified
    file_name <- paste0("/home/wangh/project1_glaucoma/analysis/TGFM/variant_info_file/", window_name, "_variant_info.txt")
    
    # Write the filtered data to the file (tab-delimited, same format as PLINK bim file)
    write.table(variants_in_window[, c("CHR", "SNP", "X", "BP", "A1", "A2")], 
                file = file_name, 
                row.names = FALSE, 
                col.names = FALSE, 
                sep = "\t", 
                quote = FALSE)
    
    cat("File created: ", file_name, "\n")
  } else {
    cat("No variants found in window: ", window_name, "\n")
  }
}

## (4)为每个窗口生成LD_matrix
# 1) 每个窗口生成SNP列表
input_dir <- "/home/wangh/project1_glaucoma/analysis/TGFM/variant_info_file"
output_dir <- "/home/wangh/project1_glaucoma/analysis/TGFM/LD_file/snp"

files <- list.files(input_dir, pattern = "*.txt", full.names = TRUE)

for (file in files) {
  # 读取文件，假设文件是无列名的制表符分隔文件
  data <- read.table(file, header = FALSE, sep = "\t")

  snp_column <- data[, 2]
  
  output_file <- file.path(output_dir, basename(file))
  write.table(snp_column, file = output_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

  cat("Processed and saved SNP column from:", file, "\n")
}

# 2) 根据每个snp列表提取hg38的bed文件
# genotype extract (for loop)
conda activate plink
cd /home/wangh/DATA/ref/LDREF_hg38/G1000_all
plink --bfile 1000G.EUR_all \
      --extract /home/wangh/project1_glaucoma/analysis/TGFM/LD_file/snp/1:816186:3816186_variant_info.txt \
      --make-bed \
      --out /home/wangh/project1_glaucoma/analysis/TGFM/LD_file/bed/1:816186:3816186_variant_info

# 3) LD calculation (for loop)
nohup bash -c 'source /home/wangh/miniconda3/bin/activate && \

conda activate plink
cd /home/wangh/project1_glaucoma/analysis/TGFM/LD_file/bed

input_dir="/home/wangh/project1_glaucoma/analysis/TGFM/LD_file/bed"
output_dir="/home/wangh/project1_glaucoma/analysis/TGFM/LD_file/ld"

# 遍历输入目录中的所有.bed文件（因为.bed文件是必需的）
for bedfile in $input_dir/*_variant_info.txt.bed; do
    # 提取文件名的前缀部分，去掉路径和扩展名
    filename_in=$(basename "$bedfile" ".bed")
    filename_out=$(basename "$bedfile" "_variant_info.txt.bed")

    # 使用 plink 计算 LD 矩阵并输出到指定目录
    plink --bfile "$input_dir/$filename_in" \
          --r2 \
          --matrix \
          --out "$output_dir/$filename_out"

    # 输出处理信息
    echo "LD matrix computed for: $filename_out"
done

' > /home/wangh/project1_glaucoma/analysis/TGFM/LD_file/ld.log 2>&1 &

## (5) 制作输入文件
conda activate TGFM
export R_LIBS_USER=/home/wangh/miniconda3/envs/tgfm_r/lib/R/library
export LD_LIBRARY_PATH=/home/wangh/miniconda3/envs/tgfm_r/lib:$LD_LIBRARY_PATH

for chr in {1..21}; do
tissue_summary="/home/wangh/project1_glaucoma/analysis2/TGFM/tissue_summary_file/tissue_summary_file$chr"
output_file="/home/wangh/project1_glaucoma/analysis2/TGFM/result/input2/$chr/qtl"
output_dir="/home/wangh/project1_glaucoma/analysis2/TGFM/result/input2/$chr"
mkdir -p $output_dir

python /home/wangh/software/TGFM/prepare_input_data_for_tgfm.py \
    --chrom $chr \
    --window-file /home/wangh/project1_glaucoma/analysis/TGFM/window_df \
    --tissue-summary-file $tissue_summary \
    --gwas-summary-file /home/wangh/project1_glaucoma/analysis/TGFM/gwas_summary_file \
    --standardize-gwas-summary-statistics \
    --out $output_file
done
    
###  
nohup bash -c 'source /home/wangh/miniconda3/bin/activate && \

conda activate TGFM
export R_HOME=${CONDA_PREFIX}/lib/R
export R_ENVIRON=${CONDA_PREFIX}/lib/R/etc/Renviron
export R_PROFILE=${CONDA_PREFIX}/lib/R/etc/Rprofile.site
cd /home/wangh/project1_glaucoma/analysis/TGFM

window_file="window_df"
gwas_summary_file="gwas_summary_file"
base_output_dir="/home/wangh/project1_glaucoma/analysis/TGFM/result/input"

# 遍历染色体 1 到 22,每个染色体同时运行49个组织的summary data
for chr in {1..22}; do
    # 生成对应染色体的 tissue-summary-file 文件路径
    tissue_summary_file="/home/wangh/project1_glaucoma/analysis/TGFM/tissue_summary1_22/tissue_summary_file${chr}"
    
    # 创建输出目录（确保目录存在），使用chr作为文件夹名
    output_dir="${base_output_dir}/chr${chr}"
    mkdir -p $output_dir  # 确保该目录存在
    
    # 运行 TGFM 脚本，传递染色体号和对应的 tissue-summary-file
    python /home/wangh/software/TGFM/prepare_input_data_for_tgfm.py \
    --chrom $chr \
    --window-file $window_file \
    --tissue-summary-file $tissue_summary_file \
    --gwas-summary-file $gwas_summary_file \
    --standardize-gwas-summary-statistics \
    --out "${output_dir}/"  # 这里确保文件输出到特定的子目录
    
    # 打印执行的命令（可选，用于调试）
    echo "Executed TGFM for chromosome $chr, output saved in $output_dir"
done

' > /home/wangh/project1_glaucoma/analysis/TGFM/result/input/tgfm_input.log 2>&1 &


### 3. run no sampling TGFM (文档建议分10个chunk，但是我不理解其意义，依然按染色体分析)
### 注意：这里的export R有别于上面两步，否则无法运行（虽然我也不知道为啥）
nohup bash -c 'source /home/wangh/miniconda3/bin/activate && \

conda activate TGFM
export R_HOME=${CONDA_PREFIX}/lib/R
export R_ENVIRON=${CONDA_PREFIX}/lib/R/etc/Renviron
export R_PROFILE=${CONDA_PREFIX}/lib/R/etc/Rprofile.site

for chr in {1..21}
do
input="/home/wangh/project1_glaucoma/analysis2/TGFM/result/input/${chr}/qtl_chr${chr}_input_data_summary.txt"
out="/home/wangh/project1_glaucoma/analysis2/TGFM/result/no_sample/${chr}/qtl"
out_dir="/home/wangh/project1_glaucoma/analysis2/TGFM/result/no_sample/${chr}"
mkdir -p $out_dir

python /home/wangh/software/TGFM/run_tgfm_without_sampling.py \
--trait-name glaucoma \
--tgfm-input-data $input \
--parallel-job-identifier 1 \
--gene-tissue-pip-threshold 0 \
--p-value-threshold 1 \
--out $out

done
' > /home/wangh/project1_glaucoma/analysis2/TGFM/result/no_sample/no_sample.log 2>&1 &

## bash loop
nohup bash -c 'source /home/wangh/miniconda3/bin/activate && \

conda activate TGFM
export R_HOME=${CONDA_PREFIX}/lib/R
export R_ENVIRON=${CONDA_PREFIX}/lib/R/etc/Renviron
export R_PROFILE=${CONDA_PREFIX}/lib/R/etc/Rprofile.site

input_base_dir="/home/wangh/project1_glaucoma/analysis/TGFM/result/input"
output_base_dir="/home/wangh/project1_glaucoma/analysis/TGFM/result/no_sample"

trait_name="glaucoma"

# 循环处理每个chr文件夹
for i in {1..22}
do

    chr_folder=$(printf "chr%d" $i)

    input_file="${input_base_dir}/${chr_folder}/_${chr_folder}_input_data_summary.txt"

    output_dir="${output_base_dir}/${chr_folder}"
    out_stem="${output_base_dir}/${chr_folder}/TGFM"
    mkdir -p "$output_dir"

    echo "Running TGFM for ${chr_folder}..."
    
    python /home/wangh/software/TGFM/run_tgfm_without_sampling.py \
    --trait-name $trait_name \
    --tgfm-input-data $input_file \
    --parallel-job-identifier 1 \
    --gene-tissue-pip-threshold 0 \
    --p-value-threshold 1 \
    --out $out_stem
done

' > /home/wangh/project1_glaucoma/analysis/TGFM/result/no_sample/tgfm_no_sample.log 2>&1 &

### 4. TGFM tissue prior

for chr in {1..21}
do

input_summary="/home/wangh/project1_glaucoma/analysis2/TGFM/tissue_summary_file/tissue_summary_file${chr}"
out="/home/wangh/project1_glaucoma/analysis2/TGFM/result/prior/qtl_${chr}"
input_nosample="/home/wangh/project1_glaucoma/analysis2/TGFM/result/no_sample/${chr}/qtl"

python /home/wangh/software/TGFM/run_tgfm_tissue_specific_prior.py \
--trait-name glaucoma \
--tissue-summary-file $input_summary \
--tgfm-parallel-job-identifier-file /home/wangh/project1_glaucoma/analysis/TGFM/result/no_sample/identify \
--tgfm-without-sampling-output $input_nosample \
--out $out

done
##
conda activate TGFM
export R_HOME=${CONDA_PREFIX}/lib/R
export R_ENVIRON=${CONDA_PREFIX}/lib/R/etc/Renviron
export R_PROFILE=${CONDA_PREFIX}/lib/R/etc/Rprofile.site

input_tissue_base="/home/wangh/project1_glaucoma/analysis/TGFM/tissue_summary1_22"
no_sample_base="/home/wangh/project1_glaucoma/analysis/TGFM/result/no_sample"
out_base="/home/wangh/project1_glaucoma/analysis/TGFM/result/tissue_prior"

for i in {1..22}
do

    chr_folder=$(printf "chr%d" $i)
    chr=$i
    
    input_tissue_file="${input_tissue_base}/tissue_summary_file$chr"
    input_nosample="${no_sample_base}/${chr_folder}/TGFM"

    out_stem="${out_base}/${chr_folder}"
    
python /home/wangh/software/TGFM/run_tgfm_tissue_specific_prior.py \
--trait-name glaucoma \
--tissue-summary-file $input_tissue_file \
--tgfm-parallel-job-identifier-file /home/wangh/project1_glaucoma/analysis/TGFM/result/no_sample/identify \
--tgfm-without-sampling-output $input_nosample \
--out $out_stem

done


### 5. TGFM final run
nohup bash -c 'source /home/wangh/miniconda3/bin/activate && \
conda activate TGFM
export R_HOME=${CONDA_PREFIX}/lib/R
export R_ENVIRON=${CONDA_PREFIX}/lib/R/etc/Renviron
export R_PROFILE=${CONDA_PREFIX}/lib/R/etc/Rprofile.site

for chr in {1..21}
do

input_summary="/home/wangh/project1_glaucoma/analysis2/TGFM/tissue_summary_file/tissue_summary_file${chr}"
input="/home/wangh/project1_glaucoma/analysis2/TGFM/result/input/${chr}/qtl_chr${chr}_input_data_summary.txt"
prior="/home/wangh/project1_glaucoma/analysis2/TGFM/result/prior/qtl_${chr}"
out="/home/wangh/project1_glaucoma/analysis2/TGFM/result/sample/${chr}/qtl"
outdir="/home/wangh/project1_glaucoma/analysis2/TGFM/result/sample/${chr}"
mkdir -p "$outdir"

python /home/wangh/software/TGFM/run_tgfm.py \
--trait-name glaucoma \
--tgfm-input-data $input \
--parallel-job-identifier 1 \
--tissue-summary-file $input_summary \
--tissue-specific-prior $prior \
--gene-tissue-pip-threshold 0 \
--p-value-threshold 1 \
--out $out

done
' > /home/wangh/project1_glaucoma/analysis2/TGFM/result/sample/sample.log 2>&1 &

## TGFM final run
nohup bash -c 'source /home/wangh/miniconda3/bin/activate && \

conda activate TGFM
export R_HOME=${CONDA_PREFIX}/lib/R
export R_ENVIRON=${CONDA_PREFIX}/lib/R/etc/Renviron
export R_PROFILE=${CONDA_PREFIX}/lib/R/etc/Rprofile.site

input_base="/home/wangh/project1_glaucoma/analysis/TGFM/result/input"
input_tissue_base="/home/wangh/project1_glaucoma/analysis/TGFM/tissue_summary1_22"
tissue_prior_base="/home/wangh/project1_glaucoma/analysis/TGFM/result/tissue_prior/NA_sub"

out_base="/home/wangh/project1_glaucoma/analysis/TGFM/result/output"

for i in {1..22}
do

    chr_folder=$(printf "chr%d" $i)
    chr=$i
    
    input_file="${input_base}/${chr_folder}/_${chr_folder}_input_data_summary.txt"
    input_tissue_file="${input_tissue_base}/tissue_summary_file$chr"
    input_tissue_prior="${tissue_prior_base}/${chr_folder}"
    
    out_file="${out_base}/${chr_folder}"
    out_stem="${out_base}/${chr_folder}/${chr_folder}"
    mkdir -p "$out_file"
    
python /home/wangh/software/TGFM/run_tgfm.py \
--trait-name glaucoma \
--tgfm-input-data $input_file \
--parallel-job-identifier 1 \
--tissue-summary-file $input_tissue_file \
--tissue-specific-prior $input_tissue_prior \
--gene-tissue-pip-threshold 0 \
--p-value-threshold 1 \
--out $out_stem

done

' > /home/wangh/project1_glaucoma/analysis/TGFM/result/output/tgfm_final.log 2>&1 &

### 数据处理###
library(biomaRt)
ensembl <- useMart("ensembl") 
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl) 

gene_tissue_ids <- TGFM_gene_tissue_PIP_0$gene

gene_tissue_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = gene_tissue_ids,
                      mart = dataset)
                      
# merge
TGFM_gene_tissue_PIP_0$gene_symbol <- gene_tissue_symbols$hgnc_symbol[match(TGFM_gene_tissue_PIP_0$gene, gene_symbols$ensembl_gene_id)]

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = "109"
)

gene_positions <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  filters = "ensembl_gene_id",
  values = TGFM_gene_PIP_05_267$gene, 
  mart = ensembl
)

TGFM_gene_PIP_05_267 <- TGFM_gene_PIP_05_267 %>% left_join(gene_positions[,c(2,4,5)],by=c("gene_symbol"="hgnc_symbol"))


### 2. 

library(dplyr)
library(tidyr)
library(stats)

input_dir <- "/home/wangh/project1_glaucoma/analysis/TGFM/result/tissue_prior"
output_dir <- "/home/wangh/project1_glaucoma/analysis/TGFM/result/tissue_prior/tissue_sig"

if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

files <- list.files(input_dir, pattern = "prior_summary.txt$", full.names = TRUE)

for (file in files) {
    chr_number <- sub("chr(\\d+)_.*", "\\1", basename(file)) 
    
    df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    df$bootstrapped_prior_distribution <- as.character(df$bootstrapped_prior_distribution)
    
    bootstrapped_prior <- strsplit(df$bootstrapped_prior_distribution, ";")
    bootstrapped_prior_df <- as.data.frame(do.call(rbind, bootstrapped_prior), stringsAsFactors = FALSE)
    bootstrapped_prior_df[] <- lapply(bootstrapped_prior_df, function(x) as.numeric(x))  # 转换为数字
    
    mean_values <- apply(bootstrapped_prior_df, 1, function(x) mean(x, na.rm = TRUE))
    se_values <- apply(bootstrapped_prior_df, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    
    mean_values[apply(bootstrapped_prior_df, 1, function(x) all(is.na(x) | x == 0))] <- 0
    se_values[apply(bootstrapped_prior_df, 1, function(x) all(is.na(x) | x == 0))] <- 1
    
    z_values <- mean_values / se_values
    
    p_values <- 2 * (1 - pnorm(abs(z_values)))
    
    bh_p_values <- p.adjust(p_values, method = "fdr")
    
    df$mean <- mean_values
    df$se <- se_values
    df$Z <- z_values
    df$P <- p_values
    df$BH_P <- bh_p_values
    df$CHR <- chr_number
    
    df <- df %>% dplyr::select(-mean_prior, -bootstrapped_prior_distribution)
    
    output_file <- file.path(output_dir, basename(file))
    write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("done！\n")

