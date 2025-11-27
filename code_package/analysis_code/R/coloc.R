#### colocalization
# more detailed instruction coud be found in website:"https://github.com/chr1swallace/coloc"
# indir: contains POAG MTAG summary and eQTL summary for each gene (within 500kb)

# install the coloc r package:
install.packages("coloc")

if(!require("remotes"))
   install.packages("remotes") # if necessary
remotes::install_github("chr1swallace/coloc@main",build_vignettes=TRUE)



# coloc analysis

library(data.table)
library(coloc)

indir <- "/coloc/coloc_retina_eqtl/input"
outdir <- "/coloc/coloc_retina_eqtl"

files_gwas <- list.files(indir, pattern = "_GWAS.txt$")
genes <- sub("_GWAS.txt", "", files_gwas)
total_genes <- length(genes)
all_results <- list()

for (i in seq_along(genes)) {
  gene <- genes[i]
  
  # read
  gwas <- fread(file.path(indir, paste0(gene, "_GWAS.txt")), sep = "\t", header = TRUE)
  qtl  <- fread(file.path(indir, paste0(gene, "_QTL.txt")), sep = "\t", header = TRUE)

  gwas <- gwas[!duplicated(SNP)]
  qtl  <- qtl[!duplicated(SNP)]

  gwas <- gwas[!(BETA == 0 | VAR == 0)]
  qtl  <- qtl[!(BETA == 0 | VAR == 0)]

  # intersect
  common_snps <- intersect(gwas$SNP, qtl$SNP)
  gwas <- gwas[SNP %in% common_snps]
  qtl  <- qtl[SNP %in% common_snps]

  if (nrow(gwas) == 0 || nrow(qtl) == 0) next
  
  # input
  gwas_coloc <- list(
    beta = gwas$BETA,
    varbeta = gwas$VAR,
    snp = gwas$SNP,
    type = "cc"
  )
  
  qtl_coloc <- list(
    beta = qtl$BETA,
    varbeta = qtl$VAR,
    snp = qtl$SNP,
    type = "quant",
    sdY = 1
  )

  #  coloc
  my.res <- coloc.abf(dataset1 = gwas_coloc, dataset2 = qtl_coloc)
  res_df <- my.res$results
  res_df$gene <- gene

  all_results[[gene]] <- res_df

  cat(sprintf("[%d/%d] %.2f%% - %s done\n", 
              i, total_genes, (i/total_genes)*100, gene))
}

if (length(all_results) > 0) {
  merged_res <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  filtered_res <- merged_res[SNP.PP.H4 > 0.5]
  filtered_file <- file.path(outdir, "retina_coloc_0.5.txt")
  fwrite(filtered_res, filtered_file, sep = "\t")
}


