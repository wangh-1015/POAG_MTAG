#### scDRS
# scdrs (version=1.0.2) software, you can download from
# https://martinjzhang.github.io/scDRS/
# gene_z.txt: Zscore file of our POAG MTAG
# scRNA.h5ad: .h5ad file of scRNA-seq data
# munge.gs: output file of "munge-gs"
# MTAG_POAG.txt: our POAG MTAG summary data
# POAG.full_score.gz: output file of "compute score"
# cell_type: metadata in scRNA.h5ad describing the cell type of each cell

# install scDRS
git clone https://github.com/martinjzhang/scDRS.git
cd scDRS
pip install -e .

# 1. munge-gs
# Select top 1,000 genes and use z-score weights

scdrs munge-gs \
--out-file munge.gs \
--zscore-file gene_z.txt \
--weight zscore \
--n-max 1000

# 2. compute score

scdrs compute-score \
--h5ad-file scRNA.h5ad \
--h5ad-species human \
--gs-file munge.gs \
--gs-species human \
--out-folder scDRS_result/POAG \
--flag-filter-data True \
--flag-raw-count True \
--n-ctrl 1000 \
--flag-return-ctrl-raw-score False \
--flag-return-ctrl-norm-score True

  
# 3. perform downstream analysis

scdrs perform-downstream \
--h5ad-file scRNA.h5ad \
--score-file scDRS_result/POAG.full_score.gz \
--out-folder scDRS_result/downstream_result \
--group-analysis cell_type \
--gene-analysis \
--flag-filter-data FALSE \
--flag-raw-count True
