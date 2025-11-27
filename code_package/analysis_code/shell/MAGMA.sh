############# MAGMA 
# MAGMA_v1.10: main director of MAGMA software, you can download from
# https://cncr-nl.ontw.stuurlui.dev/research/magma/
# MTAG_POAG.txt: our POAG MTAG summary data
# g1000: ref panel from EUR
# magma_anno.genes.annotz: output file of "gene anotation"
# magma.genes.raw:  output file of "gene-level analysis"

### 1. gene anotation 
cd MAGMA_v1.10
./magma \
--annotate window=10 MTAG_POAG.txt \
--snp-loc MTAG_POAG.txt \
--gene-loc /dbSNP38/NCBI38.gene.loc \
--out magma_anno

####### 2. gene-level analysis(p-value) 
./magma \
--bfile g1000 \
--pval MTAG_POAG.txt ncol=N \
--gene-annot magma_anno.genes.annot \
--out magma_genelevel

### 3. pathway analysis
# /home/wangh/software/MAGMA_v1.10/MSigDB/msigdb_v2024.1.Hs_GMTs/c2.all.v2024.1.Hs.entrez.gmt
cd MAGMA_v1.10
  ./magma \
  --gene-results magma.genes.raw \
	--set-annot /MSigDB/msigdb_v2024.1.Hs_GMTs/c2.all.v2024.1.Hs.entrez.gmt \
	--out MAGMA_geneset


