# TWAS  
# install FUSION, download ref panel
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip master.zip
cd fusion_twas-master

wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar xjvf LDREF.tar.bz2

wget https://github.com/gabraham/plink2R/archive/master.zip
unzip master.zip

# POAG_MTAG.txt: our POAG MTAG summary data, have four columns "SNP A1(effect allele) A2(other allele) Z"
# you need to download retina expression weight file from website: "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115828"

# TWAS analysis
cd /fusion_twas-master
for chr in {1..22}
do
    Rscript FUSION.assoc_test.R \
    --sumstats POAG_MTAG.txt \
    --weights /eyeGTEx/final.pos \
    --weights_dir /eyeGTEx/Rdata \
    --ref_ld_chr /fusion_twas-master/LDREF_hg38/1000G.EUR. \
    --chr "$chr" \
    --out TWAS_POAG."$chr".dat
done







