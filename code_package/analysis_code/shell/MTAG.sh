#### MTAG

# MTAG/mtag.py: main director of MTAG software, you can download from
# https://github.com/JonJala/mtag
# sumstats: four sources of GWAS summary data

 git clone https://github.com/omeed-maghzian/mtag.git
 cd mtag
./mtag.py -h

python mtag.py \
--sumstats MVP_POAG.txt,Finn_POAG.txt,GBMI_POAG.txt,prev_MTAG.txt \
--out MTAG_POAG \
--n_min 0.0 \
--incld_ambig_snps \
--stream_stdout \
--snp_name SNP \
--chr_name CHR \
--bpos_name BP \
--a1_name A1 \
--a2_name A2 \
--eaf_name FRQ \
--beta_name BETA \
--se_name SE \
--z_name Z \
--p_name P \
--n_name N

