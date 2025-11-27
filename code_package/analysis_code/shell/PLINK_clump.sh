### Plink clump
# you can download MAGMA from website: https://www.cog-genomics.org/plink2/
# g1000: ref panel from EUR
# MTAG_POAG.txt: our POAG MTAG summary data

plink --bfile g1000 \
--clump MTAG_POAG.txt \
--clump-p1 5e-8 \
--clump-p2 0.01 \
--clump-r2 0.01 \
--clump-kb 1000 \
--clump-snp-field SNP \
--clump-field P \
--out plink_result
