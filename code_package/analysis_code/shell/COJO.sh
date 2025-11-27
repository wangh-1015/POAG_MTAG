### COJO conditional analysis: 
# loci_495_hg38.txt have 3 columns: SNP CHR BP, which are used to perform
# COJO analysis on each locus (within 500bp of lead SNP)
# g1000: ref panel from EUR
# maf:  maf threshold

row_number=1
while read -r rs chrom pos; do
    echo "Reading: $rs $chrom $pos"
    
    output_file="locus$row_number"

    ./gcta64 \
    --bfile g1000 \
    --maf 0.001 \
    --cojo-file MTAG_POAG.txt \ # POAG MTAG summary data
    --extract-region-snp "$rs" 1000 \
    --cojo-cond snp_list_495.txt \ # conditional SNP list
    --thread-num 8 \
    --out "$output_file"
    
    ((row_number++))
done < loci_495_hg38.txt
