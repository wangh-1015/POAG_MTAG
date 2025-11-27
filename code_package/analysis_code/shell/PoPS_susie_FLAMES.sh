########## FLAMES pipeline ###########
# you can download FLAMES from website: "https://github.com/Marijn-Schipper/FLAMES"
# FLAMES pipeline need MAGMA, POPS and susie result for input!!!
# the full instruction could be found in https://github.com/Marijn-Schipper/FLAMES

### 1 MAGMA tissue analysis 
cd /MAGMA_v1.10
./magma \
--gene-results /MAGMA/Gene_analysis/magma.genes.raw \
--gene-covar /FLAMES/data/gtex_v8_ts_avg_log2TPM_entrez.txt \
--out /FLAMES/magma_tissue

### 2. PoPS analysis
cd /MAGMA_v1.10
./magma \
--bfile g1000 \
--pval MTAG_POAG_hg38.txt ncol=N \
--gene-annot /MAGMA/Gene_anotation/magma_pops.genes.annot \
--out /MAGMA/Gene_analysis/for_pops/magma_gene_pops

# pops
conda activate pops
cd /PoPS/pops-master

python pops.py \
--gene_annot_path /FLAMES/data/pops_features_full/gene_annot.txt \
--feature_mat_prefix /FLAMES/data/pops_features_full/munged_features/pops_features \
--num_feature_chunks 116 \
--magma_prefix /MAGMA/Gene_analysis/for_pops/magma_pops \
--control_features /FLAMES/data/pops_features_full/control.features \
--out_prefix /pops/POAG


### 3. SuSiE
conda activate plink

plink --bfile g1000_eur \
--extract /SuSiEx/mydata/ref/snp_list.txt \
--make-bed \
--out /SuSiEx/mydata/ref/g1000_eur

### susie

cd /SuSiEx/mydata

../bin/SuSiEx \
	--sst_file=hg38.sumstats.txt \
	--n_gwas=2317284 \
	--ref_file=./ref/g1000_eur \
	--ld_file=g1000_eur \
	--out_dir=./out \
	--out_name=SuSiEx.MTAG_locus3.cs \
	--level=0.95 \
	--pval_thresh=1e-5 \
	--chr=1 \
	--bp=87873327,88873327 \
	--snp_col=1 \
	--chr_col=2 \
	--bp_col=3 \
	--a1_col=4 \
	--a2_col=5 \
	--eff_col=6 \
	--se_col=7 \
	--pval_col=8 \
	--plink=../utilities/plink \
	--keep-ambig=True \
	--mult-step=True \
	--threads=16
	
	rm -f g1000_eur.ld.bin g1000_eur_frq.frq g1000_eur_ref.bim



# susie
nohup bash -c '

unset GREP_OPTIONS

cd /SuSiEx/mydata

SNP_FILE="/SuSiEx/mydata/loci_434_hg38.txt"

# Skip header line
sed 1d "$SNP_FILE" | while read -r snp chr bp
do
    # Calculate region ±500kb
    start=$((bp - 500000))
    end=$((bp + 500000))
    
    # Get the line number (index) for naming
    # First, find the line number of this SNP in the original file
    line_num=$(grep -n "^${snp}" "$SNP_FILE" | cut -d':' -f1)
    # Subtract 1 to account for the header
    index=$((line_num - 1))
    
    echo "Processing SNP $snp (Chr $chr: $bp) - Output will be locus_${index}.cred1"
    
    # Run SuSiEx command
    ../bin/SuSiEx \
        --sst_file=hg38.sumstats.txt \
        --n_gwas=2317284 \
        --ref_file=./ref/g1000_eur \
        --ld_file=g1000_eur \
        --out_dir=./out \
        --out_name=locus_${index}.cred1 \
        --level=0.95 \
        --pval_thresh=1e-5 \
        --chr=$chr \
        --bp=$start,$end \
        --snp_col=1 \
        --chr_col=2 \
        --bp_col=3 \
        --a1_col=4 \
        --a2_col=5 \
        --eff_col=6 \
        --se_col=7 \
        --pval_col=8 \
        --plink=../utilities/plink \
        --keep-ambig=True \
        --mult-step=True \
        --threads=16
    
    # Clean up temporary files
    rm -f g1000_eur.ld.bin g1000_eur_frq.frq g1000_eur_ref.bim
    
    echo "Completed analysis for SNP $snp"
    echo "----------------------------------------"
done

' > /SuSiEx/mydata/susiex.log 2>&1 &

### extract credible set

#!/usr/bin/env python3
import os
import glob
import pandas as pd
import sys

# Directory containing the input files
input_dir = "/SuSiEx/mydata/out"

# Get all files matching the pattern
input_files = glob.glob(os.path.join(input_dir, "locus_*.cred1.cs"))

# Track progress
total_files = len(input_files)
processed_files = 0

# Create a log file for errors
log_file = os.path.join(input_dir, "extraction_log.txt")
with open(log_file, 'w') as log:
    log.write("Extraction log:\n")

for input_file in input_files:
    # Extract the locus number from the filename
    try:
        locus_num = os.path.basename(input_file).split("_")[1].split(".")[0]
        processed_files += 1
        
        # Show progress
        sys.stdout.write(f"\rProcessing file {processed_files}/{total_files}: {os.path.basename(input_file)}")
        sys.stdout.flush()
        
        # Read the input file
        try:
            df = pd.read_csv(input_file, sep="\t", error_bad_lines=False)  # For older pandas versions
        except TypeError:
            # For newer pandas versions
            df = pd.read_csv(input_file, sep="\t", on_bad_lines='skip')
            
        # Check if the required columns exist
        required_columns = ['CS_ID', 'BP', 'REF_ALLELE', 'ALT_ALLELE', 'CS_PIP']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            with open(log_file, 'a') as log:
                log.write(f"Skipping {input_file}: Missing columns: {', '.join(missing_columns)}\n")
            continue
        
        # Get unique credible set IDs
        cs_ids = df['CS_ID'].unique()
        
        # Process each credible set
        for cs_id in cs_ids:
            # Filter dataframe for current CS_ID
            cs_df = df[df['CS_ID'] == cs_id].copy()
            
            # Skip if empty after filtering
            if cs_df.empty:
                continue
                
            # Reset index for the filtered dataframe
            cs_df = cs_df.reset_index(drop=True)
            
            # Create a new dataframe for output
            output_df = pd.DataFrame()
            
            # Add index column (1-based)
            output_df['index'] = range(1, len(cs_df) + 1)
            
            # Create the cred1 column as BP:ALT_ALLELE_REF_ALLELE
            # Check if SNP column exists
            if 'SNP' in cs_df.columns:
                output_df['cred1'] = cs_df.apply(
                    lambda row: f"{row['BP']}:{row['ALT_ALLELE']}_{row['REF_ALLELE']}" 
                    if pd.isna(row['SNP']) else f"{row['SNP']}:{row['ALT_ALLELE']}_{row['REF_ALLELE']}", 
                    axis=1
                )
            else:
                output_df['cred1'] = cs_df.apply(
                    lambda row: f"{row['BP']}:{row['ALT_ALLELE']}_{row['REF_ALLELE']}", 
                    axis=1
                )
            
            # Add the prob1 column (CS_PIP)
            output_df['prob1'] = cs_df['CS_PIP']
            
            # Create output filename
            output_file = os.path.join(input_dir, f"locus_{locus_num}_cs{cs_id}.txt")
            
            # Write to file
            output_df.to_csv(output_file, sep=" ", index=False)
            
    except Exception as e:
        with open(log_file, 'a') as log:
            log.write(f"Error processing {input_file}: {str(e)}\n")

print("\nProcessing complete!")
print(f"Processed {processed_files}/{total_files} files")
print(f"Check {log_file} for any errors")

######### add chr
library(tidyr,dplyr,data.table)
setwd("/SuSiEx/mydata/out/cred")

ref <- fread(
  "/SuSiEx/mydata/hg38.sumstats.txt"
            ,sep = "\t", header = TRUE)
ref <- ref[,c(1,2,3)]

file_names <- list.files(pattern = "locus_.*\\.txt$")

# 
for (file_name in file_names) {
  # 
  data <- read.table(file_name, header = TRUE)
  
  data <- data %>% separate(cred1, into = c("SNP", "A"), sep = ":", remove = FALSE)
  
  data <- data %>%
  left_join(ref, by = "SNP") %>% mutate(cred1=paste(CHR,BP,A,sep = ":")) %>% 
  select(index,cred1,prob1)
  
  output_file <- file.path("/SuSiEx/mydata/out/cred2", file_name)

  write.table(data, file =output_file,row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")
}


############## 4. FLAMES
# (1) annotation
conda activate FLAMES
cd /FLAMES

python /FLAMES/FLAMES-master/FLAMES.py annotate \
-o /FLAMES/annot/ \
-a /FLAMES/data/Annotation_data/ \
-p POAG.preds \
-m magma_pops.genes.out \
-mt magma_tissue.gsa.out \
-id indexfile.txt \
-b GRCh38

# (2) FLAMES score
cd /FLAMES
python /FLAMES/FLAMES-master/FLAMES.py FLAMES \
-id indexfile.txt \
-o /FLAMES




