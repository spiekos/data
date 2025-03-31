#!/bin/bash

# make output directories
mkdir -p output
mkdir -p output/egenes
mkdir -p output/sig_pairs

# Get all unique filename starts (up to the first period)
unique_starts=$(ls GTEx_Analysis_v10_eQTL | cut -d '.' -f 1 | sort | uniq)
echo "List of tissues generated!"

# Iterate through each tissue in the list
for tissue in $unique_starts; do
    gunzip -c GTEx_Analysis_v10_eQTL/$tissue.v10.eGenes.txt.gz > GTEx_Analysis_v10_eQTL/$tissue.v10.egenes.txt
    python3 scripts/format_gtex_eqtls.py \
        $tissue \
        GTEx_Analysis_v10_eQTL/$tissue.v10.egenes.txt \
        GTEx_Analysis_v10_eQTL/$tissue.v10.eQTLs.signif_pairs.parquet \
        input/gencode.v39.GRCh38.genes.only.gtf \
        input/hgnc_ensembl_mapping.txt \
        input/rsid_lookup_table.txt \
        output/egenes/$tissue.v10.egenes.csv \
        output/sig_pairs/$tissue.v10.signif_variant_gene_pairs.csv
done
