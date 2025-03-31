#!/bin/bash

# download GTEx eQTL data
curl -o GTEx_Analysis_v10_eQTL.tar https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/GTEx_Analysis_v10_eQTL.tar
tar -xf GTEx_Analysis_v10_eQTL.tar
rm GTEx_Analysis_v10_eQTL.tar

# make sure that there is no trailing info to the folder
mv GTEx_Analysis_v10_eQTL_* GTEx_Analysis_v10_eQTL

# download GTEx genes.gtf file
mkdir input
cd input
curl -o gencode.v39.GRCh38.genes.gtf https://storage.googleapis.com/adult-gtex/references/v10/reference-tables/gencode.v39.GRCh38.genes.gtf
awk '$3 == "gene"' gencode.v39.GRCh38.genes.gtf > gencode.v39.GRCh38.genes.only.gtf
rm gencode.v39.GRCh38.genes.gtf

# download GTEx lookup_table
curl -o GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz https://storage.googleapis.com/adult-gtex/references/v10/reference-tables/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz
gunzip GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz
awk -F '\t' 'NR==FNR { print $1, $7 }' GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt > rsid_lookup_table.txt
rm GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt

# download HGNC ID, gene symbol, ensembl ID mapping table from HGNC
curl -o hgnc_complete_set.txt storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
awk -F  '\t' 'NR==FNR { print $1, $2 }' hgnc_complete_set.txt > hgnc_ensembl_mapping.txt
rm hgnc_complete_set.txt