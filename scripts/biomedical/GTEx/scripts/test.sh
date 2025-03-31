#!/bin/bash

# download data commons java test tool version 0.1-alpha.1k
rm -rf tmp
mkdir -p tmp; cd tmp
wget https://github.com/datacommonsorg/import/releases/download/0.1-alpha.1k/datacommons-import-tool-0.1-alpha.1-jar-with-dependencies.jar
cd ..
mkdir -p reports
mkdir -p reports/egenes
mkdir -p reports/sig_pairs

# Get all unique filename starts (up to the first period)
unique_starts=$(ls GTEx_Analysis_v10_eQTL | cut -d '.' -f 1 | sort | uniq)
echo "List of tissues generated!"

# Iterate through each tissue in the list
for tissue in $unique_starts; do
	java -jar tmp/datacommons-import-tool-0.1-alpha.1-jar-with-dependencies.jar lint tMCFs/gtex_egenes.tmcf output/egenes/$tissue.*.egenes.csv *.mcf
	mv dc_generated reports/egenes/$tissue

	java -jar tmp/datacommons-import-tool-0.1-alpha.1-jar-with-dependencies.jar lint tMCFs/gtex_sig_pairs.tmcf output/sig_pairs/$tissue.*.signif_variant_gene_pairs.csv *.mcf
	mv dc_generated reports/sig_pairs/$tissue
done
