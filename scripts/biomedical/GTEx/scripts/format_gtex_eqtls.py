# Copyright 2024 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Author: Samantha Piekos
Date: 02/15/2022
Name: format_gtex_eqtls.py
Edited By: Samantha Piekos
Last Edited: 11/11/24
Description: Converts egenes and signif_variant_gene_pair from GTEx into their
own CSV files with paired tmcf files for import into Data Commons knowledge
graph. Additional lookup files are needed to convert Ensembl gene ids to gene
symbol and GTEx variant ids to rsIDs. Note this script needs to be run in parallel
as the adult GTEx has these paired files in 54 distinct tissues. This code was
written for v8 which was conducted using Genome Assembly GRCh38.p10. This assembly
is used to correspond the chromosomes with the appropriate genome_accession.version.
@tissue                 the tissue in which all the statistical associations between
                        genes and genetic variants were observed
@file_egenes            the GTEx egene tsv file
@file_sig_pairs         the GTEx signif_variant_gene_pairs tsv file
@file_gene_lookup       lookup table to convert ensembl gene id to gene symbol,
                        this was generated from data from NCBI Gene
@file_rsid_lookup       lookup table to convert GTEx variant ids into rsIDs, this
                        is a GTEx auxillary file
@file_output_egenes     formatted egenes csv file
@file_output_sig_pairs  formated signif_variant_gene_pairs csv file
"""

# set up environment
import csv
import pandas as pd
import sys


# declare universal variables
# convert biotype column to TypeOfGeneEnum
BIOTYPE_DICT = {
    'IG_C_gene': 'dcs:TypeOfGeneIGCGene',
    'IG_C_pseudogene': 'dcs:TypeOfGeneIGCPseudogene',
    'IG_J_gene': 'dcs:TypeOfGeneIGJGene',
    'IG_V_gene': 'dcs:TypeOfGeneIGVGene',
    'TEC': 'dcs:TypeOfGeneToBeExperimentallyConfirmed',
    'TR_C_gene': 'dcs:TypeOfGeneTRCGene',
    'TR_V_gene': 'dcs:TypeOfGeneTRVGene',
    'TR_V_pseudogene': 'dcs:TypeOfGeneTRVPseudogene',
    'lncRNA': 'dcs:TypeOfGenelncRNA',
    'miRNA': 'dcs:TypeOfGenemiRNA',
    'misc_RNA': 'dcs:TypeOfGenemiscRNA',
    'polymorphic_pseudogene': 'dcs:TypeOfGenePolymorphicPseudogene',
    'processed_pseudogene': 'dcs:TypeOfGeneProcessedPseudogene',
    'protein_coding': 'dcs:TypeOfGeneProteinCoding',
    'pseudogene': 'dcs:TypeOfGenePseudo',
    'rRNA_pseudogene': 'dcs:TypeOfGenerRNAPseudogene',
    'scRNA': 'dcs:TypeOfGenescRNA',
    'scaRNA': 'dcs:TypeOfGenescaRNA',
    'snRNA': 'dcs:TypeOfGenesnRNA',
    'snoRNA': 'dcs:TypeOfGenesnoRNA',
    'transcribed_processed_pseudogene': 'dcs:TypeOfGeneTranscribedProcessedPseudogene',
    'transcribed_unitary_pseudogene': 'dcs:TypeOfGeneTranscribedUnitaryPseudogene',
    'transcribed_unprocessed_pseudogene': 'dcs:TypeOfGeneTranscribedUnprocessedPseudogene',
    'translated_processed_pseudogene': 'dcs:TypeOfGeneTranslatedProcessedPseudogene',
    'translated_unprocessed_pseudogene': 'dcs:TypeOfGeneTranslatedUnprocessedPseudogene',
    'unitary_pseudogene': 'dcs:TypeOfGeneUnitaryPseudogene',
    'unprocessed_pseudogene': 'dcs:TypeOfGeneUnprocessedPseudogene'  
}


# use to convert +/- to StrandOrientationEnum
STRAND_ORIENTATION_DICT = {
    '+': 'dcs:StrandOrientationPositive',
    '-': 'dcs:StrandOrientationNegative'
}


# declare functions
def parse_gtf(file_path):
  """
  Reads a GTF file into a pandas DataFrame and expands the last column 
  into separate columns based on the provided format.

  Args:
    file_path: Path to the GTF file.

  Returns:
    A pandas DataFrame with the parsed GTF data.
  """

  df = pd.read_csv(file_path, sep='\t', header=None, comment='#')

  # Split the last column by semicolon and create a list of key-value pairs
  df[len(df.columns)] = df[df.columns[-1]].str.strip('; ').str.split('; ').apply(
      lambda x: {item.split(' ')[0]: item.split(' ')[1].strip('"') for item in x}
  )

  # Expand the dictionary into separate columns
  df = pd.concat([df, df[len(df.columns)-1].apply(pd.Series)], axis=1)

  # keep only the needed columns
  df = df[['gene_id', 'gene_name', 'hgnc_id']]

  return df


'''
def create_variant_rsid_dict(file_path):
  """
  Reads a GTEx variant annotation file and creates a dictionary mapping variant_id to rsID,
  skipping rows with missing rsIDs.

  Args:
    file_path: Path to the GTEx variant annotation file.

  Returns:
    A dictionary where keys are variant_ids and values are rsIDs.
  """

  variant_rsid_dict = {}
  with open(file_path, 'r') as file:
    reader = csv.reader(file, delimiter=' ')
    next(reader)  # Skip the header row

    for row in reader:
      variant_id = row[0]
      rs_id = row[1]
      if rs_id != '.':  # Skip rows where rs_id is missing
        variant_rsid_dict[variant_id] = rs_id

  return variant_rsid_dict
  '''


def create_variant_rsid_dict(file_path):
    """
    Reads a GTEx variant annotation file and creates a pandas df mapping variant_id to rsID,
    skipping rows with missing rsIDs.

    Args:
        file_path: Path to the GTEx variant annotation file.

    Returns:
        A pandas df of variant_ids and corresponding rsIDs.
    """
    df = pd.read_csv(file_path, sep=' ', names=["variant_id", "rsID"], header=0)
    df = df[df['rsID'] != '.']  # Remove rows with '.' in 'rsID'
    return df


def format_strand_orientation_enum(df):
    # Convert strand column to StrandOrientationEnum
    df['strand'] = df['strand'].map(STRAND_ORIENTATION_DICT).fillna(df['strand'])
    return df


def format_type_of_gene_enum(df):
    # Convert biotype column to TypeOfGeneEnum
    df['biotype'] = df['biotype'].map(BIOTYPE_DICT).fillna('dcs:TypeOfGeneUnknown')
    return df


def add_alternate_name(df, col1, col2, col3):
  """
  Adds an 'alternate_name' column to a DataFrame based on the comparison 
  of two existing columns.

  Args:
    df: The input DataFrame.
    col1: The name of the first column to compare.
    col2: The name of the second column to compare.

  Returns:
    The DataFrame with an added 'alternate_name' column.
  """

  df['alternate_name'] = df.apply(lambda row: row[col1] if row[col1] != row[col2] and row[col1] != row[col3] else '', axis=1)
  
  # drop the columns that are now redundan
  df = df.drop([col1, col2], axis=1)

  return df


def add_gene_symbol(df, df_gene_lookup, df_hgnc_lookup):
    # Add gene symbol to df using df_gene_lookup matching on ensembl_gene_id
    df = pd.merge(df, df_gene_lookup, on='gene_id', how='left')

    # Add gene symbol by joining on hgnc id with the hgnc lookup table
    df = pd.merge(df, df_hgnc_lookup, on='hgnc_id', how='left')

    # Fill missing 'gene_symbol' values with 'gene_name'
    df['symbol'] = df['symbol'].fillna(df['gene_name_y'])
    df['symbol'] = df['symbol'].fillna(df['gene_name_x'])

    # identify when there is a mismatch in gene_name between GTEx files
    # note the one used by egenes as an alternate name for the gene
    df = add_alternate_name(df, 'gene_name_x','gene_name_y', 'symbol')

    return df


def check_for_illegal_charc(s):
    """Checks for illegal characters in a string and prints an error statement if any are present
    Args:
        s: target string that needs to be checked
    
    """
    list_illegal = ["'", "*" ">", "<", "@", "]", "[", "|", ":", ";" " "]
    if any([x in s for x in list_illegal]):
        print('Error! dcid contains illegal characters!', s)
    return


def format_association_dcids(df, tissue, rs_col):
    # genearte dcid for GeneGeneticVariantAssociation node
    df['dcid'] = 'bio/' + df['symbol'].replace('@', '_Cluster') + '_' + df[rs_col] + '_' + tissue

    # check generated dcids for illegal characters
    df['dcid'].apply(check_for_illegal_charc)

    # generate name for GeneGeneticVariantAssociation node
    df['name'] = df['symbol'] + ' ' + df[rs_col] + ' Association In ' + tissue.replace('_', ' ')

    return df


def format_gene_genetic_variant_dcids(df, tissue, rs_col):
    # generate dcid to reference gene from the gene symbol
    df['dcid_gene'] = 'bio/' + df['symbol']
    df['dcid_gene'] = df['dcid_gene'].replace('@', '_Cluster')
    
    # generate dcid to reference the genetic variant from the rsID
    df['dcid_variant'] = 'bio/' + df[rs_col]

    # check generated dcids for illegal characters
    df['dcid_gene'].apply(check_for_illegal_charc)
    df['dcid_variant'].apply(check_for_illegal_charc)

    # genearte dcid for GeneGeneticVariantAssociation node
    format_association_dcids(df, tissue, rs_col)

    return df


def format_genomic_coordinate_dcids(df):
    # generate gene GenomicCoordinate dcids and names
    df['dcid_gene_coordinates'] = 'bio/GRCh38.p13_' + df['symbol'] + '_coordinates'
    df['name_gene_coordinates'] = 'GRCh38.p13 ' + df['symbol'] + ' Coordinates'

    # generate genetic variant GenomicPosition dcids and names
    df['dcid_variant_pos'] = 'bio/GRCh38.p13_' + df['chr'] + '_' + df['variant_pos'].astype(str) 
    df['name_variant_pos'] = 'GRCh38.p13 ' + df['chr'] + ' ' + df['variant_pos'].astype(str) 

    # check generated dcids for illegal characters
    df['dcid_gene_coordinates'].apply(check_for_illegal_charc)
    df['dcid_variant_pos'].apply(check_for_illegal_charc)

    return df


def format_egenes(df_egenes, df_gene_lookup, df_hgnc_lookup, tissue):
    # add columns that indicate tissue
    df_egenes['Tissue'] = tissue.replace('_', ' ')

    # add column containing just the ensembl id of the gene
    df_egenes['ensembl_id'] = df_egenes['gene_id'].str.split('.').str[0]

    # convert strand orientation +/- to enum
    df_egenes = format_strand_orientation_enum(df_egenes)

    # convert biotype to TypeOfGeneEnum
    df_egenes = format_type_of_gene_enum(df_egenes)

    # map gene symbol to the provided ensembl gene id
    df_egenes = add_gene_symbol(df_egenes, df_gene_lookup, df_hgnc_lookup)

    # GTEx variant ids can map to more than one rsID
    # expand the 'rs_id_dbSNP155_GRCh38p13' so there is one rsID per row
    # Split the 'rs_id_dbSNP155_GRCh38p13' column by comma
    df_egenes['rs_id_dbSNP155_GRCh38p13'] = df_egenes['rs_id_dbSNP155_GRCh38p13'].str.split(',')
    expanded_df = df_egenes.explode('rs_id_dbSNP155_GRCh38p13')

    # format dcids for each unique referenced entity
    expanded_df = format_gene_genetic_variant_dcids(expanded_df, tissue, 'rs_id_dbSNP155_GRCh38p13')

    # format genomic coordinate dcids
    expanded_df = format_genomic_coordinate_dcids(expanded_df)

    # track progress
    print(tissue.replace('_', ' '), 'eGenes File Processed!')

    return expanded_df


def format_sig_pairs(df, df_gene_dcid_lookup, df_rsid_lookup, tissue):
    # format tissue name
    df['Tissue'] = tissue.replace('_', ' ')

    # add column containing just the ensembl id of the gene
    df['ensembl_id'] = df['gene_id'].str.split('.').str[0]

    # add gene dcid and symbol using egenes df
    df = pd.merge(df, df_gene_dcid_lookup, on='gene_id', how='left')

    # add variant dcid and rsID using egenes df
    df = pd.merge(df, df_rsid_lookup, on='variant_id', how='left')

    # drop rows whose variants don't map to an rsID
    df.dropna(subset=['rsID'], inplace=True)

    # format dcid for sig pairs nodes
    df = format_gene_genetic_variant_dcids(df, tissue, 'rsID')
    
    # track progress
    print(tissue.replace('_', ' '), 'Sig Pairs File Processed!')
    return df


def main():
    # read in arguments
    tissue = sys.argv[1]
    file_egenes = sys.argv[2]
    file_sig_pairs = sys.argv[3]
    file_gene_lookup = sys.argv[4]
    file_hgnc_lookup = sys.argv[5]
    file_rsid_lookup = sys.argv[6]
    file_output_egenes = sys.argv[7]
    file_output_sig_pairs = sys.argv[8]

    # print progress
    print('Processing', tissue.replace('_', ' '), 'Files!')

    # open input files
    df_egenes = pd.read_csv(file_egenes, sep='\t')
    #df_sig_pairs = pd.read_csv(file_sig_pairs, sep='\t')
    df_sig_pairs = pd.read_parquet(file_sig_pairs, engine='pyarrow')
    df_gene_lookup = parse_gtf(file_gene_lookup)
    df_hgnc_lookup = pd.read_csv(file_hgnc_lookup, sep=' ')
    df_rsid_lookup = create_variant_rsid_dict(file_rsid_lookup)

    # format files
    df_egenes = format_egenes(df_egenes, df_gene_lookup, df_hgnc_lookup, tissue)
    df_gene_dcid_lookup = df_egenes[['gene_id', 'symbol']]
    df_sig_pairs = format_sig_pairs(df_sig_pairs, df_gene_dcid_lookup, df_rsid_lookup, tissue)

    # save formatted files
    df_egenes.to_csv(file_output_egenes, doublequote=False, escapechar='\\', index=False)
    df_sig_pairs.to_csv(file_output_sig_pairs, doublequote=False, escapechar='\\', index=False) 
    print('Finished Formatting', tissue.replace('_', ' '), 'Files!\n')


if __name__ == '__main__':
    main()
  