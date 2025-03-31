# Importing NCBI ClinVar 

## Table of Contents

1. [About the Dataset](#about-the-dataset)
    1. [Download URL](#download-url)
    2. [Database Overview](#database-overview)
    3. [Schema Overview](#schema-overview)
       1. [dcid Generation](#dcid-generation)
       2. [Enum Generation](#enum-generation)
       3. [Edges](#edges)
    4. [Notes and Caveats](#notes-and-caveats)
    5. [License](#license)
    6. [Dataset Documentation and Relevant Links](#dataset-documentation-and-relevant-links)
2. [About the Import](#about-the-import)
    1. [Artifacts](#artifacts)
    2. [Import Procedure](#import-procedure)
    3. [Test](#test)
    4. [Sample](#sample-data)

## About the Dataset
"The [Adult Genotype-Tissue Expression (GTEx)](https://www.gtexportal.org/home/downloads/adult-gtex/overview) project is a comprehensive public resource for the study of tissue-specific gene expression and regulation. Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals, primarily for molecular assays including WGS, WES, and RNA-Seq.

All data files generated in the GTEx project are available for download. Only open access data files are available here on the GTEx portal. All raw sequence data files (DNA and RNA), and the full donor metadata files, are protected access data. Please see our Protected Data Downloads page for instructions on how to request access to those data."

### Download URL
The latest version of single tissue QTL files can be downloaded from the appropriate [download page](https://www.gtexportal.org/home/downloads/adult-gtex/qtl) on GTEx. For this import we use the `GTEx_Analysis_v<version_number>_eQTL.tar`

### Database Overview


### Schema Overview

#### dcid Generation

#### Enum Generation

#### Edges

### Notes and Caveats

### License

### Dataset Documentation and Relevant Links

## About the Import

### Artifacts

### Import Procedure

Download the most recent versions of GTEx data. The current download.sh script is currently downloading GTEx v10 eQTL file released on November 7, 2024. Update the path to the preferred version if necessary:

```bash
sh download.sh
```

Generate the formatted CSV:

```bash
sh run.sh
```

### Tests

The first step of `tests.sh` is to downloads Data Commons's java -jar import tool, storing it in a `tmp` directory. This assumes that the user has Java Runtime Environment (JRE) installed. This tool is described in Data Commons documentation of the [import pipeline](https://github.com/datacommonsorg/import/). The relases of the tool can be viewed [here](https://github.com/datacommonsorg/import/releases/). Here we download version `0.1-alpha.1k` and apply it to check our csv + tmcf import. It evaluates if all schema used in the import is present in the graph, all referenced nodes are present in the graph, along with other checks that issue fatal errors, errors, or warnings upon failing checks. Please note that empty tokens for some columns are expected as this reflects the original data. The imports create the Virus nodes that are then referenced within this import. This resolves any concern about missing reference warnings concerning these node types by the test.

To run tests:

```bash
sh tests.sh
```

This will generate an output file for the results of the tests on each csv + tmcf pair

### Sample
