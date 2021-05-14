# Gene Category Enrichment Analysis including Custom Null Ensembles

[![DOI](https://zenodo.org/badge/79196471.svg)](https://zenodo.org/badge/latestdoi/79196471)

This is a Matlab toolbox for performing gene-category enrichment analysis relative to two different types of null models:
1. ___Random-gene nulls___, in which categories assessed relative to categories of the same size but annotated by the same number of random genes.
   This follows the permutation-based method of Gene Score Resampling (as implemented in [*ermineJ*](https://erminej.msl.ubc.ca/)).
2. ___Ensemble-based nulls___, in which categories are assessed relative to an ensemble of randomized phenotypes.

The toolbox was introduced in our paper:
- Fulcher et al. _Nature Communications_ (2021) [:green_book: 'Overcoming false-positive gene-category enrichment in the analysis of spatially resolved transcriptomic brain atlas data'](https://doi.org/10.1038/s41467-021-22862-1).

Instructions for performing the basic functions of these analyses are in [the wiki :notebook:](https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis/wiki).

The package is currently set up to perform enrichment on [Gene Ontology](http://geneontology.org/) (GO) Biological Process annotations, but could be modified straightforwardly to use other types of GO annotations, or even to use other annotation systems like [KEGG](https://www.genome.jp/kegg/).

Pull requests to improve the functionality and clarity of documentation are very welcome!

#### Repository Organization

The package is organized into directories as follows:

__DATA__:
1. `RawData`: all data downloaded from external sources (like GO, MouseMine, etc.)
2. `ProcessedData`: raw data processed into Matlab-readable files.

__CODE__:
1. `DataProcessing`: code required to process raw data.
2. `GeneScoreResampling`, `EnsembleEnrichment`: code to run both random-gene and randomized-phenotype enrichment analysis.
3. `ResultsComparison`: code to compare GSEA results to _ermineJ_.
4. `Peripheral`: additional code files.

To initialize this toolbox, all of these subdirectories should be added to the Matlab path by running the `startup` script.

## Running analysis

A summary of how to run an enrichment analysis with this package is describd here, but please read the [wiki :notebook:](https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis/wiki) for more detailed instructions.

### Preparation: Defining gene-to-category annotations

The first step in running an enrichment analysis is defining the set of gene categories, and the genes annotated to each category.
Results of this, using hierarchy-propagated gene-to-category annotations corresponding to GO biological processes (processed on 2019-04-17), can be downloaded from [this partner Zenodo data repository](https://doi.org/10.5281/zenodo.4460713).

Code in this repository also allows you to reprocess these annotations from raw data from GO, as described on [this wiki page](https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis/wiki/Defining-gene-to-category-annotations).
You can test this pipeline using the `term` and `term2term` tables from a mySQL download of the GO term data on 2019-04-17, which are also available in the associated [Zenodo data repository](https://doi.org/10.5281/zenodo.4460713).

### Performing Enrichment

All parameters are set using `GiveMeDefaultEnrichmentParams`, as described in the [wiki](https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis/wiki/Setting-parameters).

#### Gene-score resampling (random-gene null)

The Gene Score Resampling method assesses significance relative to a 'random-gene null', and is implemented in the `SingleEnrichment` function.
Instructions to implement this are in the [wiki](https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis/wiki/Random-gene-enrichment).

#### Ensemble enrichment

Ensemble enrichment computes the enrichment of a given phenotype relative to an ensemble of randomized phenotypes, as described in [our bioRxiv preprint](https://doi.org/10.1101/2020.04.24.058958).

This proceeds across `ComputeAllCategoryNulls` (precompute category nulls) and `EnsembleEnrichment` (evaluate significance relative to these nulls), as described in the [wiki](https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis/wiki/Ensemble-enrichment).
