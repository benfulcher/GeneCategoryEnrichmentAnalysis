# Gene Category Enrichment Analysis

This is a Matlab toolbox for performing gene category enrichment analysis relative to two different types of null models:
1. ___Random-gene nulls___, in which categories assessed relative to categories of the same size but annotated by the same number of random genes.
   This follows the permutation-based method of Gene Score Resampling (as implemented in [*ermineJ*](https://erminej.msl.ubc.ca/)).
2. ___Ensemble-based nulls___, in which categories are assessed relative to an ensemble of null phenotypes, as introduced in [this bioRxiv preprint](https://doi.org/10.1101/2020.04.24.058958).

The package is currently set up to perform enrichment on [Gene Ontology](http://geneontology.org/) (GO) Biological Process annotations, but could be modified in future to use other GO annotations, or use other annotation systems (like [KEGG](https://www.genome.jp/kegg/)).

Pull requests to improve the functionality and clarity of documentation are very welcome!

#### Repository Organization
The package is organized into directories as follows:

__DATA__:
1. `RawData`: all data downloaded from external sources (like GO, MouseMine, etc.)
2. `ProcessedData`: raw data processed into Matlab-readable files.

__CODE__:
1. `DataProcessing`: code required to process raw data.
2. `Analysis`: code to run enrichment analysis.
3. `ResultsComparison`: code to compare GSEA results to _ermineJ_.
4. `Peripheral`: additional code files.

To initialize this toolbox, all of these subdirectories should be added to the Matlab path by running the `startup` script.

## Preparation: Defining gene-to-category annotations

The first step in running an enrichment analysis is defining the set of gene categories, and the genes annotated to each category.
Results of this, using hiearchy-propagated gene-to-category annotations corresponding to GO biological processes (processed as of 2019-04-17), can be downloaded from [this figshare repository](https://figshare.com/s/71fe1d9b2386ec05f421).

Code in this repository also allows you to reprocess these annotations from raw data from GO, as described on [this wiki page](https://github.com/benfulcher/GeneSetEnrichmentAnalysis/wiki/Defining-gene-to-category-annotations).

## Analysis
### Random-Gene Null (Gene Score Resampling)
The Gene Score Resampling method for GSEA is performed using the function `SingleEnrichment`.
This requires GO Terms and annotations to be processing (as described in the Initialization section below).
To assess significance of gene scores annotated to a specific GO category, this analysis compares results to a null model in which scores are assigned to genes at random.

__INPUTS__:
* `geneScores`, a numGenes-long column vector of values that quantifies something about each gene.
* `geneEntrezIDs`, numGenes-long column vector labeling the entrez ID for each gene in geneScores.
Extra options set as fields in a parameter structure:
* `dataSource`, specifies the source of GO annotations, to be loaded using `GetFilteredGOData`.
Options are `mouse-direct` (hierarchy and annotations taken directly from GO), `human-direct` (hierarchy and annotations taken directly from GO), `mouse-GEMMA` (processed hierarchy and annotations downloaded from GEMMA).
* `processFilter`, a subset of GO processes to consider.
Default: `biological_process`.
* `sizeFilter`, filter to a subset of GO categories by their number of annotations ('size').
Default is `[5,200]` (only consider GO categories with between 5 and 200 gene annotations).
* `numSamples`, number of iterations (`1e4` is the default, can ramp up to get better significance estimates for small p-values).

__EXAMPLE USAGE__:
```matlab
enrichmentSettings = GiveMeDefaultEnrichmentParams();
GOTable = SingleEnrichment(geneScores,geneEntrezIDs,enrichmentSettings);
```

The output `GOTable` sorts GO categories by their estimated _p_-value.
Note that _p_-values are estimated according to two different methods:
1. `pValPerm`: _p_-value from a permutation test.
2. `pValZ`: _p_-value estimated from a Gaussian fit to the null distribution.

Both _p_-value estimates are corrected using the method of false discovery rate (Benjamini and Hochberg), in the corresponding columns `pValPermCorr` and `pValZCorr`.

### Ensemble Enrichment
Ensemble enrichment computes the enrichment of a given phenotype relative to an ensemble of randomized phenotypes.

It proceeds through two steps:
1. Compute the ensemble using `ComputeAllCategoryNulls`
2. Perform enrichment for a phenotype of interest relative to this null ensemble using `EnsembleEnrichment`

Default parameters for enrichment (including those relevant to ensemble enrichment) are set with `GiveMeDefaultEnrichmentParams`.

__EXAMPLE USAGE__:
```matlab
% Set parameters for the calculation
% (alter parameters within this file to ensure appropriate output filename):
enrichmentParams = GiveMeDefaultEnrichmentParams();

% Compute category score null distributions resulting from a given null phenotype ensemble:
% (saves results out to a filename, enrichmentParams.fileNameOut):
ComputeAllCategoryNulls(geneDataStruct,enrichmentParams,[],true,true);

% Compares a given phenotype to the saved null ensemble
GOTablePhenotype = EnsembleEnrichment(enrichmentParams.fileNameOut,phenotypeVector)
```
