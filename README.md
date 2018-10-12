# Gene Enrichment

This is a Matlab repository for performing gene ontology enrichment given a set of gene scores.

Some information and usage examples are below:

## Data Processing

### Retrieving GO Terms

Assumes you have an SQL database setup, with details entered into `GetGOTerms`

Get biological process GO terms:
```matlab
GOTerms = GetGOTerms('biological_process')
```

### Generating formatted files processing raw GO annotations

* `ReadDirectAnnotationFile` to read in raw annotations from GO (or `ReadGEMMAAnnotationFile` to read in processed annotations from GEMMA)
* `propagateHierarchy` to propagate annotations up the hierarchy

## Analyses
### Conventional enrichment using `SingleEnrichment`
Conventional enrichment results can be obtained using `SingleEnrichment`.
This analysis considers a null model in which scores can be assigned to individual genes at random.

INPUTS:
* `geneScores`, a numGenes-long column vector of values that quantifies something about each gene.
* `geneEntrezIDs`, numGenes-long column vector labeling the entrez ID for each gene in geneScores.
* `dataSource`, specifies the source of GO annotations, to be loaded using `GetFilteredGOData`. Options are `mouse-direct` (hierarchy and annotations taken directly from GO), `human-direct` (hierarchy and annotations taken directly from GO), `mouse-GEMMA` (processed hierarchy and annotations downloaded from GEMMA).
* `processFilter`, what GO processes to consider. Default is `biological_process`.
* `sizeFilter`, filter GO categories by size. Default is `[5,200]` (only consider categories with between 5 and 200 annotations).
* `numIters`, number of iterations (`1e4` is the default, can ramp up to get better significance estimates for small p-values).

```matlab
GOTable = SingleEnrichment(geneScores,geneEntrezIDs,dataSource,processFilter,sizeFilter,numIters)
```
