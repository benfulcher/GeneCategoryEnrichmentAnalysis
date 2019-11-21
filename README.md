# Gene Set Enrichment Analysis by Gene Score Resampling

This is a Matlab repository for performing Gene Set Enrichment Analysis (GSEA) on Gene Ontology (GO) Terms.
Given a set of scores assigned to genes, we use the permutation-based method of Gene Score Resampling (GSR), following the [*ermineJ*](https://erminej.msl.ubc.ca/) software package.
The package supports conventional enrichment, which generates a null distribution for some summary statistic across the scores assigned to genes annotated to a given GO Term by randomizing the assignment of genes to GO categories.
These null distributions depend on the size of the GO Category.

There are two steps in performing an analysis:
1. Initialize the package with the latest GO hierarchy and gene annotations (or download processed results from [figshare](https://figshare.com/s/71fe1d9b2386ec05f421)).
2. Run an enrichment analysis.

## Initialization

We first describe how to retrieve and process data from GO in a form that facilitates GSEA.
Note that this step can be skipped by downloading pre-processed results from [figshare](https://figshare.com/s/71fe1d9b2386ec05f421).

This involves the following steps:
1. Retrieve and process the GO hierarchy data.
2. Retrieve and process the annotations of genes to GO Terms.
3. Iteratively propagate gene-to-Term annotations from child to parent up the GO hierarchy.

### Retrieving the GO-Term Hierarchy
There are a number of routes to [downloading the GO Term hierarchy](http://geneontology.org/page/download-ontology).
We used the **termdb** [mySQL database dump](http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz), and link to this database using a mySQL java connector, as implemented in the [Matlab_mySQL repository](https://github.com/benfulcher/Matlab_mySQL).

The data is provided in raw form as `go-basic.obo` (`basic` file ensures that annotations can be propagated), and you can also download it as a [database](ftp://ftp.geneontology.org/go/www/GO.downloads.database.shtml).

Assumes you have an SQL database setup, with details entered into `GetGOTerms`.

Get biological process GO terms and save the filtered set of GO terms to file:
```matlab
GOTerms = GetGOTerms('biological_process',true);
```
Saves out to `ProcessedData/GOTerms_BP.mat`.

### Generating formatted files for Matlab, processed from raw GO annotations



#### Downloading raw GO annotation data

You can download annotation files direct from the [GO website](http://current.geneontology.org/products/pages/downloads.html).
For _Mus musculus_, this yields: `mgi.gaf`.
For _Homo sapiens_, it is: `goa_human.gaf`.
Data used here is from the 2019-04-17 release.

#### Read in from the GO annotation file:
```matlab
ReadDirectAnnotationFile('mouse')
```
Saves processed data in the form `GOAnnotationDirect-mouse.mat`, in the `ProcessedData` directory.

Note that processed annotations from [GEMMA](https://gemma.msl.ubc.ca/annots/) can alternatively be read in using `ReadGEMMAAnnotationFile`.

#### Propagate annotations up through the hierarchy
Annotations are made at the lowest level of the GO term hierarchy.
It therefore makes sense to propagate direct annotations up to all the parent terms.
This can be achieved using:

```matlab
propagateHierarchy('mouse','biological_process');
```
Loads in from the previous step (e.g., `GOAnnotationDirect-mouse.mat`) and saves output as `GOAnnotationDirect-mouse-biological_process-Prop.mat`.
This information can then be read in for enrichment or other analyses.

Note that this script requires access to a mySQL database containing the hierarchical GO term information.

## Analyses
### Conventional enrichment using `SingleEnrichment`
Conventional enrichment results can be obtained using `SingleEnrichment`.
This analysis considers a null model in which scores can be assigned to individual genes at random.

INPUTS:
* `geneScores`, a numGenes-long column vector of values that quantifies something about each gene.
* `geneEntrezIDs`, numGenes-long column vector labeling the entrez ID for each gene in geneScores.
Extra options set as fields in a parameter structure:
* `dataSource`, specifies the source of GO annotations, to be loaded using `GetFilteredGOData`. Options are `mouse-direct` (hierarchy and annotations taken directly from GO), `human-direct` (hierarchy and annotations taken directly from GO), `mouse-GEMMA` (processed hierarchy and annotations downloaded from GEMMA).
* `processFilter`, what GO processes to consider. Default is `biological_process`.
* `sizeFilter`, filter GO categories by size. Default is `[5,200]` (only consider categories with between 5 and 200 annotations).
* `numSamples`, number of iterations (`1e4` is the default, can ramp up to get better significance estimates for small p-values).

```matlab
enrichmentSettings = struct();
enrichmentSettings.dataSource = 'mouse-direct';
enrichmentSettings.processFilter = 'biological_process';
enrichmentSettings.sizeFilter = [5,100];
enrichmentSettings.numSamples = 1e4;
GOTable = SingleEnrichment(geneScores,geneEntrezIDs,enrichmentSettings);
```
