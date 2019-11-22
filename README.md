# Gene Set Enrichment Analysis by Gene Score Resampling

This is a Matlab repository for performing Gene Set Enrichment Analysis (GSEA) on Gene Ontology (GO) Terms.
Given a set of scores assigned to genes, we use the permutation-based method of Gene Score Resampling (GSR), following the [*ermineJ*](https://erminej.msl.ubc.ca/) software package.
The package supports conventional enrichment, which generates a null distribution for some summary statistic across the scores assigned to genes annotated to a given GO Term by randomizing the assignment of genes to GO categories.
These null distributions depend on the size of the GO Category.

The package is currently set up to perform enrichment on GO Biological Processes.

There are two steps in performing an analysis:
1. Initialize the package with the latest GO hierarchy and gene annotations (or download processed results from [figshare](https://figshare.com/s/71fe1d9b2386ec05f421)).
2. Run an enrichment analysis.

___Note:___
The [figshare repository](https://figshare.com/s/71fe1d9b2386ec05f421) uses GO annotation files from the 2019-04-17 release.

## Analysis
### Gene Score Resampling
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
enrichmentSettings = struct();
enrichmentSettings.dataSource = 'mouse-direct';
enrichmentSettings.processFilter = 'biological_process';
enrichmentSettings.sizeFilter = [5,200];
enrichmentSettings.numSamples = 1e4;
GOTable = SingleEnrichment(geneScores,geneEntrezIDs,enrichmentSettings);
```

## Initialization

We first describe how to retrieve and process data from GO in a form that facilitates GSEA.
Note that this step can be skipped by downloading pre-processed results from [figshare](https://figshare.com/s/71fe1d9b2386ec05f421).

This involves the following steps:
1. Retrieve and process the GO hierarchy data.
2. Retrieve and process the annotations of genes to GO Terms.
3. Iteratively propagate gene-to-Term annotations from child to parent up the GO hierarchy.

### Processing the GO-Term Hierarchy
#### Downloading the data
There are a number of routes to [downloading the GO Term hierarchy](http://geneontology.org/page/download-ontology).
We used the **termdb** [mySQL database dump](http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz), and linked to this database from Matlab using a mySQL java connector.
Code for achieving this (e.g., in the [Matlab_mySQL repository](https://github.com/benfulcher/Matlab_mySQL)) is a dependency for this package.

___Note___: the data is also provided in raw form as `go-basic.obo` (the `basic` file ensures that annotations can be propagated), and you can also download the data as a [database](ftp://ftp.geneontology.org/go/www/GO.downloads.database.shtml).

#### Reading and Processing:
1. Set up downloaded `termdb` mySQL database, and put connection details in `ConnectMeDatabase`.
2. Retrieve Biological Process GO Terms, and save the filtered set of terms to a .mat file:

```matlab
GOTerms = GetGOTerms('biological_process',true);
```
Saves out to `ProcessedData/GOTerms_BP.mat`.

### Processing GO Term Annotations
Now that we have the GO Terms in Matlab format, we next need data on which genes are annotated to which GO Terms.

#### Downloading raw GO annotation data

Annotation files should be downloaded directly from the [GO website](http://current.geneontology.org/products/pages/downloads.html).
* For _Mus musculus_, the annotation file is `mgi.gaf`.
* For _Homo sapiens_, the annotation file is `goa_human.gaf`.

The appropriate annotation file(s) should be placed in the `RawData` directory.

#### Processing data from the annotation file

Each line in the [annotation file](http://geneontology.org/page/go-annotation-file-formats) represents an association between a gene product and a GO term with a certain evidence code, and the reference to support the association.
The `ReadDirectAnnotationFile` function reads in all of this raw data, and processes it into a Matlab table, with a row for each GO Category, including information about the category and the genes that are annotated to it.

Before this can be run, it requires a mapping from MGI gene identifiers to NCBI Entrez gene identifiers.
In mouse, this is achieved by taking data from [_MouseMine_](http://www.mousemine.org).
```bash
python3 MGI_NCBI_downloadall.py
```
This saves the required gene identifier mapping to `ALL_MGI_ID_NCBI.csv`.
In the case of human data, we mapped onto gene symbols from processed gene-expression data from the [Allen Human Brain Atlas](https://human.brain-map.org/).

```matlab
ReadDirectAnnotationFile('mouse')
```

Saves processed data as `GOAnnotationDirect-mouse.mat` (or `GOAnnotationDirect-human.mat`), in the `ProcessedData` directory.

___Note (NOT RECOMMENDED)___: Annotations processed from [GEMMA](https://gemma.msl.ubc.ca/annots/) can alternatively be read using `ReadGEMMAAnnotationFile`.

#### Propagate annotations up through the hierarchy
Annotations are made at the lowest level of the GO term hierarchy.
Annotations at a lower level of the hierarchy apply to all parent terms.
For performing enrichment, we therefore need to iteratively propagate direct annotations up the hierarchy, using `is_a` (child-parent) relationships from the `term2term` table from the GO Term database.

For mouse biological processes, this is achieved using:
```matlab
propagateHierarchy('mouse','biological_process');
```
The code takes processed data (e.g., `GOAnnotationDirect-mouse.mat`) and saves propagated output as `GOAnnotationDirect-mouse-biological_process-Prop.mat`.
This information can then be read in for enrichment or other analyses.

Note that this script requires access to a mySQL database containing the hierarchical GO term information.
