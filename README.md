# Gene Enrichment

This is a Matlab repository for performing gene ontology enrichment given a set of gene scores.

Some information and usage examples are below:

## Data Processing

### Retrieving GO Terms

Assumes you have an SQL database setup, with details entered into `GetGOTerms`

Get biological process GO terms:
```matlab
GOTerms = GetGOTerms('biological_process');
```

### Generating formatted files for Matlab, processed from raw GO annotations

#### Downloading GO term hierarchy information

[Downloading the GO Term hierarchy](http://geneontology.org/page/download-ontology) has a number of routes.
The data is provided in raw form as `go-basic.obo` (`basic` file ensures that annotations can be propagated), but you can also download it as a [database](ftp://ftp.geneontology.org/go/www/GO.downloads.database.shtml).
We have used the **termdb** [mySQL database dump](http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz).


#### Downloading raw GO annotation data

You can download annotation files direct from the [GO website](http://current.geneontology.org/products/pages/downloads.html).
For _Mus musculus_, this yields: `mgi.gaf`.
For human, it is: `goa_human.gaf`.

#### Read in from the GO annotation file:
```matlab
ReadDirectAnnotationFile
```
Alternative is to use `ReadGEMMAAnnotationFile` to read in processed annotations from GEMMA.

Saves processed data in the form `GOAnnotationDirect-mouse.mat`, in the `ProcessedData` directory.

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
