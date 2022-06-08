[Gene Ontology (GO)](http://geneontology.org/docs/download-ontology/) no longer provides their ontology as MySQL database dumps, meaning that to use new GO releases with `GeneCategoryEnrichmentAnalysis`, ontologies must be converted from raw .obo format into a database compatible with the SQL commands used in this package.

The code in this directory can be used to convert a `go-basic.obo` file to an sqlite .db file, which is compatible with the sqlite-adapted [DataProcessing/](https://github.com/alyssadai/GeneCategoryEnrichmentAnalysis/tree/sqlite/DataProcessing) scripts in the current repository.

These scripts are dependent on .txt versions of the `term` and `term2term` tables from the GO MySQL database schema, which can be obtained for your `go-basic.obo` file of interest using [OboToTerm](https://github.com/sgrote/OboToTerm).

### Software requirements

- sqlite3
- Python

### Running the code

1. Run `tsv_to_csv.py` to convert `term.txt` and `term2term.txt` (obtained using OboToTerm) to csvs for easy importing into tables in an sqlite db.
2. Follow the commands in `create_sqlitedb_from_obo.sql` to convert `term.csv` and `term2term.csv` into an sqlite db that can then be used to generate gene-category annotations as outlined in [this wiki](https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis/wiki/Defining-gene-to-category-annotations).
   * Note: These commands were adapted from an sqlite3-compatible dump (using [mysql2sqlite](https://github.com/dumblob/mysql2sqlite)) of `GODaily_2021-01-25.sql`, which is the MySQL download of the GO data used in the original Fulcher et al., 2021 paper and provided in the [partner Zenodo data repo](https://zenodo.org/record/4460714). `create_sqlitedb_from_obo.sql` also outputs a .sql dump of the generated sqlite db which can be compared against `GODaily_2021-01-25.sql` for major schema discrepancies as an additional sanity check.
