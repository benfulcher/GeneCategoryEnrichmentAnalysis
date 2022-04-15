/* 
Alyssa Dai 2022 
Description: Coerces recent GO releases in go-basic.obo form into a sqlite db usable with https://github.com/alyssadai/GeneCategoryEnrichmentAnalysis 
(adapted https://github.com/benfulcher/GeneCategoryEnrichmentAnalysis).

Requirements: term.csv, term2term.csv files for .obo file of interest, converted from tsvs generated using https://github.com/sgrote/OboToTerm
Usage: 
1. Run `sqlite3 GO_YYYY-MM-DD.db` to create a new database for your ontology, then run the commands below in sqlite
*/

DROP TABLE IF EXISTS term;
CREATE TABLE `term` (
   `id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
   `name` varchar(255) NOT NULL DEFAULT '',
   `term_type` varchar(55) NOT NULL,
   `acc` varchar(255) NOT NULL,
   `is_obsolete` integer NOT NULL DEFAULT '0',
   `is_root` integer NOT NULL DEFAULT '0',
   `is_relation` integer NOT NULL DEFAULT '0',
   UNIQUE (`acc`),
   UNIQUE (`id`)
   );
.mode csv
.import term.csv term
CREATE TABLE `term2term` (
   `id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
   `relationship_type_id` integer NOT NULL,
   `term1_id` integer NOT NULL,
   `term2_id` integer NOT NULL,
   `complete` integer NOT NULL DEFAULT '0',
   UNIQUE (`term1_id`,`term2_id`,`relationship_type_id`)
   );
.mode csv
.import term2term.csv term2term
CREATE INDEX "idx_term_t1" ON "term" (`name`);
CREATE INDEX "idx_term_t2" ON "term" (`term_type`);
CREATE INDEX "idx_term_t3" ON "term" (`acc`);
CREATE INDEX "idx_term_t4" ON "term" (`id`,`acc`);
CREATE INDEX "idx_term_t5" ON "term" (`id`,`name`);
CREATE INDEX "idx_term_t6" ON "term" (`id`,`term_type`);
CREATE INDEX "idx_term_t7" ON "term" (`id`,`acc`,`name`,`term_type`);
CREATE INDEX "idx_term2term_tt1" ON "term2term" (`term1_id`);
CREATE INDEX "idx_term2term_tt2" ON "term2term" (`term2_id`);
CREATE INDEX "idx_term2term_tt3" ON "term2term" (`term1_id`,`term2_id`);
CREATE INDEX "idx_term2term_tt4" ON "term2term" (`relationship_type_id`);
.output GO_YYYY-MM-DD.sql /*dump newly created sqlite db back into a SQL file, for comparison with .sql file used in Fulcher 2021 paper*/
.dump
.exit
