% abagen lh-mirror function was set on -> different more stable gene
% expression data are used here.
% New T values obtained from A3a_lm.R using multiple linear regression as
% method to correct for sex. The Tmaps are corrected for scanner (using combat)
% and sex (using MLR in R).


%% Setup your gene category data
clear;clc;

% setup all parameters and directories
startup; 

%% Ensemble enrichment relative to an ensemble of null phenotypes as suggested by Fulcher 2021
%---------------------------------------------------------------------------
% 1. Prepare null models
%---------------------------------------------------------------------------
%read in entrez_id.csv and entrez_id_expression.csv generated with conda
%env abagen > A1_gene_preprocessing.ipynb as table
% define geneDataStruct needed for ComputeAllCategoryNulls
geneDataStruct.expressionMatrix = readmatrix('entrez_id_expression_lh-mirror.csv');
geneDataStruct.entrezIDs = readmatrix('entrez_ids_lh-mirror.csv');

%---------------------------------------------------------------------------
%% 2. Compute category-score null distribution for every gene category
%---------------------------------------------------------------------------
% null distributions relative for category scores relative to the specified
% phenotype ensemble are computed...
%enrichmentParams = GiveMeDefaultEnrichmentParams();
%GOTableNull = ComputeAllCategoryNulls(geneDataStruct,enrichmentParams,[],true,true);

Null = load('Nulls_lhmirror_genesize10-1000_randomMap_5000.mat');
enrichmentParams=Null.enrichmentParams;

%---------------------------------------------------------------------------
%% 3. Compute enrichment of a given phenotype
%---------------------------------------------------------------------------
% assess significance of the scores obtained for a given phenotype relative
% to these nulls...
sigThresh = 0.05;

phenotypeVectorThickness = readmatrix('thickness_tvals_multipleRegression_combatCorr.csv')';
phenotypeVectorThickness = phenotypeVectorThickness(2,:)';

GOTablePhenotypeThickness = EnsembleEnrichment(geneDataStruct,enrichmentParams.fileNameOut,phenotypeVectorThickness);
numSigThick = sum(GOTablePhenotypeThickness.pValPermCorr < sigThresh);


phenotypeVectorGyrification = readmatrix('gyrification_tvals_multipleRegression_combatCorr.csv')';
phenotypeVectorGyrification = phenotypeVectorGyrification(2,:)';

GOTablePhenotypeGyrification = EnsembleEnrichment(geneDataStruct,enrichmentParams.fileNameOut,phenotypeVectorGyrification);
numSigGyr = sum(GOTablePhenotypeGyrification.pValPermCorr < sigThresh);


