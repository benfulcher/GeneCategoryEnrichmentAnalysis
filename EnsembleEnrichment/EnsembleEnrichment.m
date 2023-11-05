function GOTablePhenotype = EnsembleEnrichment(geneDataStruct,fileNullEnsembleResults,phenotypeVector)
% EnsembleEnrichment  Compute enrichment across GO categories for a given null model.
%
% Assumes that nulls have been precomputed using ComputeAllCategoryNulls.
%
%---INPUTS:
% - geneDataStruct: structure containing expression data from which to evaluate
%                   the phenotypeVector (cf. ComputeAllCategoryNulls)
% - fileNullEnsembleResults: file where the precomputed nulls are stored
% - phenotypeVector: a vector of the spatial phenotype map to be tested
%
% The parameters are taken from the null ensemble enrichment file (params)
%
%---OUTPUT:
% - GOTablePhenotype: a table with p-values estimated from the null ensemble

%-------------------------------------------------------------------------------
% Process inputs and set defaults:
%-------------------------------------------------------------------------------
if nargin < 1
    error('You must specify a file containing the gene data');
end
if nargin < 2
    error('You must specify a file containing the precomputed ensemble nulls');
end
if nargin < 3
    error('You must provide a phenotype vector');
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Load null distributions into GOTableNull
%-------------------------------------------------------------------------------
preComputedNulls = load(fileNullEnsembleResults,'GOTable','enrichmentParams');
GOTableNull = preComputedNulls.GOTable;

%-------------------------------------------------------------------------------
% Now compute scores for the input phenotype using the same parameter settings
% as for the null distribution computation:
%-------------------------------------------------------------------------------
enrichmentParams = preComputedNulls.enrichmentParams;
enrichmentParams.whatEnsemble = 'customSpecified'; % allows us to feed in our custom phenotype
GOTablePhenotype = ComputeAllCategoryNulls(geneDataStruct,enrichmentParams,phenotypeVector,false,false);

% Check that we have the same GO category IDs in both cases:
numCategories = height(GOTablePhenotype);
if ~(height(GOTableNull)==numCategories) && ~all(GOTableNull.GOID==GOTablePhenotype.GOID)
    error('Error matching GO Categories to precomputed null data...');
end

%-------------------------------------------------------------------------------
% Estimate p-values:
%-------------------------------------------------------------------------------
if isfield(enrichmentParams,'whatTail')
    whatTail = enrichmentParams.whatTail;
else
    whatTail = 'right';
    fprintf(1,'Right-tailed test by default: larger values of %s %s correlation are interesting.\n',...
                    enrichmentParams.aggregateHow,enrichmentParams.whatCorr);
end

GOTablePhenotype = EstimatePVals(GOTableNull.categoryScores,...
                    [GOTablePhenotype.categoryScores{:}],whatTail,GOTablePhenotype);

% Sort by Gaussian-approximation p-values:
%GOTablePhenotype = sortrows(GOTablePhenotype,'pValZ','ascend');
GOTablePhenotype = sortrows(GOTablePhenotype,'pValPermCorr','ascend');

% Give a basic output about significance using pValZCorr:
sigThresh = 0.05;
%numSig = sum(GOTablePhenotype.pValZCorr < sigThresh);
numSig = sum(GOTablePhenotype.pValPermCorr < sigThresh);
%fprintf(1,'%u significant categories at p_Z_corr < %.2f\n',numSig,sigThresh);
fprintf(1,'%u significant categories at p_Perm_corr < %.2f\n',numSig,sigThresh);
end
