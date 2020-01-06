function GOTablePhenotype = EnsembleEnrichment(fileNullEnsembleResults,phenotypeVector)
% EnsembleEnrichment  Compute enrichment in different GO categories according to
%                       a given null model.
%
% Assumes that nulls have been precomputed using ComputeAllCategoryNulls.
%
%---INPUTS:
% phenotypeVector: a vector of the spatial phenotype map to be tested
%
% The params are taken from the results of the null ensemble enrichment file
%
%---OUTPUT:
% GOTablePhenotype: a table with p-values estimated from the null ensemble

%-------------------------------------------------------------------------------
% Process inputs and set defaults:
%-------------------------------------------------------------------------------
if nargin < 1
    error('You must specify a file containing the precomputed ensemble nulls');
end
if nargin < 2
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
GOTablePhenotype = ComputeAllCategoryNulls(enrichmentParams,phenotypeVector,false,false);


% Check that we have the same GO category IDs in both cases:
if ~(height(GOTableNull)==height(GOTablePhenotype)) && ~all(GOTableNull.GOID==GOTablePhenotype.GOID)
    error('Error matching GO Categories to precomputed null data...');
end
numCategories = height(GOTablePhenotype);

%-------------------------------------------------------------------------------
% Estimate p-values:
%-------------------------------------------------------------------------------
GOTablePhenotype = EstimatePVals(GOTableNull.categoryScores,...
                        [GOTablePhenotype.categoryScores{:}],'right',GOTablePhenotype);
GOTablePhenotype = sortrows(GOTablePhenotype,'pValZ','ascend');

% Give a basic output about significance using pValZCorr:
numSig = sum(GOTablePhenotype.pValZCorr < enrichmentParams.sigThresh);
fprintf(1,'%u significant categories at pZ_corr < %.2f\n',numSig,params.e.sigThresh);

end
