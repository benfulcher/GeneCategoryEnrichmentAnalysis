function GOTablePhenotype = EnsembleEnrichment(params,phenotypeVector,whatNullModel)
% EnsembleEnrichment  Compute enrichment in different GO categories according to
%                       a given null model.
% Assumes that nulls have been precomputed using ComputeAllCategoryNulls.
%
%
%---INPUTS:
% params: a structure of relevant parameters
% phenotypeVector: a vector of the spatial phenotype map to be tested
% whatNullModel: the null model to use for enrichment
%
% params:
% params.whatCorr ('Spearman').
%
%---OUTPUT:
% GOTablePhenotype: a table with p-values estimated from the null ensemble

%-------------------------------------------------------------------------------
% Process Inputs and Set Defaults:
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultEnsembleParams(params);
end
if nargin < 2
    myPhenotype = 'degree';
end
if nargin < 3
    whatNullModel = 'randomMap'; % 'spatialLag'
end

% Pull out fields from the params structure:
if isfield(params,'whatCorr')
    whatCorr = params.whatCorr;
else
    whatCorr = 'Spearman';
end
if isfield(params,'aggregateHow')
    aggregateHow = params.aggregateHow;
else
    aggregateHow = 'mean';
end
if isfield(params,'numNullSamples')
    numNullSamples = params.numNullSamples;
else
    numNullSamples = 40000;
end
if isfield(params,'whatSpecies')
    whatSpecies = params.whatSpecies;
else
    whatSpecies = 'mouse';
end
if isfield(params,'structFilter')
    structFilter = params.structFilter;
else
    structFilter = 'all';
end

%-------------------------------------------------------------------------------
% Load null distributions: GOTableNull
%-------------------------------------------------------------------------------
% Check for precomputed null results (from running ComputeAllCategoryNulls):
fileNullEnsembleResults = sprintf('RandomNull_%u_%s-%s_%s_%s_%s.mat',numNullSamples,...
                    humanOrMouse,structFilter,whatNullModel,whatCorr,aggregateHow);
fileNullEnsembleResults = fullfile('EnsembleEnrichment','NullEnsembles',fileNullEnsembleResults);

preComputedNulls = load(fileNullEnsembleResults);
GOTableNull = preComputedNulls.GOTable;

if isfield(preComputedNulls,'params')
    error('UNUSABLE OLD FILE: %s\n',fileNullEnsembleResults);
else
    warning('***ASSUME THAT KEY PARAMETERS MATCH THE LOADED FILE: %s\n',fileNullEnsembleResults);
    params = preComputedData.params;
end
% (The categoryScores variable is the distribution of null samples for each GO category)

%-------------------------------------------------------------------------------
% Now compute scores for the real phenotype
%-------------------------------------------------------------------------------
GOTablePhenotype = ComputeAllCategoryNulls(params,1,myPhenotype,whatCorr,aggregateHow,false,false);

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

numSig = sum(GOTablePhenotype.pValZCorr < params.e.sigThresh);
fprintf(1,'%u significant categories at pZ_corr < %.2f\n',numSig,params.e.sigThresh);

[geneData,geneInfo,structInfo] = LoadMeG(params.g);
ListCategories(geneInfo,GOTablePhenotype,20,'pValZ');

end
