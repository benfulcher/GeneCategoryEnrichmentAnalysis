function GOTable = SingleEnrichment(geneScores,geneEntrezIDs,params)
% SingleEnrichment   Perform gene category enrichment under a random-gene null
%
% This follows the permutation-based method of Gene Score Resampling, as
% implemented in ermineJ (https://erminej.msl.ubc.ca/), which assesses the
% significance of the scores assigned to genes in a given category of genes
% relative to an ensemble of categories containing the same number of random
% genes.
%
%---INPUTS:
% * geneScores, a numGenes-long column vector of values that quantifies
%               something about each gene.
% * geneEntrezIDs, numGenes-long column vector labeling the entrez ID for each
%                  gene in geneScores.
% * params, a structure with the following fields:
%     - dataSource, specifies the source of GO annotations, to be loaded using
%                   `GetFilteredGOData`. Options are `mouse-direct` (hierarchy
%                   and annotations taken directly from GO), `human-direct`
%                   (hierarchy and annotations taken directly from GO),
%                   `mouse-GEMMA` (processed hierarchy and annotations
%                   downloaded from GEMMA).
%     - processFilter, what GO processes to consider. Default: `biological_process`.
%     - sizeFilter, filter GO categories by size. Default: `[5,200]`
%                   (i.e., consider categories with between 5 and 200 annotations).
%     - numSamples, number of permutation iterations for null. Default: 1e4.
%                   Can ramp up to get better significance estimates for small p-values).
%     - whatTail, whether higher ('right') or lower ('left') scores are 'better'
%                   (when computing p-values from null).
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs and set defaults:
if nargin < 3
    params = GiveMeDefaultEnrichmentParams();
end

% Convert out of structure:
dataSource = params.dataSource;
processFilter = params.processFilter;
sizeFilter = params.sizeFilter;
numSamples = params.numNullSamples;
if isfield(params,'whatTail');
    whatTail = params.whatTail;
else
    whatTail = 'right';
end
if isfield(params,'whatTail');
    aggregateHow = params.aggregateHow;
else
    aggregateHow = 'mean';
end


%-------------------------------------------------------------------------------
numGenes = length(geneScores);
if numGenes ~= length(geneEntrezIDs)
    error('Inconsistent numbers of genes scored (%u) and with entrez IDs provided (%u)',...
            numGenes,length(geneEntrezIDs));
end

% Retrieve GO annotations:
GOTable = GetFilteredGOData(dataSource,processFilter,sizeFilter,geneEntrezIDs);
numGOCategories = height(GOTable);
if all(GOTable.size==0)
    error('No annotations found in %s after applying filters',dataSource);
end

%-------------------------------------------------------------------------------
uniqueSizes = unique(GOTable.size);
numSizes = length(uniqueSizes);
fprintf(1,'Gene-score resampling for %u iterations across %u category sizes (%u-%u)\n',...
                numSamples,numSizes,min(uniqueSizes),max(uniqueSizes));

%-------------------------------------------------------------------------------
% Compute the mean score for each GO category:
categoryScores = nan(numGOCategories,1);
for i = 1:numGOCategories
    matchMe = ismember(geneEntrezIDs,GOTable.annotations{i});
    if sum(matchMe) > 0
        geneScoresHere = geneScores(matchMe);
        categoryScores(i) = AggregateScores(geneScoresHere,aggregateHow);
    end
end

%-------------------------------------------------------------------------------
% Generate a null distribution (through permutation testing) for each category size
nullDistribution = PermuteForNullDistributions(geneScores,uniqueSizes,numSamples);
% nullDistribution: numSizes x numSamples

%-------------------------------------------------------------------------------
% Compute p-values
pValPerm = nan(numGOCategories,1); % Discrete, count-based estimate from size-stratified permutation
pValZ = nan(numGOCategories,1); % Gaussian approximation
GOCategorySizes = GOTable.size; % Speeds up parfor loop
parfor i = 1:numGOCategories
    if ~isnan(categoryScores(i))
        nullForSize = nullDistribution(GOCategorySizes(i)==uniqueSizes,:);
        switch whatTail
        case 'right'
            pValPerm(i) = mean(nullForSize >= categoryScores(i));
            pValZ(i) = 1 - normcdf(categoryScores(i),mean(nullForSize),std(nullForSize));
        case 'left'
            pValPerm(i) = mean(nullForSize <= categoryScores(i));
            pValZ(i) = normcdf(categoryScores(i),mean(nullForSize),std(nullForSize));
        otherwise
            error('Unknown tail setting, ''%s''',whatTail)
        end
    end
end

%-------------------------------------------------------------------------------
% FDR-correction to p-values:
pValPermCorr = mafdr(pValPerm,'BHFDR','true');
pValZCorr = mafdr(pValZ,'BHFDR','true');

%-------------------------------------------------------------------------------
% Update the GO table:
GOTable.meanScore = categoryScores;
GOTable.pValZ = pValZ; % Estimated p-value from Gaussian approximation to null
GOTable.pValPerm = pValPerm; % Permutation test
GOTable.pValZCorr = pValZCorr;
GOTable.pValPermCorr = pValPermCorr;

%-------------------------------------------------------------------------------
% Sort:
GOTable = sortrows(GOTable,{'pValPerm','pValZ'},{'ascend','ascend'});

end
