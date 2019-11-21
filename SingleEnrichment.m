function GOTable = SingleEnrichment(geneScores,geneEntrezIDs,params);
% Enrichment under a random-gene null
%
%---INPUTS:
% * geneScores, a numGenes-long column vector of values that quantifies
%               something about each gene.
% * geneEntrezIDs, numGenes-long column vector labeling the entrez ID for each
%                  gene in geneScores.
% * params, structure with the following fields:
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
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs and set defaults:
if nargin < 3
    params = struct();
    params.dataSource = 'mouse-direct';
    params.processFilter = 'biological_process';
    params.sizeFilter = [5,100];
    params.numSamples = 1e4;
end

% Convert out of structure:
dataSource = params.dataSource;
processFilter = params.processFilter;
sizeFilter = params.sizeFilter;
numSamples = params.numSamples;

%-------------------------------------------------------------------------------
if length(geneScores)~=length(geneEntrezIDs)
    error('Inconsistent numbers of genes scored (%u) and with entrez IDs provided (%u)',...
                        length(geneScores),length(geneEntrezIDs));
end
numGenes = length(geneScores);

% Retrieve GO annotations:
GOTable = GetFilteredGOData(dataSource,processFilter,sizeFilter,geneEntrezIDs);
numGOCategories = height(GOTable);
if all(GOTable.size==0)
    error('No annotations found in %s after applying filters',dataSource);
end

%-------------------------------------------------------------------------------
uniqueSizes = unique(GOTable.size);
numSizes = length(uniqueSizes);
fprintf(1,'Gene score resampling for %u iterations across %u category sizes (%u-%u)\n',...
                        numSamples,numSizes,min(uniqueSizes),max(uniqueSizes));

%-------------------------------------------------------------------------------
% Compute the mean score for each GO category:
categoryScores = nan(numGOCategories,1);
for i = 1:numGOCategories
    matchMe = ismember(geneEntrezIDs,GOTable.annotations{i});
    if sum(matchMe) > 0
        categoryScores(i) = nanmean(geneScores(matchMe));
    end
end

%-------------------------------------------------------------------------------
% Generate a null distribution (through permutation testing) for each category size
nullDistribution = PermuteForNullDistributions(geneScores,uniqueSizes,numSamples);
% nullDistribution: numSizes x numSamples

%-------------------------------------------------------------------------------
% Compute p-values (bigger scores are better)
pVal = nan(numGOCategories,1); % Discrete, count-based estimate from size-stratified permutation
pValZ = nan(numGOCategories,1); % Gaussian approximation
parfor i = 1:numGOCategories
    if ~isnan(categoryScores(i))
        nullForSize = nullDistribution(GOTable.size(i)==uniqueSizes,:);
        pVal(i) = mean(categoryScores(i) < nullForSize);
        pValZ(i) = 1 - normcdf(categoryScores(i),mean(nullForSize),std(nullForSize));
    end
end

%-------------------------------------------------------------------------------
% FDR 'correction' to both methods of estimating p-values:
pValCorr = mafdr(pVal,'BHFDR','true');
pValZCorr = mafdr(pValZ,'BHFDR','true');

%-------------------------------------------------------------------------------
% Update the GO table:
GOTable.meanScore = categoryScores;
GOTable.pValZ = pValZ; % Estimated p-value from Gaussian approximation to null
GOTable.pValPerm = pVal; % Permutation test
GOTable.pValZCorr = pValZCorr;
GOTable.pValPermCorr = pValCorr;

%-------------------------------------------------------------------------------
% Sort:
GOTable = sortrows(GOTable,{'pValPerm','pValZ'},{'ascend','ascend'});

end
