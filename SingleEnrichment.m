function GOTable = SingleEnrichment(geneScores,geneEntrezIDs,params);
% Enrichment under a random-gene null
%
%---INPUTS:
% * `geneScores`, a numGenes-long column vector of values that quantifies something about each gene.
% * `geneEntrezIDs`, numGenes-long column vector labeling the entrez ID for each gene in geneScores.
% * `params` structure with the following fields:
%     - `dataSource`, specifies the source of GO annotations, to be loaded using `GetFilteredGOData`. Options are `mouse-direct` (hierarchy and annotations taken directly from GO), `human-direct` (hierarchy and annotations taken directly from GO), `mouse-GEMMA` (processed hierarchy and annotations downloaded from GEMMA).
%     - `processFilter`, what GO processes to consider. Default is `biological_process`.
%     - `sizeFilter`, filter GO categories by size. Default is `[5,200]` (only consider categories with between 5 and 200 annotations).
%     - `numSamples`, number of permutation iterations for null (`1e4` is the default, can ramp up to get better significance estimates for small p-values).
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
% numSizes x numSamples

%-------------------------------------------------------------------------------
% Compute p-values (bigger scores are better)
pVal = nan(numGOCategories,1);
parfor i = 1:numGOCategories
    if ~isnan(categoryScores(i))
        nullForSize = nullDistribution(GOTable.size(i)==uniqueSizes,:);
        pVal(i) = mean(categoryScores(i) < nullForSize);
    end
end

%-------------------------------------------------------------------------------
% FDR correct:
pValCorr = mafdr(pVal,'BHFDR','true');

%-------------------------------------------------------------------------------
% Update the GO table:
GOTable.meanScore = categoryScores;
GOTable.pVal = pVal;
GOTable.pValCorr = pValCorr;

%-------------------------------------------------------------------------------
% Sort:
GOTable = sortrows(GOTable,{'pVal','pValCorr'},{'ascend','ascend'});

end
