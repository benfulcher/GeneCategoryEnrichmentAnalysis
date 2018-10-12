function GOTable = SingleEnrichment(geneScores,geneEntrezIDs,dataSource,processFilter,sizeFilter,numSamples)
% Enrichment under a random-gene null
%
%---INPUTS:
% * `geneScores`, a numGenes-long column vector of values that quantifies something about each gene.
% * `geneEntrezIDs`, numGenes-long column vector labeling the entrez ID for each gene in geneScores.
% * `dataSource`, specifies the source of GO annotations, to be loaded using `GetFilteredGOData`. Options are `mouse-direct` (hierarchy and annotations taken directly from GO), `human-direct` (hierarchy and annotations taken directly from GO), `mouse-GEMMA` (processed hierarchy and annotations downloaded from GEMMA).
% * `processFilter`, what GO processes to consider. Default is `biological_process`.
% * `sizeFilter`, filter GO categories by size. Default is `[5,200]` (only consider categories with between 5 and 200 annotations).
% * `numSamples`, number of permutation iterations for null (`1e4` is the default, can ramp up to get better significance estimates for small p-values).
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs and set defaults:
if nargin < 3
    dataSource = 'mouse-direct';
end
if nargin < 4
    processFilter = 'biological_process';
end
if nargin < 5
    sizeFilter = [5,100];
end
if nargin < 6
    numSamples = 1e4;
end

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
% Compute the mean score for within each category:
categoryScores = nan(numGOCategories,1);
for j = 1:numGOCategories
    matchMe = ismember(geneEntrezIDs,GOTable.annotations{j});
    if sum(matchMe) == 0
        continue
    end
    categoryScores(j) = nanmean(geneScores(matchMe));
end

%-------------------------------------------------------------------------------
% Generate a null distribution (through permutation testing) for each category size
nullDistribution = PermuteForNullDistributions(geneScores,uniqueSizes,numSamples);
% numSizes x numSamples


keyboard

%-------------------------------------------------------------------------------
% Compute p-values (bigger scores are better)
pVals = nan(numGOCategories,1);
parfor i = 1:numGOCategories
    if ~isnan(categoryScores(i))
        nullForSize = nullDistribution(GOTable.size(i)==uniqueSizes,:);
        pVals(i) = mean(categoryScores(i) < nullForSize);
    end
end

%-------------------------------------------------------------------------------
% FDR correct:
pVals_corr = mafdr(pVals,'BHFDR','true');

%-------------------------------------------------------------------------------
% Update the GO table:
GOTable.pVal = pVals;
GOTable.pValCorr = pVals_corr;
GOTable.meanScore = categoryScores;

%-------------------------------------------------------------------------------
% Sort:
GOTable = sortrows(GOTable,{'pVal','pValCorr'},{'ascend','ascend'});

end
