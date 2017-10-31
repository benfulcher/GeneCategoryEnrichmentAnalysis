function [GOTable,geneEntrezAnnotations] = SingleEnrichment(geneScores,geneEntrezIDs,dataSource,processFilter,sizeFilter,numIters)
%-------------------------------------------------------------------------------
% Do an ermineJ style analysis for a given set of entrezIDs and scores
%-------------------------------------------------------------------------------

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
    numIters = 10000;
end
%-------------------------------------------------------------------------------
numGenes = length(geneScores);

% Retrieve GO annotations:
GOTable = GetFilteredGOData(dataSource,processFilter,sizeFilter,geneEntrezIDs);
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Compute the mean score for within each category:
categoryScores = nan(numGOCategories,1);
for j = 1:numGOCategories
    matchMe = ismember(geneEntrezIDs,GOTable.annotations{j});
    if sum(matchMe) <= 1
        continue
    end
    categoryScores(j) = nanmean(geneScores(matchMe));
end

%-------------------------------------------------------------------------------
% Generate a null distribution
uniqueSizes = unique(sizeGOCategories);
numSizes = length(uniqueSizes);
nullDistribution = zeros(numSizes,numIters);
fprintf(1,'Gene score reasmpling for %u iterations across %u category sizes (%u-%u)\n',...
                        numIters,numSizes,min(uniqueSizes),max(uniqueSizes));
parfor j = 1:numIters
    rp = randperm(numGenes); % takes a millisecond to compute this (put outside inner loop)
    for i = 1:numSizes
        nullDistribution(i,j) = nanmean(geneScores(rp(1:uniqueSizes(i))));
    end
end

%-------------------------------------------------------------------------------
% Compute p-values
pVals = zeros(numGOCategories,1);
parfor i = 1:numGOCategories
    % Bigger is better:
    pVals(i) = mean(categoryScores(i) < nullDistribution(uniqueSizes==sizeGOCategories(i),:));
end

%-------------------------------------------------------------------------------
% FDR correct:
pVals_corr = mafdr(pVals,'BHFDR','true');

%-------------------------------------------------------------------------------
% Update the GO table:
GOTable.pVal = pVals;
GOTable.pVal_corr = pVals_corr;
GOTable.meanScore = categoryScores;

%-------------------------------------------------------------------------------
% Sort:
[GOTable,ix] = sortrows(GOTable,{'pVal','pVal_corr'},{'ascend','ascend'});

end
