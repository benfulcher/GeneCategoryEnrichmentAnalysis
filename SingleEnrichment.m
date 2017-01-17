function GOTable = SingleEnrichment(geneScores,geneEntrezIDs,processFilter,sizeFilter,numIters)
%-------------------------------------------------------------------------------
% Do an ermineJ style analysis for a given set of entrezIDs and scores
%-------------------------------------------------------------------------------

if nargin < 3
    processFilter = 'biological_process';
end
if nargin < 4
    sizeFilter = [5,200];
end
if nargin < 5
    numIters = 10000;
end
%-------------------------------------------------------------------------------
numGenes = length(geneScores);

% Retrieve GO annotations:
[GOTable,geneEntrezAnnotations] = GetFilteredGOData(processFilter,sizeFilter,geneEntrezIDs);
sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
% Add to table:
GOTable.size = sizeGOCategories;
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Compute the mean score for within each category:
categoryScores = nan(numGOCategories,1);
for j = 1:numGOCategories
    matchMe = ismember(geneEntrezIDs,geneEntrezAnnotations{j});
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
for j = 1:numIters
    rp = randperm(numGenes); % takes a millisecond to compute this (put outside inner loop)
    for i = 1:numSizes
        nullDistribution(i,j) = nanmean(geneScores(rp(1:uniqueSizes(i))));
    end
end

%-------------------------------------------------------------------------------
% Compute p-values
pVals = zeros(numGOCategories,1);
for i = 1:numGOCategories
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

end
