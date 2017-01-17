function [summaryTable,allGOLabels,allGONames,allGOIDs,ix_runs] = PrepareSummaryTable(enrichmentTables,doReorder)
% Takes in enrichment tables from RunErmineJ
%-------------------------------------------------------------------------------

if nargin < 2
    doReorder = false;
end

% 1. What are the different GO IDs implicated
GOIDs = cellfun(@(x)x.GOID,enrichmentTables,'UniformOutput',0);
allGOIDs = unique(vertcat(GOIDs{:}));
numGOIDs = length(allGOIDs);
% map to names
allGONames = cell(numGOIDs,1);
for i = 1:numGOIDs
    isHere = cellfun(@(x)ismember(allGOIDs{i},x),GOIDs);
    isHere = find(isHere,1,'first');
    thisRow = strcmp(enrichmentTables{isHere}.GOID,allGOIDs{i});
    allGONames{i} = enrichmentTables{isHere}.GOName{thisRow};
end
allGOLabels = arrayfun(@(x)sprintf('%s (%s)',allGONames{x},allGOIDs{x}),...
                            1:numGOIDs,'UniformOutput',false);

% 2. Prepare output
numRuns = length(enrichmentTables);
summaryTable = nan(numGOIDs,numRuns);
for i = 1:numRuns
    if isempty(enrichmentTables{i}), continue; end
    [~,ia,ib] = intersect(allGOIDs,enrichmentTables{i}.GOID);
    summaryTable(ia,i) = enrichmentTables{i}.pVal(ib);
end

%-------------------------------------------------------------------------------
% Order GO categories by relevance:
if doReorder
    summaryTableSat = summaryTable;
    summaryTableSat(isnan(summaryTable)) = max(summaryTable(:));
    propSig = mean(summaryTableSat,2);
    [~,ix_GOcat] = sort(propSig,'descend');
    propSig = mean(summaryTableSat,1);
    [~,ix_runs] = sort(propSig,'descend');
    % Now reorder
    allGOLabels = allGOLabels(ix_GOcat);
    allGONames = allGONames(ix_GOcat);
    allGOIDs = allGOIDs(ix_GOcat);
    summaryTable = summaryTable(ix_GOcat,ix_runs);
else
    ix_runs = 1:numRuns;
end

end
