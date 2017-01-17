function [GOTable,geneEntrezAnnotations] = GetFilteredGOData(whatFilter,sizeFilter,ourEntrez)

if nargin < 1
    whatFilter = 'biological_process';
end
if nargin < 2
    sizeFilter = [5,200];
end
if nargin < 3
    ourEntrez = [];
end

% Get GO annotation data (processed):
load('GOAnnotation.mat','allGOCategories','geneEntrezAnnotations');

% Get GO ontology details
GOTable = GetGOTerms(whatFilter);

%-------------------------------------------------------------------------------
% Filter
%-------------------------------------------------------------------------------
% Filter by ontology details:
[~,ia,ib] = intersect(GOTable.GOID,allGOCategories);
fprintf(1,'Filtering to %u annotated GO categories related to %s\n',length(ia),whatFilter);
GOTable = GOTable(ia,:);
allGOCategories = allGOCategories(ib);
geneEntrezAnnotations = geneEntrezAnnotations(ib);

% Filter by category size:
numGOCategories = length(allGOCategories);
if isempty(ourEntrez)
    fprintf(1,'Filtering on actual annotated size of GO category\n');
    sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
else
    fprintf(1,'Filtering on size of GO category after matching to our genes\n');
    % Make sure all annotations match those in our set of entrez_ids
    geneEntrezAnnotations = cellfun(@(x)x(ismember(x,ourEntrez)),geneEntrezAnnotations,'UniformOutput',false);
    sizeGOCategories = cellfun(@length,geneEntrezAnnotations);
end
isGoodSize = (sizeGOCategories >= sizeFilter(1)) & (sizeGOCategories <= sizeFilter(2));
allGOCategories = allGOCategories(isGoodSize);
geneEntrezAnnotations = geneEntrezAnnotations(isGoodSize);
GOTable = GOTable(isGoodSize,:);
sizeGOCategories = sizeGOCategories(isGoodSize);
numGOCategories = length(allGOCategories);
fprintf(1,'Filtered to %u categories with between %u and %u annotations\n',...
                numGOCategories,sizeFilter(1),sizeFilter(2));

end
