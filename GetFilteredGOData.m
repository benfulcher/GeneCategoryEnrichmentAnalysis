function [GOTable,geneEntrezAnnotations] = GetFilteredGOData(whatSource,whatFilter,sizeFilter,ourEntrez)

if nargin < 1
    whatSource = 'direct'; % Direct annotations from GO
    % whatSource = 'GEMMA'; % Annotations derived from GEMMA
end
if nargin < 2
    whatFilter = 'biological_process';
end
if nargin < 3
    sizeFilter = [5,200];
end
if nargin < 4
    ourEntrez = [];
end

%-------------------------------------------------------------------------------
% Load processed GO annotation data (i.e., direct annotations propagated up the hierarchy):
% cf. propagateHierarchy to map files generated from ReadDirectAnnotationFile or ReadGEMMAAnnotationFile
switch whatSource
case 'direct'
    fileNameLoad = 'GOAnnotationDirectProp.mat';
case 'GEMMA'
    fileNameLoad = 'GOAnnotationGEMMAProp.mat';
otherwise
    error('Unknown annotation source: ''%s''',whatSource);
end
load('GOAnnotationDirect.mat','allGOCategories','geneEntrezAnnotations');

% Get GO ontology details
if strcmp(whatFilter,'biological_process') && exist('GOTerms_BP.mat','file')
    load('GOTerms_BP.mat','GOTable');
else
    % Retrieve from mySQL database:
    GOTable = GetGOTerms(whatFilter);
end

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
GOTable.size = sizeGOCategories;
fprintf(1,'GO categories have between %u and %u annotations\n',min(sizeGOCategories),max(sizeGOCategories));
isGoodSize = (sizeGOCategories >= sizeFilter(1)) & (sizeGOCategories <= sizeFilter(2));
allGOCategories = allGOCategories(isGoodSize);
geneEntrezAnnotations = geneEntrezAnnotations(isGoodSize);
GOTable = GOTable(isGoodSize,:);
sizeGOCategories = sizeGOCategories(isGoodSize);
numGOCategories = length(allGOCategories);
fprintf(1,'Filtered to %u categories with between %u and %u annotations\n',...
                numGOCategories,sizeFilter(1),sizeFilter(2));

end
