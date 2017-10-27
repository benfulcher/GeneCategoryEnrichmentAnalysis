function [GOTable,geneEntrezAnnotationsFull] = GetFilteredGOData(whatSource,whatFilter,sizeFilter,ourEntrez)

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
    fileNameLoad = sprintf('GOAnnotationDirect-%s-Prop.mat',whatFilter);
case 'GEMMA'
    fileNameLoad = 'GOAnnotationGEMMAProp.mat';
otherwise
    error('Unknown annotation source: ''%s''',whatSource);
end
load(fileNameLoad,'GOTerms','geneEntrezAnnotationsFull');
fprintf(1,'Loaded annotations from %s\n',fileNameLoad);

% Get GO ontology details
if strcmp(whatFilter,'biological_process') && exist('GOTerms_BP.mat','file')
    GOTermsFile = 'GOTerms_BP.mat';
    load(GOTermsFile,'GOTable');
    fprintf(1,'Loaded biological process GO terms from %s\n',GOTermsFile);
else
    % Retrieve from mySQL database:
    GOTable = GetGOTerms(whatFilter);
end

%-------------------------------------------------------------------------------
% Filter
%-------------------------------------------------------------------------------
% Filter by ontology details:
[~,ia,ib] = intersect(GOTable.GOID,GOTerms.GOID);
fprintf(1,'Filtering to %u annotated GO categories related to %s\n',length(ia),whatFilter);
GOTable = GOTable(ia,:);
GOTerms = GOTerms(ib,:);
geneEntrezAnnotationsFull = geneEntrezAnnotationsFull(ib);

% Filter by category size:
numGOCategories = height(GOTerms);
if isempty(ourEntrez)
    fprintf(1,'Filtering on actual annotated size of GO category\n');
    sizeGOCategories = cellfun(@length,geneEntrezAnnotationsFull);
    fprintf(1,'%u GO categories have no annotations :-/\n',...
                    sum(sizeGOCategories==0));
else
    fprintf(1,'Filtering on size of GO category after matching to our %u genes\n',length(ourEntrez));
    % Make sure all annotations match those in our set of entrez_ids
    geneEntrezAnnotationsFull = cellfun(@(x)x(ismember(x,ourEntrez)),geneEntrezAnnotationsFull,...
                                        'UniformOutput',false);
    sizeGOCategories = cellfun(@length,geneEntrezAnnotationsFull);
    fprintf(1,'%u GO categories have no annotations matching our %u genes\n',...
                    sum(sizeGOCategories==0),length(ourEntrez));
end
GOTable.size = sizeGOCategories;
fprintf(1,'GO categories have between %u and %u annotations\n',min(sizeGOCategories),max(sizeGOCategories));
isGoodSize = (sizeGOCategories >= sizeFilter(1)) & (sizeGOCategories <= sizeFilter(2));
GOTable = GOTable(isGoodSize,:);
GOTerms = GOTerms(isGoodSize,:);
geneEntrezAnnotationsFull = geneEntrezAnnotationsFull(isGoodSize);
sizeGOCategories = sizeGOCategories(isGoodSize);

%-------------------------------------------------------------------------------
numGOCategories = length(geneEntrezAnnotationsFull);
fprintf(1,'Filtered to %u categories with between %u and %u annotations\n',...
                numGOCategories,sizeFilter(1),sizeFilter(2));

end
