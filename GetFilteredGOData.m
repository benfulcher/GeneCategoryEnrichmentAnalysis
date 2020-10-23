function GOTable = GetFilteredGOData(whatSource,whatFilter,sizeFilter,ourEntrez)
% GetFilteredGOData       Returns a GO table of gene annotations
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs:
if nargin < 1
    whatSource = 'mouse-direct'; % Direct annotations from GO
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
% cf. propagateHierarchy to map files generated from ReadDirectAnnotationFile (or ReadGEMMAAnnotationFile)
switch whatSource
case 'mouse-direct'
    fileNameLoad = sprintf('GOAnnotationDirect-mouse-%s-Prop.mat',whatFilter);
case 'human-direct'
    fileNameLoad = sprintf('GOAnnotationDirect-human-%s-Prop.mat',whatFilter);
case 'mouse-GEMMA'
    fileNameLoad = 'GOAnnotationGEMMAProp.mat';
otherwise
    error('Unknown annotation source: ''%s''',whatSource);
end

load(fileNameLoad,'GOTerms');
fprintf(1,'Loaded annotated GO Terms from %s\n',fileNameLoad);

%-------------------------------------------------------------------------------
% Get full GO ontology details:
GOTermsFile = 'GOTerms_BP.mat';
if strcmp(whatFilter,'biological_process') && exist(GOTermsFile,'file')
    load(GOTermsFile,'GOTable');
    fprintf(1,'Loaded biological process GO terms from file: %s\n',GOTermsFile);
else
    % Retrieve from mySQL database:
    GOTable = GetGOTerms(whatFilter);
    fprintf(1,'Loaded biological process GO terms from SQL database\n');
end

%-------------------------------------------------------------------------------
% Filter categories
%-------------------------------------------------------------------------------
% Filter by ontology details:
[~,ia,ib] = intersect(GOTable.GOID,GOTerms.GOID);
fprintf(1,'Filtering to %u annotated GO categories related to %s\n',length(ia),whatFilter);
GOTable = GOTable(ia,:);
GOTerms = GOTerms(ib,:);

% Filter by category size:
numGOCategories = height(GOTerms);
if isempty(ourEntrez)
    fprintf(1,'Filtering GO categories on their number of gene annotations\n');
    sizeGOCategories = GOTerms.size;
    fprintf(1,'%u GO categories have no annotations :-/\n',...
                    sum(sizeGOCategories==0));
else
    fprintf(1,'Filtering on size of GO category after matching to our %u genes\n',length(ourEntrez));
    % Make sure all annotations match those in our set of entrez_ids
    GOTerms.annotations = cellfun(@(x)x(ismember(x,ourEntrez)),GOTerms.annotations,...
                                        'UniformOutput',false);
    sizeGOCategories = cellfun(@length,GOTerms.annotations);
    GOTerms.size = sizeGOCategories;
    fprintf(1,'%u GO categories have no annotations matching our %u genes\n',...
                    sum(sizeGOCategories==0),length(ourEntrez));
end

GOTable.size = GOTerms.size;
GOTable.annotations = GOTerms.annotations;
fprintf(1,'GO categories have between %u and %u annotations\n',...
            min(sizeGOCategories),max(sizeGOCategories));
isGoodSize = (sizeGOCategories >= sizeFilter(1)) & (sizeGOCategories <= sizeFilter(2));
GOTable = GOTable(isGoodSize,:);
fprintf(1,'Filtered to %u GO categories with between %u and %u annotations\n',...
                height(GOTable),sizeFilter(1),sizeFilter(2));


%-------------------------------------------------------------------------------
% Check that all annotations in GOTable are column vectors:
isRowVector = find(cellfun(@(x)size(x,2),GOTable.annotations)~=1);
numRowVectorAnnotation = length(isRowVector);
for i = 1:numRowVectorAnnotation
    GOTable.annotations{isRowVector(i)} = GOTable.annotations{isRowVector(i)}';
end
fprintf(1,'Transposed %u category annotation vectors from row-vector -> column vector\n',numRowVectorAnnotation);

end
