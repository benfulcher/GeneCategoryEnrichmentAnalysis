function ReadDirectAnnotationFile(filePathRead)
% Can download annotation files direct from the GO website
% http://geneontology.org/page/download-annotations
% Each line in the file represents a single association between a gene product
% and a GO term with a certain evidence code and the reference to support the link.
% (http://geneontology.org/page/go-annotation-file-formats)
%
% This function processes these raw annotation files -> .mat file
%-------------------------------------------------------------------------------
if nargin < 1
    filePathRead = 'mus_muscus_annotation.mgi';
end
%-------------------------------------------------------------------------------

fprintf(1,'Reading data from %s...',filePathRead);
fid = fopen(filePathRead,'r');
% 1.DB, 2.DB Object ID, 3.DB ObjectSymbol, 4.Qualifier, 5.GOID, 6.DBReference(s),
% 7.Evidence Code, 8.With/from, 9.Aspect (GO DAG: F,P,C), 10.DBObjectName,
% 11.DBObject Synonym(s), 12.DBObject Type, 13.Taxon, 14.Date, 15.Assigned By
% 16.Annotation Extension, 17.Gene Product Form ID
C = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','CommentStyle','!','EndOfLine','\n');
fclose(fid);

% Represent as a Matlab table object (of form entrez,acronym,name,GO-annotations)
annotationTable = table();
annotationTable.MGI_ID = C{2};
annotationTable.acronym = C{3};
annotationTable.qualifier = categorical(C{4});
annotationTable.GO = C{5};
annotationTable.evidenceCode = categorical(C{7});
annotationTable.geneName = C{10};
annotationTable.geneOrProtein = categorical(C{12});
fprintf(1,' Data loaded\n');

%-------------------------------------------------------------------------------
% Exclude NOTs:
% (could check later that they're not added to the labeled GO terms, after propagation...?)
isNOT = (annotationTable.qualifier=='NOT');
fprintf(1,'%u NOT annotations ignored\n',sum(isNOT));
annotationTable = annotationTable(~isNOT,:);

%-------------------------------------------------------------------------------
% Exclude ND evidence codes
% (cf. http://geneontology.org/page/guide-go-evidence-codes)
% "Note: The ND evidence code, unlike other evidence codes, should be considered as
% a code that indicates curation status/progress than as method used to derive an
% annotation." [http://geneontology.org/page/nd-no-biological-data-available]
isND = (annotationTable.evidenceCode=='ND');
fprintf(1,'%u ND annotations ignored\n',sum(isND));
annotationTable = annotationTable(~isND,:);

%-------------------------------------------------------------------------------
% Map genes -> entrez IDs:
uniqueGenes = unique(annotationTable.acronym);
numUniqueGenes = length(uniqueGenes);
geneEntrez = zeros(numUniqueGenes,1);
for i = 1:numUniqueGenes
    geneEntrez(i) = GiveMeEntrezID(uniqueGenes{i},'mouse');
end

%-------------------------------------------------------------------------------
% Get list of genes annotated to each GO category:
allGOCategories = unique(annotationTable.GO);
numGOCategories = length(allGO);
fprintf(1,'%u GO terms represented\n',numGOCategories);
geneAnnotations = cell(numGOCategories,1);
for i = 1:numGOCategories
    geneAnnotations{i} = annotationTable.acronym(strcmp(annotationTable.GO,allGOCategories{i}));
end

hasGOAnn = cellfun(@(x)~isempty(x),allGO);
allGO = allGO(hasGOAnn);
allEntrez = annotationTable.entrez_id(hasGOAnn);

% Split on pipe:
allGOSplit = regexp(allGO,'\|','split');
toNumber = @(GOCell) cellfun(@(x)str2num(x(4:end)),GOCell,'UniformOutput',true);
allGOSplitNum = cellfun(@(x)toNumber(x),allGOSplit,'UniformOutput',false);
allGOTogether = horzcat(allGOSplitNum{:});
allGOCategories = unique(allGOTogether);
numGOCategories = length(allGOCategories);

% Now, for each GO term, list the gene entrez that are annotated to it:
geneEntrezAnnotations = cell(numGOCategories,1);
parfor i = 1:numGOCategories
    geneEntrezAnnotations{i} = allEntrez(cellfun(@(x)ismember(allGOCategories(i),x),allGOSplitNum));
    % fprintf(1,'%u/%u\n',i,numGOCategories);
end

%-------------------------------------------------------------------------------
% Save to file:
fileNameSave = 'GOAnnotationGEMMA.mat';
save(fileNameSave,'annotationTable','allGOCategories','geneEntrezAnnotations');
fprintf(1,'Saved to %s\n',fileNameSave);


end
