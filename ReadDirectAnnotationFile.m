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
annotationTable.acronym = C{3};
annotationTable.qualifier = C{4};
annotationTable.GO = C{4};
fprintf(1,' Data loaded\n');

%-------------------------------------------------------------------------------
% Get entrez IDs for each GO category
allGO = annotationTable.GO;
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


%-------------------------------------------------------------------------------
% function GOIDs = toNumber(GOCell)
%     GOIDs = cellfun(@(x)str2num(x(4:end)),GOCell,'UniformOutput',true);
% end

end
