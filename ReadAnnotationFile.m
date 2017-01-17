function ReadAnnotationFile(filePath)

if nargin < 1
    filePath = '/Users/benfulcher/ermineJ.data/Generic_mouse_ncbiIds_noParents_Nov2016.an.txt';
end
%-------------------------------------------------------------------------------

fid = fopen(filePath,'r');
% headings = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s',1,'Delimiter','\t','CollectOutput',1);
headings = textscan(fid,'%s %s %s %s %s %s',1,'Delimiter','\t','CollectOutput',1,'CommentStyle','#');
C = textscan(fid,'%u %s %s %s %s %s','Delimiter','\t','CommentStyle','#','EndOfLine','\n');
fclose(fid);

annotationTable = table();
annotationTable.entrez_id = C{1};
annotationTable.acronym = C{2};
annotationTable.name = C{3};
annotationTable.GO = C{4};

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
filePath = fullfile('DataOutputs','GOAnnotation.mat');
save(filePath,'annotationTable','allGOCategories','geneEntrezAnnotations');
fprintf(1,'Saved to %s\n',filePath);


%-------------------------------------------------------------------------------
% function GOIDs = toNumber(GOCell)
%     GOIDs = cellfun(@(x)str2num(x(4:end)),GOCell,'UniformOutput',true);
% end

end
