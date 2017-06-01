function ReadAnnotationFile(filePath)

if nargin < 1
    filePath = '/Users/aurina/ermineJ.data/Generic_worm_noParents2017.an.txt';
end
%-------------------------------------------------------------------------------

fprintf(1,'Reading data from %s...',filePath);
fid = fopen(filePath,'r');
% headings = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s',1,'Delimiter','\t','CollectOutput',1);
headings = textscan(fid,'%s %s %s %s %s %s',1,'Delimiter','\t','CollectOutput',1,'CommentStyle','#');
C = textscan(fid,'%s %s %s %s %s %s','Delimiter','\t','CommentStyle','#','EndOfLine','\n');
fclose(fid);

annotationTable = table();
annotationTable.ProbeName = C{1};
annotationTable.acronym = C{2};
annotationTable.name = C{3};
annotationTable.GO = C{4};
fprintf(1,' Data loaded\n');
% So now we have annotations for each gene (as rows of the annotation file)

%-------------------------------------------------------------------------------
% Get gene acronyms IDs for each GO category
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

% Now, for each GO term, list the gene entrez that are directly annotated to it:
geneEntrezAnnotations = cell(numGOCategories,1);
parfor i = 1:numGOCategories
    geneAcronymAnnotations{i} = allacronym(cellfun(@(x)ismember(allGOCategories(i),x),allGOSplitNum));
    % fprintf(1,'%u/%u\n',i,numGOCategories);
end

%-------------------------------------------------------------------------------
% Save to file:
filePath = fullfile('Data/ermineJdata/','GOAnnotation.mat');
save(filePath,'annotationTable','allGOCategories','geneAcronymAnnotations');
fprintf(1,'Saved to %s\n',filePath);

fprintf(1,'Need to run propagateHierarchy...!\n');

%-------------------------------------------------------------------------------
% function GOIDs = toNumber(GOCell)
%     GOIDs = cellfun(@(x)str2num(x(4:end)),GOCell,'UniformOutput',true);
% end

end
