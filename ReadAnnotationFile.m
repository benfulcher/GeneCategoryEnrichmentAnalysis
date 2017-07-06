function annotationTable = ReadAnnotationFile(filePath)

if nargin < 1
    filePath = '/Users/aurina/ermineJ.data/Generic_worm_2017.an.txt';
end
%-------------------------------------------------------------------------------

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

%-------------------------------------------------------------------------------
% Get gene acronyms IDs for each GO category
allGO = annotationTable.GO;
hasGOAnn = cellfun(@(x)~isempty(x),allGO);
allGO = allGO(hasGOAnn);
allacronym = annotationTable.acronym(hasGOAnn);
% Split on pipe:
allGOSplit = regexp(allGO,'\|','split');
toNumber = @(GOCell) cellfun(@(x)str2num(x(4:end)),GOCell,'UniformOutput',true);
allGOSplitNum = cellfun(@(x)toNumber(x),allGOSplit,'UniformOutput',false);
allGOTogether = horzcat(allGOSplitNum{:});
allGOCategories = unique(allGOTogether);
numGOCategories = length(allGOCategories);

% Now, for each GO term, list the gene acronym that are annotated to it:
geneAcronymAnnotations = cell(numGOCategories,1);
parfor i = 1:numGOCategories
    geneAcronymAnnotations{i} = allacronym(cellfun(@(x)ismember(allGOCategories(i),x),allGOSplitNum));
    % fprintf(1,'%u/%u\n',i,numGOCategories);
end

%-------------------------------------------------------------------------------
% Save to file:
filePath = fullfile('Data/ermineJdata/','GOAnnotationWithParents.mat');
save(filePath,'annotationTable','allGOCategories','geneAcronymAnnotations');
fprintf(1,'Saved to %s\n',filePath);


%-------------------------------------------------------------------------------
% function GOIDs = toNumber(GOCell)
%     GOIDs = cellfun(@(x)str2num(x(4:end)),GOCell,'UniformOutput',true);
% end

end
