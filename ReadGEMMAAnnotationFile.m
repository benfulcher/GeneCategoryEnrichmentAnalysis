function ReadGEMMAAnnotationFile(filePathRead)
% Reads in a GO annotation file (noParents from GEMMA) and outputs the result
% to a Matlab file: GOAnnotation.mat
% Details of the annotation file format are here:
% http://erminej.chibi.ubc.ca/help/input-files/gene-annotations/
% Description of the format

% The file is tab-delimited text. Comma-delimited files or Excel spreadsheets (for example) are not supported.
% There is a one-line header included in the file for readability.
% ---The first column contains the probe identifier. The probe IDs must exactly
% match the ones you provide in your Gene score file.
% Any probes not having an entry will be ignored.
% If you are not using probes, this will probably contain gene symbols.
% The main requirement here is that it matches the identifiers you provide in your input data files.
% ---The second column usually contains a gene symbol. This should not be blank.
% If the gene name is not known, a sequence identifier or arbitrary code can be
% used instead. This is used to determine whether a gene has more than one probe,
% as well as providing information for display purposes.
% ---The third column contains the gene name (or description). This can be blank.
% It is only used for display purposes.
% ---The fourth column contains a delimited list of GO identifiers.
% These include the “GO:” prefix. Thus they read “GO:00494494” and not “494494”.
% The ids within this field can be delimited by spaces, commas, or pipe (‘|’) symbols.
% This field can be blank if there are no GO annotations (or if you aren’t using GO).
%-------------------------------------------------------------------------------
if nargin < 1
    filePathRead = '/Users/benfulcher/ermineJ.data/Generic_mouse_ncbiIds_noParents_Nov2016.an.txt';
end
%-------------------------------------------------------------------------------

fprintf(1,'Reading data from %s...',filePathRead);
fid = fopen(filePathRead,'r');
% headings = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s',1,'Delimiter','\t','CollectOutput',1);
headings = textscan(fid,'%s %s %s %s %s %s',1,'Delimiter','\t','CollectOutput',1,'CommentStyle','#');
C = textscan(fid,'%u %s %s %s %s %s','Delimiter','\t','CommentStyle','#','EndOfLine','\n');
fclose(fid);

% Represent as a Matlab table object (of form entrez,acronym,name,GO-annotations)
annotationTable = table();
annotationTable.entrez_id = C{1};
annotationTable.acronym = C{2};
annotationTable.name = C{3};
annotationTable.GOID = C{4};
fprintf(1,' Data loaded\n');

%-------------------------------------------------------------------------------
% Get entrez IDs for each GO category
allGO = annotationTable.GOID;
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
