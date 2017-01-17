function filePath = writeErmineJFile(whatData,geneMeasures,theGeneEntrez,columnName)
% Writes a tab-delimited input file for ermineJ
% Ben Fulcher, 2014-10-14
% ------------------------------------------------------------------------------

if nargin < 4
    columnName = 'geneMeanP';
end

%-------------------------------------------------------------------------------
fileName = sprintf('ermineJ_%s.txt',whatData);
filePath = fullfile('DataOutputs','ermineJ',fileName);

% Delete the enrichment file if it already exists:
if exist(filePath,'file')==2
    delete(filePath);
    fprintf(1,'Deleted %s\n',filePath);
end

%-------------------------------------------------------------------------------
fid = fopen(filePath,'w');
numGenes = length(theGeneEntrez);

% 1. The header
fprintf(fid,'%s\t%s\n','Gene',columnName);

% 2. The gene list with p-values.
for i = 1:numGenes
    fprintf(fid,'%u\t%.8g\n',theGeneEntrez(i),geneMeasures(i));
    % fprintf(fid,'%s\t%f\n',theGeneStruct(i).gene_acronym,meanP(i));
end

fclose(fid);

% ------------------------------------------------------------------------------
% Display result:
fprintf(1,'\nWrote %s to file for ermineJ.\n\n',filePath);

end
