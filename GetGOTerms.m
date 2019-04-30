function GOTable = GetGOTerms(whatFilter,doSave)
% From a mySQL database, retrieves GO terms for a given process
% e.g., GOTerms = GetGOTerms('biological_process');
%-------------------------------------------------------------------------------

if nargin < 1
    whatFilter = 'biological_process';
end
if nargin < 2
    doSave = false;
end
%-------------------------------------------------------------------------------

dbc = ConnectMeDatabase();

% e.g., Get all biological_process-tagged GO categories
selectText = sprintf('SELECT acc,name FROM term WHERE term_type LIKE ''%s''',whatFilter);
tableResults = mysql_dbquery(dbc,selectText);
SQL_closedatabase(dbc);

%-------------------------------------------------------------------------------

GOTable = table();
GOTable.GOIDlabel = tableResults(:,1);
GOTable.GOName = tableResults(:,2);

% Filter by GO categories with valid labels:
isGOCat = cellfun(@(x)strcmp(x(1:3),'GO:'),GOTable.GOIDlabel);
GOTable(~isGOCat,:) = [];

GOTable.GOID = cellfun(@(x)str2num(x(4:end)),GOTable.GOIDlabel);

% Sort:
GOTable = sortrows(GOTable,'GOID');

%===============================================================================
% Save to file
if doSave
    if strcmp(whatFilter,'biological_process')
        GOTermFile = fullfile('ProcessedData','GOTerms_BP.mat');
        save(GOTermFile,'GOTable');
        fprintf(1,'Saved filtered GO table to %s\n',GOTermFile);
    else
        fprintf(1,'Don''t know how to save data that isn''t a simple BP filter...\n');
    end
end


end
