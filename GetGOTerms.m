function GOTable = GetGOTerms(whatFilter)
% From a mySQL database, retrieves GO terms for a given process
% e.g., GOTerms = GetGoTerms('biological_process');
%-------------------------------------------------------------------------------

if nargin < 1
    whatFilter = 'biological_process';
end

dbc = SQL_opendatabase;

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

end
