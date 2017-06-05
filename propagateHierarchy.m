function propagateHierarchy(whatFilter)

if nargin < 1
    whatFilter = 'biological_process';
end

load('GOAnnotation.mat','annotationTable','allGOCategories','geneEntrezAnnotations');

GOTable = GetGOTerms(whatFilter);

%-------------------------------------------------------------------------------
% mySQL query for hierarchy:
%-------------------------------------------------------------------------------
connSettings.hostname = 'localhost:1234';
connSettings.dbname = 'GODaily';
connSettings.username = 'benfulcher';
connSettings.password = 'ben';
dbc = SQL_opendatabase(connSettings);

% e.g., Get all biological_process-tagged GO categories
% selectText = sprintf('SELECT term1_id,term2_id FROM term2term WHERE relationship_type_id = 1');
% selectText1 = sprintf(['SELECT term2term.term1_id,term2term.term2_id FROM term2term INNER JOIN',...
%             ' term ON term2term.term1_id=term.id WHERE term2term.relationship_type_id=1',...
%             ' AND term.term_type LIKE "%s"'],whatFilter);
% selectText2 = sprintf(['SELECT term2term.term1_id,term2term.term2_id FROM term2term INNER JOIN',...
%             ' term ON term2term.term2_id=term.id WHERE term2term.relationship_type_id=1',...
%             ' AND term.term_type LIKE "%s"'],whatFilter);
% tableResults1 = mysql_dbquery(dbc,selectText1);
% tableResults2 = mysql_dbquery(dbc,selectText2);
selectText1 = ['SELECT term2term.id,term.acc FROM term INNER JOIN term2term ON ',...
                'term.id=term2term.term1_id WHERE term2term.relationship_type_id IN (1,25,27) ',...
                'AND term.acc LIKE "GO:%"'];
hierarchyRel1 = mysql_dbquery(dbc,selectText1);
selectText2 = ['SELECT term2term.id,term.acc FROM term INNER JOIN term2term ON ',...
                'term.id=term2term.term2_id WHERE term2term.relationship_type_id IN (1,25,27) ',...
                'AND term.acc LIKE "GO:%"'];
hierarchyRel2 = mysql_dbquery(dbc,selectText2);
%1==is_a, 25==part_of, 27==regulates
SQL_closedatabase(dbc);

id1 = vertcat(hierarchyRel1{:,1});
id2 = vertcat(hierarchyRel2{:,1});
[idBoth,ix1,ix2] = intersect(id1,id2);
hierarchyRel = [vertcat(cellfun(@(x)str2num(x(4:end)),hierarchyRel1(ix1,2))),...
                vertcat(cellfun(@(x)str2num(x(4:end)),hierarchyRel2(ix2,2)))];
% So we have pairwise *is_a* interactions in hierarchyRel

%===============================================================================
% Now filter on terms that exist in our set of BP GO terms
yeahBP = ismember(hierarchyRel,GOTable.GOID);
isBP = all(yeahBP,2);
fprintf(1,'%.2f%% of hierarchical relationships are related to our categories\n',mean(isBP)*100);
hierarchyRel = hierarchyRel(isBP,:);
numGOTerms = height(GOTable);
fprintf(1,'So we have %u hierachical relationships between %u GO terms\n',length(hierarchyRel),numGOTerms);

% Convert to pairwise matrix
hierarchyMatrix = zeros(numGOTerms,numGOTerms);
for i = 1:size(hierarchyRel,1)
    hierarchyMatrix(GOTable.GOID==hierarchyRel(i,1),GOTable.GOID==hierarchyRel(i,2)) = 1;
end

% Propagate:
numGOCategories = length(allGOCategories);
for j = 1:numGOCategories
    % Get parents of category using the full GO hierarchy
    parentIDs = PropagateUp(find(GOTable.GOID==allGOCategories(j)),hierarchyMatrix);
    % Filter to those with annotated terms:
    idx = find(ismember(allGOCategories,GOTable.GOID(parentIDs)));
    % Add terms to parent
    for k = 1:length(idx) % loop over parents
        geneEntrezAnnotations{idx(k)} = union(geneEntrezAnnotations{idx(k)},geneEntrezAnnotations{j});
        % Add annotations of child to all hierarchical parents
    end
    fprintf(1,'%u/%u\n',j,numGOCategories);
end

%-------------------------------------------------------------------------------
% Save
hierarchyMatrix = sparse(hierarchyMatrix);
fileNameSave = 'GOAnnotationProp.mat';
save(fileNameSave,'annotationTable','allGOCategories','geneEntrezAnnotations','hierarchyMatrix');
fprintf(1,'Saved to %s\n',fileNameSave);

end
