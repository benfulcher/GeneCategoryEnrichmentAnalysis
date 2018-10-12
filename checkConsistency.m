% checkConsistency: check how consistent results are from our random-gene null
% enrichment method, and ermineJ's enrichment results
%-------------------------------------------------------------------------------


[GOTable,geneEntrezAnnotations] = GetFilteredGOData('biological_process',[0,100],geneEntrezIDs);

i = 3;
match = strcmp(ermineJResults.GOID,GOTable.GOIDlabel{i});
fprintf(1,'%s:%s -- %u genes from ermineJ, %u from matlab-based analysis\n',...
                        GOTable.GOIDlabel{i},...
                        GOTable.GOName{i},...
                        ermineJResults.numGenes(match),...
                        GOTable.size(i));
ermineJResults.geneMembers{match}
matlabGenes = cell(GOTable.size(i),1);
for j = 1:GOTable.size(i)
    ij = geneInfo.entrez_id==geneEntrezAnnotations{i}(j);
    matlabGenes{j} = geneInfo.acronym{ij};
end
matlabGenes = sort(matlabGenes);
for j = 1:GOTable.size(i)
    fprintf(1,'%s,',matlabGenes{j});
end
fprintf(1,'\n');
