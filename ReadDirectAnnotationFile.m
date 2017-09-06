function ReadDirectAnnotationFile(filePathRead)
% Can download annotation files direct from the GO website
% http://geneontology.org/page/download-annotations
% Each line in the file represents a single association between a gene product
% and a GO term with a certain evidence code and the reference to support the link.
% (http://geneontology.org/page/go-annotation-file-formats)
%
% This function processes these raw annotation files -> .mat file
%-------------------------------------------------------------------------------
% if nargin < 1
    filePathRead = 'mus_muscus_annotation.mgi';
% end
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
annotationTable.MGI_ID = C{2};
annotationTable.acronym = C{3};
annotationTable.qualifier = categorical(C{4});
annotationTable.GO = C{5};
annotationTable.evidenceCode = categorical(C{7});
annotationTable.geneName = C{10};
annotationTable.geneOrProtein = categorical(C{12});
fprintf(1,' Data loaded\n');

%-------------------------------------------------------------------------------
% Exclude NOTs:
% (could check later that they're not added to the labeled GO terms, after propagation...?)
isNOT = (annotationTable.qualifier=='NOT');
fprintf(1,'%u NOT annotations ignored\n',sum(isNOT));
annotationTable = annotationTable(~isNOT,:);

%-------------------------------------------------------------------------------
% Exclude ND evidence codes
% (cf. http://geneontology.org/page/guide-go-evidence-codes)
% "Note: The ND evidence code, unlike other evidence codes, should be considered as
% a code that indicates curation status/progress than as method used to derive an
% annotation." [http://geneontology.org/page/nd-no-biological-data-available]
isND = (annotationTable.evidenceCode=='ND');
fprintf(1,'%u ND annotations ignored\n',sum(isND));
annotationTable = annotationTable(~isND,:);

%-------------------------------------------------------------------------------
% Organize the remaining genes:
uniqueGenes = unique(annotationTable.MGI_ID);
numUniqueGenes = length(uniqueGenes);

%-------------------------------------------------------------------------------
% Map genes -> entrez IDs:
% Load mapping (generated from mousemine: cf. MGI_NCBI_downloadall.py)
fprintf(1,'Loading MGI->NCBI mapping from file...');
MGI_NCBI_Map = ImportMGI_NCBI_Map(true);
fprintf(1,' Loaded.\n');
geneEntrez = nan(numUniqueGenes,1);
for i = 1:numUniqueGenes
    isHere = strcmp(MGI_NCBI_Map.MGIID,uniqueGenes{i});
    if sum(isHere)==0
        fprintf(1,'Nothing for %s\n',uniqueGenes{i});
    elseif sum(isHere)==1
        geneEntrez(i) = MGI_NCBI_Map.NCBIGeneNumber(isHere);
    else
        error('whoa')
    end
end

%-------------------------------------------------------------------------------
% Filter out MGI annotations for genes with no Entrez ID:
MGI_hasNoEntrez = uniqueGenes(isnan(geneEntrez));
fprintf(1,'%u genes have no Entrez ID\n',sum(isnan(geneEntrez)));
annotationsNoEntrez = ismember(annotationTable.MGI_ID,MGI_hasNoEntrez);
fprintf(1,'%u annotations for genes/proteins with no NCBI (Entrez) ID are being ignored\n',...
                        sum(annotationsNoEntrez));
annotationTable = annotationTable(~annotationsNoEntrez,:);

%-------------------------------------------------------------------------------
% Filter uniqueGenes to just be those with EntrezIDs:
uniqueGenes = uniqueGenes(~isnan(geneEntrez));
geneEntrez = geneEntrez(~isnan(geneEntrez));
numUniqueGenes = length(uniqueGenes);

%-------------------------------------------------------------------------------
% Assign Entrez_ID to each annotation, and add as a column in the annotation table:
% (slow because matching on strings... :-/ Could convert to integers for faster matching)
fprintf(1,'Mapping MGI IDs to Entrez IDs across the annotation table (%u genes -> %u annotations)...',...
                        numUniqueGenes,height(annotationTable));
annotationEntrez = nan(height(annotationTable),1);
for i = 1:numUniqueGenes
    annotationsHere = strcmp(annotationTable.MGI_ID,uniqueGenes{i});
    annotationEntrez(annotationsHere) = geneEntrez(i);
end
annotationTable.EntrezID = annotationEntrez;
fprintf(1,'Done.\n');

% Hopefully this is zero:
fprintf(1,'%u = 0!!\n',sum(isnan(annotationTable.EntrezID));

%-------------------------------------------------------------------------------
% Now get list of genes annotated to each GO category:
% (slow because matching on strings; could convert GO IDs to integers for faster matching)
allGOCategories = unique(annotationTable.GO);
numGOCategories = length(allGOCategories);
fprintf(1,'%u GO terms represented... Annotating to Entrez IDs...',numGOCategories);
geneEntrezAnnotations = cell(numGOCategories,1);
parfor i = 1:numGOCategories
    geneEntrezAnnotations{i} = annotationTable.EntrezID(strcmp(annotationTable.GO,allGOCategories{i}));
end
fprintf(1,' Annotated.\n');

%-------------------------------------------------------------------------------
% Filter out categories that have no annotations
hasGOAnn = cellfun(@(x)~isempty(x),geneEntrezAnnotations);
allGOCategories = allGOCategories(hasGOAnn);
fprintf(1,'Filtered out %u GO categories with no annotations\n',sum(~hasGOAnn));

%-------------------------------------------------------------------------------
% Save to file:
fileNameSave = 'GOAnnotationDirect.mat';
save(fileNameSave,'annotationTable','allGOCategories','geneEntrezAnnotations');
fprintf(1,'Saved to %s\n',fileNameSave);

end
