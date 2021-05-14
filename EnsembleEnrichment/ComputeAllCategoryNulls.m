function GOTable = ComputeAllCategoryNulls(geneDataStruct,enrichmentParams,phenotypeVector,saveOut,beVerbose)
% ComputeAllCategoryNulls   Compute an ensemble-based null distribution for all
%                               GO categories.
%
% Can also compute for a specific phenotype by specifying a non-empty phenotypeVector
% All parameters of the caculation are set in the params structure
%   (easiest way is to set default values using GiveMeDefaultEnrichmentParams)
%
% INPUTS:
% ---geneDataStruct should be a structure containing:
%       - expressionMatrix
%       - entrezIDs (labeling genes as columns)
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs and set defaults:
%-------------------------------------------------------------------------------
if nargin < 1
    error('You must provide processed gene-expression data');
end
if nargin < 2
    enrichmentParams = GiveMeDefaultEnrichmentParams();
end
if nargin < 3
    % If a specific phenotype is specified, repeat the calculation on this single phenotype:
    phenotypeVector = [];
end
if nargin < 4
    saveOut = true;
end
if nargin < 5
    beVerbose = true;
end

%-------------------------------------------------------------------------------
% Load gene-expression data:
geneData = geneDataStruct.expressionMatrix;
entrezIDs = geneDataStruct.entrezIDs;
numAreas = size(geneData,1);

%-------------------------------------------------------------------------------
% Get a generic GO Table:
GOTable = GetFilteredGOData(enrichmentParams.dataSource,...
                                enrichmentParams.processFilter,...
                                enrichmentParams.sizeFilter,...
                                entrezIDs);
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Unpack out structure for parfor?
whatCorr = enrichmentParams.whatCorr;
numNullSamples = enrichmentParams.numNullSamples;

%-------------------------------------------------------------------------------
% Get randomized phenotypes to use as your null:
switch enrichmentParams.whatEnsemble
case 'customSpecified'
    % Specify a custom phenotype and just run the full calculation on this:
    nullMaps = phenotypeVector;
    numNullSamples = 1;
    fprintf(1,'Computing category scores for the spatial phenotype provided...\n');
case 'randomMap'
    % Generate as many random maps as null samples:
    nullMaps = rand(numAreas,numNullSamples);
    fprintf(1,'Computing category scores for %u random %u-region phenotypes\n',...
                                numNullSamples,numAreas);
case 'customEnsemble'
    % Get the pre-computed surrogate data from a comma-delimited text file:
    % nullMaps = dlmread(enrichmentParams.dataFileSurrogate,',',1,1);
    % Get the pre-computed surrogate data from the variable 'nullMaps' in a .mat file:
    load(enrichmentParams.dataFileSurrogate,'nullMaps');
    fprintf(1,'Computing category scores for %u custom-loaded %u-region phenotypes\n',...
                                    size(nullMaps,2),size(nullMaps,1));
    fprintf(1,'(Null spatial maps loaded from %s)\n',enrichmentParams.dataFileSurrogate);
otherwise
    error('Unknown null type: ''%s''',whatEnsemble);
end

%-------------------------------------------------------------------------------
% Prepare for category-wise agglomeration by first computing the correlation of
% each gene with a given spatial map (or a null ensemble of spatial maps)
allAnnotatedGenesEntrez = unique(vertcat(GOTable.annotations{:}));
numAnnotatedTotal = length(allAnnotatedGenesEntrez);
allAnnotatedGenesEntrez = intersect(allAnnotatedGenesEntrez,entrezIDs);
numAnnotatedGenes = length(allAnnotatedGenesEntrez);
fprintf(1,'Of %u annotated genes, %u match genes we have expression data for\n',...
                    numAnnotatedTotal,numAnnotatedGenes);

geneScores = nan(numAnnotatedGenes,numNullSamples);
fprintf(1,'Computing null distributions to %u null phenotypes for all %u genes annotated to GO categories\n',...
                            numNullSamples,numAnnotatedGenes);
for g = 1:numAnnotatedGenes
    % Get this gene's expression vector:
    matchMe = (entrezIDs==allAnnotatedGenesEntrez(g));
    if sum(matchMe)~=1
        fprintf(1,'How the heck did I get %u matches for gene entrez %u?!\n',...
                            sum(matchMe),allAnnotatedGenesEntrez(g));
    end
    expressionVector = geneData(:,matchMe);

    % The correlation to compute for this gene:
    theCorr_fun = @(x) corr(x,expressionVector,'type',whatCorr,'rows','pairwise');
    if numNullSamples==1
        geneScores(g) = theCorr_fun(nullMaps(:,1));
    else
        parfor n = 1:numNullSamples
            geneScores(g,n) = theCorr_fun(nullMaps(:,n));
        end
    end
end

%-------------------------------------------------------------------------------
% Agglomerate gene scores by GO category
%-------------------------------------------------------------------------------
categoryScores = cell(numGOCategories,1);
for i = 1:numGOCategories
    if beVerbose
        fprintf(1,'\n\n\nLooking in at Category %u/%u. %s:%s (%u)\n',...
            i,numGOCategories,GOTable.GOIDlabel{i},GOTable.GOName{i},GOTable.size(i));
    end

    % Match genes for this category:
    categoryEntrez = GOTable.annotations{i};
    matchMe = find(ismember(allAnnotatedGenesEntrez,categoryEntrez));
    numGenesCategory = length(matchMe);

    if beVerbose
        fprintf(1,'%u/%u genes from this GO category have matching records in the expression data\n',...
                    numGenesCategory,length(theGenesEntrez));
    end

    % Retrieve the distribution of gene scores across phenotypes:
    scoresHere = geneScores(matchMe,:);

    %---------------------------------------------------------------------------
    % Aggregate gene-wise scores into an overall GO category score (for each phenotype)
    categoryScores{i} = AggregateScores(scoresHere,enrichmentParams.aggregateHow);
end

%-------------------------------------------------------------------------------
% Assign to the table
GOTable.categoryScores = categoryScores;

%-------------------------------------------------------------------------------
% Save results to .mat file
if saveOut
    fprintf(1,'Saving %s nulls from %u iterations to ''%s''\n',...
                    enrichmentParams.whatEnsemble,...
                    numNullSamples,...
                    enrichmentParams.fileNameOut);
    save(enrichmentParams.fileNameOut,'GOTable','enrichmentParams','-v7.3');
end

end
