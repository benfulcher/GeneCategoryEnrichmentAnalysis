function GOTable = ComputeAllCategoryNulls(geneDataStruct,enrichmentParams,phenotypeVector,saveOut,beVerbose)
% ComputeAllCategoryNulls   Compute an ensemble-based null distribution for all
%                               GO categories.
%
% Can also compute for a specific phenotype by specifying a non-empty phenotypeVector
% All parameters of the caculation are set in the params structure
%   (easiest way is to set default values using GiveMeDefaultEnsembleParams)
%
% geneDataStruct should be a structure containing:
%   - expressionMatrix
%   - entrezIDs (labeling genes as columns)
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
numAreas = size(geneData,1);
entrezIDs = geneDataStruct.entrezIDs;

%-------------------------------------------------------------------------------
% Get a generic GO Table:
GOTable = GetFilteredGOData(enrichmentParams.dataSource,...
                                enrichmentParams.processFilter,...
                                enrichmentParams.sizeFilter,...
                                entrezIDs);
numGOCategories = height(GOTable);

%-------------------------------------------------------------------------------
% Get random vectors from real genes to use as null spatial maps:
switch enrichmentParams.whatEnsemble
case 'customSpecified'
    % Specify a custom phenotype and just run the full calculation on this:
    nullMaps = phenotypeVector;
    enrichmentParams.numNullSamples = 1;
    fprintf(1,'Computing category scores for the spatial phenotype provided\n');
case 'randomMap'
    % Generate as many random maps as null samples:
    nullMaps = rand(numAreas,enrichmentParams.numNullSamples);
    fprintf(1,'Computing category scores for %u random %u-region phenotypes\n',...
                                enrichmentParams.numNullSamples,numAreas);
case 'customEnsemble'
    % Get the pre-computed surrogate data from a comma-delimited text file:
    nullMaps = dlmread(enrichmentParams.dataFileSurrogate,',',1,1);
    fprintf(1,'Computing category scores for %u custom-loaded %u-region phenotypes\n',...
                                    size(nullMaps,2),size(nullMaps,1));
    fprintf(1,'(Null spatial maps loaded from %s)\n',enrichmentParams.dataFileSurrogate);
otherwise
    error('Unknown null type: ''%s''',whatEnsemble);
end

%-------------------------------------------------------------------------------
% Correlation of genes with a given spatial map (or null ensemble of spatial maps):
categoryScores = cell(numGOCategories,1);
for i = 1:numGOCategories
    if beVerbose
        fprintf(1,'\n\n\nCategory %u/%u\n',i,numGOCategories);

        fprintf(1,'Looking in at %s:%s (%u)\n',GOTable.GOIDlabel{i},...
                            GOTable.GOName{i},GOTable.size(i));
    end

    % Match genes for this category:
    theGenesEntrez = GOTable.annotations{i};
    matchMe = find(ismember(entrezIDs,theGenesEntrez));
    geneDataCategory = geneData(:,matchMe);
    numGenesCategory = length(matchMe);

    if beVerbose
        fprintf(1,'%u/%u genes from this GO category have matching records in the expression data\n',...
                            length(matchMe),length(theGenesEntrez));
    end

    % Compute the distribution of gene-category scores for correlation with the null maps:
    scoresHere = nan(numGenesCategory,enrichmentParams.numNullSamples);
    for k = 1:numGenesCategory
        expressionVector = geneDataCategory(:,k);
        if enrichmentParams.numNullSamples==1
            scoresHere(k) = corr(nullMaps(:,1),expressionVector,'type',enrichmentParams.whatCorr,'rows','pairwise');
        else
            parfor j = 1:enrichmentParams.numNullSamples
                scoresHere(k,j) = corr(nullMaps(:,j),expressionVector,'type',enrichmentParams.whatCorr,'rows','pairwise');
            end
        end
    end

    % Aggregate gene-wise scores into an overall GO category score
    switch enrichmentParams.aggregateHow
    case 'mean'
        categoryScores{i} = nanmean(scoresHere,1);
    case 'median'
        categoryScores{i} = nanmedian(scoresHere,1);
    otherwise
        error('Unknown aggregation option: ''%s''',enrichmentParams.aggregateHow);
    end
end

%-------------------------------------------------------------------------------
% Assign to the table:
GOTable.categoryScores = categoryScores;

%-------------------------------------------------------------------------------
% Save results to .mat file
if saveOut
    fprintf(1,'Saving %s nulls from %u iterations to ''%s''\n',...
                    enrichmentParams.whatEnsemble,...
                    enrichmentParams.numNullSamples,...
                    enrichmentParams.fileNameOut);
    save(enrichmentParams.fileNameOut,'GOTable','enrichmentParams','-v7.3');
end

end
