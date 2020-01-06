function GOTable = ComputeAllCategoryNulls(params,phenotypeVector,saveOut,beVerbose)
% ComputeAllCategoryNulls   Compute an ensemble-based null distribution for all
%                               GO categories.
%
% Can also compute for a specific phenotype by specifying a non-empty phenotypeVector
% All parameters of the caculation are set in the params structure
%   (easiest way is to set default values using GiveMeDefaultEnsembleParams)
%
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs and Set Defaults:
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultEnsembleParams();
end
if nargin < 2
    % If a specific phenotype is specified, repeat the calculation on this single phenotype:
    phenotypeVector = [];
end
if nargin < 3
    saveOut = true;
end
if nargin < 4
    beVerbose = true;
end

%-------------------------------------------------------------------------------
% Get real data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);
numAreas = height(structInfo);

%-------------------------------------------------------------------------------
% Get a generic GO Table:
GOTable = GiveMeGOData(params,geneInfo.entrez_id);
numAreas = height(structInfo);
numGOCategories = height(GOTable);
numGenesReal = height(geneInfo);

%-------------------------------------------------------------------------------
% Get random vectors from real genes to use as null spatial maps:
switch params.whatEnsemble
case 'customSpecified'
    % Specify a custom phenotype and just run the full calculation on this:
    nullMaps = phenotypeVector;
    params.numNullSamples = 1;
case 'randomMap'
    % Generate as many random maps as null samples:
    nullMaps = rand(numAreas,numNullSamples);
case 'customEnsemble'
    % Get the pre-computed surrogate data:
    nullMaps = dlmread(dataFileSurrogate,',',1,1);
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
    matchMe = find(ismember(geneInfo.entrez_id,theGenesEntrez));
    geneDataCategory = geneData(:,matchMe);
    numGenesCategory = length(matchMe);

    if beVerbose
        fprintf(1,'%u/%u genes from this GO category have matching records in the expression data\n',...
                            length(matchMe),length(theGenesEntrez));
    end

    % Compute the distribution of gene-category scores for correlation with the null maps:
    scoresHere = nan(numGenesCategory,params.numNullSamples);
    for k = 1:numGenesCategory
        expressionVector = geneDataCategory(:,k);
        if params.numNullSamples==1
            scoresHere(k) = corr(nullMaps(:,1),expressionVector,'type',params.whatCorr,'rows','pairwise');
        else
            parfor j = 1:params.numNullSamples
                scoresHere(k,j) = corr(nullMaps(:,j),expressionVector,'type',params.whatCorr,'rows','pairwise');
            end
        end
    end

    % Aggregate gene-wise scores into an overall GO category score
    switch params.aggregateHow
    case 'mean'
        categoryScores{i} = nanmean(scoresHere,1);
    case 'median'
        categoryScores{i} = nanmedian(scoresHere,1);
    otherwise
        error('Unknown aggregation option: ''%s''',params.aggregateHow);
    end
end

%-------------------------------------------------------------------------------
GOTable.categoryScores = categoryScores;

%-------------------------------------------------------------------------------
% Save results to .mat file
if saveOut
    fileNameOut = sprintf('RandomNull_%u_%s-%s_%s_%s_%s.mat',params.numNullSamples,whatSpecies,...
                                    params.g.structFilter,params.whatEnsemble,whatCorr,aggregateHow);
    fileNameOut = fullfile('DataOutputs',fileNameOut);
    save(fileNameOut,'GOTable','params','-v7.3');
    fprintf(1,'Results of %u iterations saved to %s\n',params.numNullSamples,fileNameOut);
end

end
