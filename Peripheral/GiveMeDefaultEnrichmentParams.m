function enrichmentParams = GiveMeDefaultEnrichmentParams()

% Store everything you need in this structure:
enrichmentParams = struct();

%-------------------------------------------------------------------------------
% Generic details about the enrichment:
% (sufficient for conventional enrichment)
%-------------------------------------------------------------------------------

% What high-level GO categories to filter:
enrichmentParams.processFilter = 'biological_process';

% What range of GO category sizes to include:
enrichmentParams.sizeFilter = [10,200];

% How many samples to compute per GO category (to estimate the null distribution):
enrichmentParams.numNullSamples = 4e4;

% Display categories with corrected p-value below this threshold:
enrichmentParams.sigThresh = 0.05;

% Higher or lower scores are 'better' (compute p-value from right or left tail):
enrichmentParams.whatTail = 'right'; % higher scores (e.g., from aggregateHow) are 'better'.

%-------------------------------------------------------------------------------
% Additional details specific to the ensemble-based enrichment:
% (on the basis of spatial correlations)
%-------------------------------------------------------------------------------
% What type of correlation to use
enrichmentParams.whatCorr = 'Spearman'; % 'Pearson', 'Spearman'

% How to aggregate scores within a GO category:
enrichmentParams.aggregateHow = 'mean'; % 'mean', 'absmean', 'median', 'absmedian'

% What type of null:
enrichmentParams.whatEnsemble = 'randomMap'; % 'randomMap', 'customEnsemble'

% Specify a custom data file in the case of running 'customEnsemble' enrichment:
% (file containing the matrix of null phenotypes):
enrichmentParams.dataFileSurrogate = [];

% Map parameters on to an appropriate file name to save results to:
% (can include different mappings for different combinations of parameter variation)
enrichmentParams.fileNameOut = sprintf('PhenotypeNulls_%s_%u.mat',...
                                enrichmentParams.whatEnsemble,...
                                enrichmentParams.numNullSamples);

% Store all null computations in a consistent directory:
enrichmentParams.fileNameOut = fullfile('EnsembleEnrichment','NullEnsembles',...
                                            enrichmentParams.fileNameOut);

end
