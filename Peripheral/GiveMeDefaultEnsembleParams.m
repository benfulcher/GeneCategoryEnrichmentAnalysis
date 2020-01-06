function params = GiveMeDefaultEnsembleParams()

params = struct();

% What type of correlation to use
params.whatCorr = 'Spearman'; % 'Pearson', 'Spearman'

% How to agglomerate scores within a GO category:
params.aggregateHow = 'mean'; % 'mean'

% How many null samples to compute per GO category:
params.numNullSamples = 40000;

% What type of null:
params.whatEnsemble = 'randomMap';

% Specify a custom data file in the case of running 'spatialLag' enrichment:
params.dataFileSurrogate = [];

end
