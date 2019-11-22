function params = GiveMeDefaultEnsembleParams()

params = struct();

% What type of correlation to use
params.whatCorr = 'Spearman'; % 'Pearson', 'Spearman'

% How to agglomerate scores within a GO category:
params.aggregateHow = 'mean'; % 'mean'

% How many null samples to compute per GO category:
params.numNullSamples = 40000;

% What type of null:
params.whatNullType = 'randomMap';

% What data to use:
params.whatSpecies = 'mouse';
params.structFilter = 'all';

end
