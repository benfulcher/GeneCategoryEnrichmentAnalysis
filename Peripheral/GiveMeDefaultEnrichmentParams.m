function params = GiveMeDefaultEnrichmentParams()
% GiveMeDefaultEnrichmentParams   Returns default GSR parameters
%-------------------------------------------------------------------------------

params = struct();

% Where to get annotation information from
params.dataSource = 'mouse-direct';

% What high-level GO categories to filter:
params.processFilter = 'biological_process';

% What range of GO category sizes to include:
params.sizeFilter = [5,200];

% How many samples to generate to estimate the null distribution:
params.numSamples = 1e4;

end
