function categoryScore = AggregateScores(geneScores,aggregateHow)
% Aggregates a set of scores of genes in a category in a category-wide score
%
% geneScores are aggregated down columns (assuming each phenotype is a column)
%-------------------------------------------------------------------------------

if nargin < 2
    aggregateHow = 'mean';
end

%---------------------------------------------------------------------------
% Aggregate gene-wise scores into an overall GO category score (for each phenotype)
switch aggregateHow
case 'mean'
    categoryScore = nanmean(geneScores,1);
case 'absmean'
    categoryScore = nanmean(abs(geneScores),1);
case 'median'
    categoryScore = nanmedian(geneScores,1);
case 'absmedian'
    categoryScore = nanmedian(abs(geneScores),1);
otherwise
    error('Unknown aggregation option: ''%s''',aggregateHow);
end

end
