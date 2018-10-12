function nullDistribution = PermuteForNullDistributions(geneScores,uniqueSizes,numSamples)
% Generate a null distribution for each size of GO category
%-------------------------------------------------------------------------------

numSizes = length(uniqueSizes);
numGenes = length(geneScores);

nullDistribution = zeros(numSizes,numSamples);
for j = 1:numSamples
    rp = randperm(numGenes);
    for i = 1:numSizes
        nullDistribution(i,j) = nanmean(geneScores(rp(1:uniqueSizes(i))));
    end
end

%-------------------------------------------------------------------------------
% Inspect to check:
% Nulls should be wider for smaller categories (greater sampling variance)
% (hSmall should have greater variance than hBig):
% f = figure('color','w'); hold on
% hSmall = histogram(nullDistribution(1,:));
% hBig = histogram(nullDistribution(end,:));
% legend([hSmall,hBig],{num2str(uniqueSizes(1)),num2str(uniqueSizes(end))});

end
