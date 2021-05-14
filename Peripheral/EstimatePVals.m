function GOTable = EstimatePVals(nullScores,realScores,whatTail,GOTable)
% EstimatePVals     Estimates p-values for each GO category from a given null distribution
%
% Computes using both a permutation-based and fitted Gaussian approximation.
%
% INPUTS:
% - nullScores (cell): a set of null scores for each GO category
% - realScores (vector): scores for the same set of GO categories obtained from real data
% - whatTail ('right' or 'left'): whether to use a right- or left-tailed test.
%           'right': tests whether each realScores is greater than the null.
%           'left': tests whether each realScores is less than the null.
% - GOTable (table): the table to save the computations back to.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Checks on matching nullScores, realScores, GOTable:
%-------------------------------------------------------------------------------
numCategories = height(GOTable);
assert(numCategories==length(nullScores))
assert(numCategories==length(realScores))

%-------------------------------------------------------------------------------
% Estimate category p-values, looping over categories
%-------------------------------------------------------------------------------
pValPerm = zeros(numCategories,1);
pValZ = zeros(numCategories,1);
for i = 1:numCategories
    scoreHere = realScores(i);
    nullHere = nullScores{i};
    switch whatTail
    case 'right'
        pValPerm(i) = mean(nullHere >= scoreHere);
        pValZ(i) = 1 - normcdf(scoreHere,mean(nullHere),std(nullHere));
    case 'left'
        pValPerm(i) = mean(nullHere <= scoreHere);
        pValZ(i) = normcdf(scoreHere,mean(nullHere),std(nullHere));
    otherwise
        error('Unknown tail setting, ''%s''',whatTail)
    end
end

%-------------------------------------------------------------------------------
% Multiple hypothesis correction using Benjamini-Hochberg (FDR):
%-------------------------------------------------------------------------------
pValPermCorr = mafdr(pValPerm,'BHFDR',true,'showPlot',false);
pValZCorr = mafdr(pValZ,'BHFDR',true,'showPlot',false);

%-------------------------------------------------------------------------------
% Assign values to categories of GOTable:
%-------------------------------------------------------------------------------
GOTable.pValZ = pValZ;
GOTable.pValZCorr = pValZCorr;
GOTable.pValPerm = pValPerm;
GOTable.pValPermCorr = pValPermCorr;

end
