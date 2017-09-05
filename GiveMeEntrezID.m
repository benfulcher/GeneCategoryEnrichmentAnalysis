function entrezID = GiveMeEntrezID(geneAbbreviation,theSpecies)
% Uses http://mygene.info to retrieve an entrezID from a given gene name
%-------------------------------------------------------------------------------

% Check inputs, set defaults:
if nargin < 2
    theSpecies = 'mouse'; % 'mouse', 'human'
end
%-------------------------------------------------------------------------------

queryText = sprintf('http://mygene.info/v3/query?q=%s&fields=entrezgene&species=%s&entrezonly=1',geneAbbreviation,theSpecies);
data = webread(queryText);

switch length(data.hits)
case 1
    entrezID = data.hits.entrezgene;
case 2
    if iscell(data.hits)
        hasEntrez = cellfun(@(x)isfield(x,'entrezgene'),data.hits);
        entrezIDs = [data.hits{hasEntrez}.entrezgene];
    else
        % structure
        entrezIDs = [data.hits.entrezgene];
    end
    if length(entrezIDs)==1
        entrezID = entrezIDs;
    elseif all(entrezIDs==entrezIDs(1))
        entrezID = entrezIDs(1);
    else
        warning('%u inconsistent entrez IDs for %s',length(entrezIDs),geneAbbreviation);
    end
otherwise
    warning('Error getting entrez ID from mygene.info for %s',geneAbbreviation);
end


end
