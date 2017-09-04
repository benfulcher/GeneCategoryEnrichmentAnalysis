function entrezID = GiveMeEntrezID(geneAbbreviation)
% Uses http://mygene.info to retrieve an entrezID from a given gene name
%-------------------------------------------------------------------------------

queryText = sprintf('http://mygene.info/v3/query?q=%s&fields=entrezgene&species=mouse&entrezonly=1',geneAbbreviation);
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
