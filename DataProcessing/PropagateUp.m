function parentIDs = PropagateUp(idHere,hierarchyMatrix,beVocal,GOTerms,downInstead)
% Start at an ID and propagate up the GO hierarchy, annotating matches on the way
%-------------------------------------------------------------------------------

% Set defaults:
if nargin < 3
    % Can turn on talking to check performance
    beVocal = false;
end
if nargin < 4
    GOTerms = [];
end
if nargin < 5
    downInstead = false;
end

%-------------------------------------------------------------------------------
% Get direct parents:
if ~downInstead
    parentsHere = find(hierarchyMatrix(:,idHere));
else
    % (actually children!:)
    parentsHere = find(hierarchyMatrix(idHere,:));
end
numParents = length(parentsHere);
parentIDs = parentsHere;

%-------------------------------------------------------------------------------
if beVocal
    fprintf(1,'--%s--\n',GOTerms.GOName{idHere});
    fprintf(1,'%u parents:',numParents);
    for j = 1:numParents
        fprintf(1,'%s, ',GOTerms.GOName{parentsHere(j)});
    end
    fprintf(1,'\n');
end
if numParents == 0
    if beVocal
        fprintf(1,'Top of hierarchy at %u\n',idHere);
    end
else
    if beVocal
        fprintf(1,'Going up for %u parents\n',numParents);
    end
    for j = 1:numParents
        % Give parents intersection of their annotations and child's
        % geneEntrezAnnotations(parents(j)) = intersect(geneEntrezAnnotations(parents(j)),geneEntrezAnnotations(idHere));
        parentIDs = [parentIDs;PropagateUp(parentsHere(j),hierarchyMatrix,beVocal,GOTerms,downInstead)];
    end
end
parentIDs = unique(parentIDs);

end
