function parentIDs = PropagateUp(idHere,hierarchyMatrix)
    % start at an ID and propagate up the hierarchy, annotating matches on the way
    % Get parents:
    beVocal = false;
    parentsHere = find(hierarchyMatrix(:,idHere));
    numParents = length(parentsHere);
    parentIDs = parentsHere;
    % if beVocal
    %     fprintf(1,'--%s--\n',GOTable.GOName{idHere});
    %     fprintf(1,'%u parents:',numParents);
    %     for j = 1:numParents
    %         fprintf(1,'%s, ',GOTable.GOName{parentsHere(j)});
    %     end
    %     fprintf(1,'\n');
    % end
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
            parentIDs = [parentIDs;PropagateUp(parentsHere(j),hierarchyMatrix)];
        end
    end
    parentIDs = unique(parentIDs);
end
