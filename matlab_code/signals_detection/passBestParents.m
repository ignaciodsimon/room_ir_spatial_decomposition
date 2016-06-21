function poolOfParentsAndChildren = passBestParents(sortedIndividualPool, currentBestPool)

% Pass 2 parents through to next generation

passedParent = 1;

    for i = 1 : passedParent
        currentBestPool(size(currentBestPool,1)/2*i,:) = sortedIndividualPool(i,:);
    end
    
    poolOfParentsAndChildren = currentBestPool;
end