function theBestPool = pickTheBestPool(childrenPoolMutated, copiedPoolFull)
%This forms a pool that is a combination of 2 sub pools
%     amountOfBestForNextGen = 2;
    
    theBestPool = zeros(size(childrenPoolMutated));

    halfWayPoint = size(theBestPool,1)/2;

    bestOfChildrenPoolMutated = childrenPoolMutated(1:ceil(size(childrenPoolMutated,1)/2),:);

    bestOfCopiedPoolFull = copiedPoolFull(1:ceil(size(copiedPoolFull,1)/2),:);

    theBestPool = [bestOfChildrenPoolMutated; bestOfCopiedPoolFull];
    
%     theBestPool = copiedPoolFull(1:amountOfBestForNextGen,:);
                                        
end