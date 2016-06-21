function parentPool = EMSelection(sortedIndividualPool)
% This function executes a parent selection routine on the sortedIndividualPool.
% This will form the parentPool. 

EMS = 1;


parentPool = zeros(size(sortedIndividualPool,1) / 2 + 1, size(sortedIndividualPool,2));

    if EMS
        FITTEST_PARENTS = 2;

        % Choose randomly between adjacent pairs
        k = 0;
        for i = 1 : size(parentPool,1) - 1

            choosePair = randi([1,2],1);
            parentPool(i + 1,:) = sortedIndividualPool(k + choosePair,:);
            k = k + 2;

        end


        % Make sure the #FITTEST_PARENTS is moved to parentPool
        parentPool(1:FITTEST_PARENTS,:) = sortedIndividualPool(1:FITTEST_PARENTS,:);
%         disp(['Fittest parent: ' num2str(sortedIndividualPool(1,:))])
%         disp(['Fittest parents: ' num2str(sortedIndividualPool(2,:))])
    end




end 