function noisedCopy = bestFitNoised(copiedPoolFull, sortedIndividualPool)
% rng('shuffle');
amountOfIndividuals = 1;
bestFit = sortedIndividualPool(1:amountOfIndividuals,:);
    for individualIndex = 1 : amountOfIndividuals
        for geneIndex = 1 : size(sortedIndividualPool,2)
            if randi([1 2]) == 1
                signOfChange = 1;
            else
                signOfChange = -1;
            end

            amountOfChange = rand(1)/20 * pi / 180;

            bestFit(individualIndex, geneIndex) = bestFit(individualIndex, geneIndex) + amountOfChange * signOfChange;
        end
    end
    copiedPoolFull(4: 3 + amountOfIndividuals,:) = bestFit;
    noisedCopy = copiedPoolFull;
    
end