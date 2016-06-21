function [sortedFitnessEvaluation, sortedIndividualPool]  = rankAndSort(fitnessEvaluation, individualPool)
%This function will rank and sort the individual pool in ascending order
%based on the fitness evaluation of the individual pool

sortedIndividualPool = zeros(size(individualPool));

[sortedFitnessEvaluation, indexMatrix] = sort(fitnessEvaluation,1,'ascend');


%Now rearrange the individualPool
for i = 1 : size(fitnessEvaluation,1)
    sortedIndividualPool(i,:) = individualPool(indexMatrix(i),:);
end


end