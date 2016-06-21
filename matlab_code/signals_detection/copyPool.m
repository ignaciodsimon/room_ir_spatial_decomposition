function [poolM2] = copyPool(sortedIndividualPool)
%This function will create two pools that consists of the best individuals
%before crossover is applied and just a copy of the individualPool,
%respectivily. These two pools should then be mutated a lot to escape from
%potential local minima:
amountOfChange = zeros(1,size(sortedIndividualPool,2));
PROBABILITY_MUTATION = 0.6;


% -- Form the pools
% sizeOfPoolM1 = [size(sortedIndividualPool,1) - randi([1 size(sortedIndividualPool,1)]), size(sortedIndividualPool,2)];
% poolM1 = sortedIndividualPool(sizeOfPoolM1(2),sizeOfPoolM1(1));
poolM2 = sortedIndividualPool;

    % -- Mutate the pools
    for geneIndex = 1 : size(poolM2,2)
        for individualIndex = 1 : size(poolM2,1)
            if rand(1) < PROBABILITY_MUTATION
                
                for sourceIndex = 1 : size(sortedIndividualPool,2)
                    amountOfChange(sourceIndex) = 1/(1 + var(sortedIndividualPool(:,sourceIndex))) * pi;
                    % disp(['Amount of change is:..' num2str(amountOfChange)])
                    % disp(['Variance is:..........' num2str(var(sortedIndividualPool(:,sourceIndex)))])
                end

    %             individualChooserPoolM1 = randi(size(poolM1,1),1);
    %             geneChooserPoolM1 = randi(size(poolM1,2),1);

%                 individualChooserPoolM2 = randi(size(poolM2,1),1);
%                 geneChooserPoolM2 = randi(size(poolM2,2),1);

    %             poolM1(individualChooserPoolM1,geneChooserPoolM1) = sortedIndividualPool(individualChooserPoolM1,geneChooserPoolM1) + ...
    %                                                                     (-amountOfChange + (amountOfChange+amountOfChange).*rand(1));

%                 poolM2(individualChooserPoolM2,geneChooserPoolM2) = sortedIndividualPool(individualChooserPoolM2,geneChooserPoolM2) + ...
%                                                                         (-amountOfChange + (amountOfChange+amountOfChange).*rand(1));
%             poolM2(individualIndex,geneIndex) = sortedIndividualPool(individualIndex,geneIndex) + ...
%                                                                         (-amountOfChange(geneIndex) + (amountOfChange(geneIndex) ...
%                                                                         +amountOfChange(geneIndex)).*rand(1));
                if randi([1 2]) == 1
                    signOfChange = 1;
                else
                    signOfChange = -1;
                end
            
            poolM2(individualIndex,geneIndex) = sortedIndividualPool(individualIndex,geneIndex) + amountOfChange(geneIndex) * signOfChange;

            else
    %         poolM1 = poolM1;
            poolM2(individualIndex,geneIndex) = poolM2(individualIndex,geneIndex);
            end
        end
    end
    
%         % -- Perform search space control
    for gene = 1 : size(poolM2,2)
        for individual = 1 : size(poolM2,1)
            if max(poolM2(individual, gene)) > pi
                poolM2(individual, gene) = pi;
%                 disp('over')
            elseif max(poolM2(individual, gene)) < -pi
                poolM2(individual, gene) = -pi;
%                 disp('under')
            end
        end
    end


end