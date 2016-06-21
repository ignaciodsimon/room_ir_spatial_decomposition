function childrenPoolMutated =  mutation(childrenPool)
%This function will mutate a random gene in a random individual
PROBABILITY_MUTATION = 0.05;
amountOfChange = pi * 0.5/180;
childrenPoolMutated = childrenPool;
 

    % -- Mutate the childrenPool
    for i = 1 : size(childrenPool,1)
        if rand(1) < PROBABILITY_MUTATION
%             individualChooser = randi(size(childrenPool,1),1);
            individualChooser = i;
%             if individualChooser == 1
%                 disp('Emperor mutated')
%             end
%             geneChooser = randi(size(childrenPool,2),1);
            geneChooser = randi(size(childrenPool,2),1);

            childrenPoolMutated(individualChooser,geneChooser) = childrenPool(individualChooser,geneChooser) + ...
                                                                    (-amountOfChange + (amountOfChange+amountOfChange).*rand(1));
        else
        childrenPoolMutated = childrenPool;
        end
    end

    % -- Perform search space control
    for gene = 1 : size(childrenPoolMutated,2)
        for individual = 1 : size(childrenPoolMutated,1)
            if max(childrenPoolMutated(individual, gene)) > pi
                childrenPoolMutated(individual, gene) = pi;
%                 disp('over')
            elseif max(childrenPoolMutated(individual, gene)) < -pi
                childrenPoolMutated(individual, gene) = -pi;
%                 disp('under')
            end
        end
    end
    

%     childrenPoolMutated(1,:) = childrenPool(1,:);
    
end