function childrenPool = crossOver(parentPool)
%This function mixes the genes of the individual according the the cross
%over method
EMS = 0;
BMW = 1;

childrenPool = zeros(size(parentPool,1) * 2 - 2, size(parentPool,2));

%Perform crossover

CROSS_OVER_POINT = ceil(size(parentPool,2)/2);

    if EMS
    %It is always the emperor individual who cross over with lesser individuals

    emperor = parentPool(1,:);

    k = 1;

        for i = 1 : size(parentPool,1)/2

            if size(parentPool,2) < 2
                child1 = emperor;
                child2 = parentPool(i+1);
            else
                child1 = [emperor(1:CROSS_OVER_POINT) parentPool(i + 1, CROSS_OVER_POINT + 1:size(parentPool,2))];
                child2 = [parentPool(i + 1, 1:CROSS_OVER_POINT) emperor(CROSS_OVER_POINT + 1:size(emperor,2))];
            end

            childrenPool(k,:) = child1;
            childrenPool(k + 1,:) = child2;


            k = k + 2;

        end
    end
    
    if BMW
        
        k = 1;
        fatCounter = 0;

        for individualIndex = 1 : size(parentPool,1)/2
            currentParentFit = parentPool(individualIndex,:);
            currentParentFat = parentPool(size(parentPool,1) - fatCounter,:);
            
            if size(parentPool,2) < 2
                child1 = currentParentFit;
                child2 = currentParentFat;
            else
                child1 = [currentParentFit(1:CROSS_OVER_POINT) currentParentFat(CROSS_OVER_POINT + 1 : size(currentParentFat,2))];
                child2 = [currentParentFat(1:CROSS_OVER_POINT) currentParentFit(CROSS_OVER_POINT + 1 : size(currentParentFit,2))];
            end

            childrenPool(k,:) = child1;
            childrenPool(k + 1,:) = child2;

            k = k + 2;
            fatCounter = fatCounter + 1;
        end
    end
        

end