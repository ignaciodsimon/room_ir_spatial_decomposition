function windowContents = getTimeWindowContent(inputIRs, centralWindowTime, windowSize)

    % Obtain the start and end markers
    initPoint = round(centralWindowTime - windowSize/2);
    if initPoint < 1
        initPoint = 1;
        disp('Warning: Limiting the beginning of the time window!')
    end
    endPoint = round(centralWindowTime + windowSize/2);
    if endPoint > size(inputIRs, 2);
        endPoint = size(inputIRs, 2);
        disp('Warning: Limiting the ending of the time window!')
    end

    % Return the extracted data
    windowContents = inputIRs(:, initPoint : endPoint);

end
