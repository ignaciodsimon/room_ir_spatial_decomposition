function filteredIRs = averageIRShortly(inputData)

    filterLength = 5;
    filterTimes = 3;

    % Average all sensors inputs
    filteredIRs = abs(sum(inputData, 1));

    % Low pass filter the averaged signal
    for index = 1 : filterTimes
        filteredIRs = conv(filteredIRs, ones(filterLength, 1) / filterLength);
    end
    filteredIRs = filteredIRs(round(filterTimes * filterLength / 2) : length(filteredIRs));

end
