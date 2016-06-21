function detetedAmount = detectAmountOfSignals(inputSignals)

    % Compute a low pass average of all IRs
    windowPeaksCount = zeros(size(inputSignals, 1), 1);
    for index = 1 : size(inputSignals, 1)

        % Estimate the slope
        filteredIR = averageIRShortly(inputSignals(index,:));

        % Normalize and vertical stretch
        filteredIR = (filteredIR - min(filteredIR)) / max(filteredIR - min(filteredIR));

        % Find the peaks in the average
        [tempIRPeaksValues, tempIRPeaksPlaces] = findpeaks(filteredIR);
        tempIRPeaksValues = tempIRPeaksValues / max(tempIRPeaksValues);

        tempPeaksCount = 0;
        for index2 = 1 : length(tempIRPeaksPlaces)
            if tempIRPeaksValues(index2) > 0.02
                tempPeaksCount = tempPeaksCount + 1;
            end
        end
        windowPeaksCount(index) = tempPeaksCount;      
    end

    % Return the value that appears most
    detetedAmount = mode(windowPeaksCount);
end

function filteredIRs = averageIRShortly(inputSignal)

    filterLength = 5;
    filterTimes = 3;

    % Average all sensors inputs
    filteredIRs = abs(sum(inputSignal, 1));

    % Low pass filter the averaged signal
    for index = 1 : filterTimes
        filteredIRs = conv(filteredIRs, ones(filterLength, 1) / filterLength);
    end
    filteredIRs = filteredIRs(round(filterTimes * filterLength / 2) : length(filteredIRs));

end
