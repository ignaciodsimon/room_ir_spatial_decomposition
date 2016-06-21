function filteredPeaks = findPeaksOnIrs(inputIRs)

    % Find the peaks in the average
    [peaksValues, peaksPlaces] = findpeaks(inputIRs);
    peaksValues = peaksValues / max(peaksValues);

    % TODO: This has to come from some analysis of the peaks / signals ...
    peakThreshold = 10^(-25/20);

    % Filter the peaks according to their amplitude
    filteredPeaks = [];
    peakVector = zeros(size(inputIRs));
    for index = 1 : length(peaksPlaces)
        if peaksValues(index) > peakThreshold
            peakVector(peaksPlaces(index)) = peaksValues(index);
            filteredPeaks(size(filteredPeaks, 1) + 1, :) = [peaksPlaces(index) peaksValues(index)];
        end
    end

end