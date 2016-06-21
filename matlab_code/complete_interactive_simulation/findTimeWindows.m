function timeWindowsCenterTime = findTimeWindows(impulseResponses)

    addpath('../ir_decomposition/');

    % Low-pass filter and detect peaks on IRs
    filteredIRs = averageIRs(impulseResponses);
    filteredPeaks = findPeaksOnIrs(filteredIRs);

    % Central times for all detected windows
    timeWindowsCenterTime = filteredPeaks(:, 1);

end
