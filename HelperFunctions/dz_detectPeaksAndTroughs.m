function [allLocsSorted, allPointsSorted] = dz_detectPeaksAndTroughs(waveform,timeVector)
            [peaks, peakLocs] = findpeaks(waveform); % Peaks
            [troughs, troughLocs] = findpeaks(-waveform); % Troughs (negate to use findpeaks)
            troughs = -troughs; % Restoring original values

            %% Combine peaks and troughs
            allPoints = [peaks; troughs];
            allLocs = [peakLocs; troughLocs];
            allLocs_t = timeVector(allLocs);

            %% Sort locations to maintain order
            [allLocsSorted, sortIdx] = sort(allLocs_t);
            allPointsSorted = allPoints(sortIdx);
end