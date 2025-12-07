function [cch,binCenters] = dz_autoCorr(t1, t2, binSize, numBins)
    % Ensure unique spike times to avoid duplicates
    t1 = unique(t1);  
    t2 = unique(t2);  

    % Define bin edges (ensuring symmetry)numBins
    total=(2*numBins)+1;
    half=(total-1)/2;
    %edges = linspace(-numBins * binSize, numBins * binSize, 2*numBins +2); 
    % use this if that  doesn't work

    % % Ensure edges are correctly calculated
    edges = (-half:half + 1) * binSize; 
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;  % Ensure proper centering
    centerShift = binCenters(round(end/2));  % Find what MATLAB thinks is the center
    binCenters = binCenters - abs(centerShift);   % Correct by shifting all bins

    % Initialize ACG
    cch = zeros(1, length(binCenters));

    % Iterate over spikes in t1 to compute ACG
    for i = 1:length(t1)
        diffs = t2(t2 ~= t1(i)) - t1(i);  % Exclude self-comparison
        cch = cch + histcounts(diffs, edges);
    end

    % Find the true center bin dynamically
    [~, centerIdx] = min(abs(binCenters)); % Find index closest to 0
    centerBinValue = binCenters(centerIdx); % Actual center bin value
    binCenters = binCenters - centerBinValue;  

    %%  Plotting
    % figure;
    % stairs(binCenters * 1000, cch, 'b', 'LineWidth', 1.5); % Convert to ms
    % hold on;
    % 
    % % Calculate Refractory Period Boundaries dynamically based on center bin
    % refractoryPeriod = 2e-3; % 2 ms
    % % xrp0 = centerBinValue - 2e-3;  % Lower boundary from actual center
    % % xrp1 = centerBinValue + 2e-3;  % Upper boundary from actual center
    % 
    % 
    % xrp0 =  - refractoryPeriod;  % Lower boundary from actual center
    % xrp1 =  refractoryPeriod;  % Upper boundary from actual center
    % 
    % ylimVal = max(cch) * 1.2;
    % plot([xrp0 * 1000, xrp0 * 1000], [0 ylimVal], 'r--', 'LineWidth', 1.2);
    % plot([xrp1 * 1000, xrp1 * 1000], [0 ylimVal], 'r--', 'LineWidth', 1.2);
    % 
    % % Labels
    % title('Autocorrelogram');
    % xlabel('Lag (ms)');
    % ylabel('Normalized Spike Count');
    % hold off;
end
