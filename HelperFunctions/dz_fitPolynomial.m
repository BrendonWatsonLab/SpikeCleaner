function bestFittedValues = dz_fitPolynomial(bestwf, timeVector)
    % Reduce standard deviation for Gaussian smoothing to prevent over-smoothing
    sigma = 2; % Lower value preserves peaks while smoothing
    kernelSize = 10; % Window size for smoothing

    %% Gaussian kernel
    x = linspace(-kernelSize/2, kernelSize/2, kernelSize);
    gaussKernel = exp(-x.^2 / (2 * sigma^2));
    gaussKernel = gaussKernel / sum(gaussKernel); % Normalize

    %% Gaussian smoothing
    smoothedValues = conv(bestwf, gaussKernel, 'same');

    %% Moving Average for extra smoothing
    windowSize = 5; % Adjust for more/less smoothing
    bestFittedValues = movmean(smoothedValues, windowSize);

    %% Plotting
    % figure;
    % plot(timeVector, bestwf, 'b--', 'LineWidth', 3);
    % hold on;
    % plot(timeVector, bestFittedValues, 'r', 'LineWidth', 2);
    % title('Gaussian + Moving Average Smoothed Waveform');
    % xlabel('Time (s)');
    % ylabel('Amplitude');
    % legend('Original Waveform', 'Gaussian + Moving Average Smoothed Fit');
    % hold off;
end
