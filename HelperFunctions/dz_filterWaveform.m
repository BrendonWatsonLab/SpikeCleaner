function dz_filterWaveform(inputFile, outputFile, samplingRate)
    % Function to  filter waveforms:
    % Parameters:
    % inputFile      - Input .mat file containing waveforms ('wf')
    % outputFile     - Output .mat file to save results
    % samplingRate   - Sampling frequency (in Hz)
  
    % Load waveform data
    if exist(inputFile, 'file') == 2
        load(inputFile, 'wf');
    else
        error('Input file %s does not exist.', inputFile);
    end

    % Initialize storage for results
    clippedWaveforms = cell(length(wf), 1); %highpass
    meanWaveforms = cell(length(wf), 1);
    bestWaveforms = cell(length(wf), 1);
    channelCorrelations = cell(length(wf), 1);
    % Smoothing parameters (Moving Average Filter)
    windowSize =9;
  
    %% Implementation in phy: highpass at 150Hz
    % Cutoff Frequency (in Hz)
    fc = 150; 
    % Normalized Cutoff Frequency (between 0 and 1)
    Wn = fc / (samplingRate / 2); 
    %Wn = (fc / samplingRate) * 2; 
    order=3;
    % Design High-Pass Butterworth Filter
    [b, a] = butter(order, Wn, 'high');


    % Loop through each cluster and process waveforms
    for ix = 1:length(wf)
        % Check if the waveform for this cluster is empty
        if isempty(wf{ix})
            clippedWaveforms{ix}=[];
            meanWaveforms{ix} = [];
            bestWaveformsChannel{ix} = [];
            bestWaveforms{ix} = [];
            channelCorrelations{ix} = [];
            continue;
        end

        %%  Extract waveforms for the current cluster
        waveforms = wf{ix};
        waveforms1=wf{ix}; %saving unfiltered
        [numSamples, numChannels] = size(waveforms);
       

        %% High-Pass Filter to remove high-frequency noise
        for ch = 1:numChannels
            waveforms(:, ch) = filtfilt(b, a,  waveforms(:, ch)); %iir
        end
        clippedWaveforms{ix} = waveforms;%just highpass
        
        waveforms1=waveforms; 
   
        %% MOVING AVERAGE WINDOW  :             
        for ch = 1:numChannels
            waveforms1(:, ch) = movmean(waveforms1(:, ch), windowSize);
        end
        meanWaveforms{ix} = waveforms1;%filtered highpass + movemean  smoothening

        %% Find the channel with the highest amplitude waveform
        amplitudes = max(waveforms, [], 1) - min(waveforms, [], 1);
        [~, maxCh] = max(amplitudes);
        bestWaveformsChannel{ix} = maxCh;  % best  wf  channel

        %% Calculate correlations between the best waveform and all other channels
        bestWaveform1 = waveforms(:, maxCh);
        bestWaveforms{ix}=bestWaveform1;

        % bestchannelforcorr=waveforms1(:, maxCh);  %%using mean  wfs -movmean
        % corrValues = zeros(1, numChannels);
        % for ch = 1:numChannels
        %     corrValues(ch) = corr(bestchannelforcorr, waveforms1(:, ch));
        % end
        % channelCorrelations{ix} = corrValues;


        % Define neighborhood:  instead of just  the max channel use max+/-3 to find corr with rest of the channels
        % store only the highestcorr from that
        neigh = 3;
        refChannels = maxCh + (-neigh:neigh);
        refChannels(refChannels < 1 | refChannels > numChannels) = []; % keep valid
        correlation=cell(length(refChannels),1);
        for i=1:numel(refChannels)
            x=refChannels(i);
            refWf = waveforms1(:, x);  % 81x1
            corrValues = zeros(1, numChannels);     % 1x128
            corrValues = corr(refWf, waveforms1);
            corrValues(:,x) = 0; %zeroing out the self correlation
            correlation{i} = corrValues;
        end    
        
        corre=cell2mat(correlation);
        cors=max(corre);
        % Sum of correlations per ref channel
        rowSums = sum(corre, 2);         % 1 sum per row/ref channel
        
        % Find the index of the ref channel with highest total correlation
        [maxScore, bestRefIdx] = max(rowSums);     % bestRefIdx = row in corrMa
        %full matrix
        channelCorrelations{ix} = cors;   % max channel amoungst the max+/-3chans

    end

    results.clippedWaveforms=clippedWaveforms; %high
    results.meanWaveforms = meanWaveforms; %filtered highpass + movemean  smoothening
    results.bestWaveformsChannel = bestWaveformsChannel;
    results.bestWaveforms = bestWaveforms;
    results.channelCorrelations = channelCorrelations;
 
    % Save all results to the output file
    save(outputFile, '-struct', 'results');
    disp(['Clipped waveforms and analysis results saved to ', outputFile]);
end
