function [wf] = dz_getWaveform(datfil, clu, ts, clustersToCheck,parameters)
    disp("Entered into sm_getWaveform");
    disp(['Size of clu: ', num2str(size(clu))]);
    disp(['Size of ts: ', num2str(size(ts))]);
    disp(['Size of clusters to analyze: ', num2str(size(clustersToCheck))]);
    %matFilename = outputFileName;
       
    nChannels=parameters.nChannels;
    samplingRate=parameters.samplingRate;
   
    fid = fopen(datfil, 'r');
    if fid == -1
        error('Failed to open the .dat file.');
    end
  

   
    % Parameters for waveform extraction
    
    uclu=unique(clu);
    nClusters=length(uclu);% Number of clusters

    %% Initializations
    wf = cell(nClusters, 1);
    % Parameters for waveform extraction
    preTimeWindow =2e-3;  %2 ms before spike (in seconds)
    postTimeWindow =2e-3;   %2 ms after spike (in seconds)
    
    % Convert time window to sample indices
    nSamplesBefore = round(preTimeWindow * samplingRate);  % Convert 2 ms to samples
    nSamplesAfter = round(postTimeWindow * samplingRate);   % Convert 2 ms to samples
    totalSamples = nSamplesBefore+nSamplesAfter+1;

   
    % Loop through each cluster to extract waveforms
    for i = 1:length(clustersToCheck)
        actualclusterid = clustersToCheck(i);  % Get the current cluster ID
        %spikeIndices = find(clu == uclu(actualclusterid));  % Find spike indices for this cluster id at indices in uclu
        
        %%only changing it for plotting split
        spikeIndicesAll = find(clu == uclu(actualclusterid)); 
        
        %limit only if large
        if length(spikeIndicesAll) > 10000  
            % Randomly selection of the spikes
            nUse = round(length(spikeIndicesAll) / 10);
            spikeIndices = spikeIndicesAll(randperm(length(spikeIndicesAll), nUse));
        else
            spikeIndices = spikeIndicesAll;
        end

        if isempty(spikeIndices)
            disp(['Skipping cluster ', num2str(actualclusterid), ' as there are no spikes in this cluster']);
            continue;
        end  
        
        disp(['Processing Cluster ', num2str(actualclusterid)]);
        
        % Preallocate waveforms array to store waveforms for all spikes and all channels
        waveforms = zeros(totalSamples, nChannels, length(spikeIndices));  % 3D array: samples x channels x spikes
       
        %% Extracting waveforms for each spike in this cluster
        for j = 1:length(spikeIndices)
            fprintf('Processing Spike %d/%d in cluster %d\n', j, length(spikeIndices), actualclusterid);
            idx = spikeIndices(j);  % Current spike index
    
            %%add start and end time stamps and find those timestamps by sample index in dat file,basically convert timestamps to index with respect to dat file
            spikeTimeInSecs = ts(idx);  % Get the spike time in seconds
   
            spikeTimeIdx = round(spikeTimeInSecs * samplingRate);

            % Define start and end indices for extracting waveform data around the spike
            startIdx = max(1, spikeTimeIdx - nSamplesBefore);  % Prevent out of bounds
            endIdx = min(spikeTimeIdx + nSamplesAfter, spikeTimeIdx + nSamplesAfter);  % Prevent out of bounds
            numSamplesToRead = endIdx - startIdx + 1;

            % % Move the file pointer to the start of the desired section
            fseek(fid, (startIdx - 1) * nChannels * 2, 'bof');
            rawBlock = fread(fid, [nChannels, numSamplesToRead], 'int16')';
            extractedData = rawBlock;

            
            % Apply padding to keep the size consistent
            extractedSize = size(extractedData, 1);  % Size of the extracted data
            
            if extractedSize < totalSamples 
                % Pad the data if the extracted size is smaller than expected
                paddedData = zeros(totalSamples, nChannels);
                if startIdx == 1  % Spike near the beginning of the file, pad at the end
                    paddedData(1:extractedSize, :) = extractedData;
                else  % Spike near the end of the file, pad at the beginning
                    paddedData(end-extractedSize+1:end, :) = extractedData;
                end
                waveforms(:, :, j) = paddedData;  % Store the padded waveform for this spike
            else
                waveforms(:, :, j) = extractedData;  % Store the extracted waveform for this spike
            end
        end
        %%take the mean of the waveforms across all spikes
        meanWaveform = mean(waveforms, 3); %Time samples*nChannels
        wf{actualclusterid} = meanWaveform;

    end    
end    
