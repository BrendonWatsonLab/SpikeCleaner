%%%By Diksha Zutshi
function dz_Curate(basename,datfil, clufile,thresholds)
        
    %%Algorithm flow:
    %1 Check for firing rate for non-biological clusters-noise
    %2 Check for noise
    %3 Check for cluster to be a MUA cluster
    %4 Assign as Single Unit

    
    disp('Inside assessNoise');
    disp(['Data file: ', datfil]);
    
    %% reading
    [acgEvaluationMode, minHW, minAmp,maxAmp, minSlope, firingThreshold,acgallthreshold,acgmaxthreshold,correlationthreshold]=thresholds{:};
    % Check if clufile is empty
    if isempty(clufile)
        error('clufile cannot be empty');
    end
    
    %% Extract the directory path from clufile
    [a, ~] = fileparts(clufile);
    %%add  folders to path
    %addpath([b filesep 'HelperFunctions']);%%helper function folder


    %% ACTIVE CHANENLS:  to only use the active channels
    % === Load rez to access channel map ===
    param=fullfile(a,'parameters.mat');
    load(param, 'parameters');

    nChannels = parameters.nChannels;
    fs=parameters.samplingRate;
    %parameters={totalChannels,fs};

    
    %% Load spike times
    tsfil = [a filesep 'spike_times.npy'];
    ts = readNPY(tsfil);
    ts = double(ts) / fs;
    %%reading clu file
    clu = readNPY(clufile);
    if length(clu) ~= length(ts)
        error('clu IDs do not match length of times');
    end


    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize variables
    originalClusters = unique(clu);
    uclu = originalClusters;
    nclu = length(uclu);  
    good = true(nclu, 1);
    noiseReasons = cell(nclu, 1);
    Reasons = cell(nclu, 1); %%%this will hold reasons why a cluster isn't Good.
    passReason=[];
    Half_width=[];
    Slope=[];
    percentHighCorr=[];
    amplitude=[];   
    corrCell = cell(nclu, 1);
    halfwidth=cell(nclu, 1);%%gives back all the halfwidths 
    slope=cell(nclu, 1);
    acgs = cell(nclu, 1);%to save all the acgs
    amplitudes = cell(nclu, 1);%to save all the acgs

   
    %% LOW FIRING RATE:
    %%%RECORDING SIZE
    % Calculate the recording duration
    fileInfo = dir(datfil); % Get file information
    fileSize = fileInfo.bytes; % Get file size in bytes
    totalSamples1 = fileSize / (2 * nChannels); % Calculate the total number of samples
    recordingDuration = totalSamples1 / fs; % Duration in seconds
    recordingDurationHours = recordingDuration / 3600; % Duration in hours     
    evaluateLowFiringRates(recordingDuration,firingThreshold)%low firing rate cluster  elimination
    
    %% NOISE CHECK: BASED ON WAVEFORMS


    good1=find(good);%%check  for all the rest of the clusters
    clustersToCheck=good1;
    
    %%EXTRACTION OF Waveforms    
    outputsDir = fullfile(a, 'Outputs');  % path to Outputs folder

    % Create Outputs folder if it doesn't exist
    if ~exist(outputsDir, 'dir')
        mkdir(outputsDir);
    end
    outputFileName = fullfile(outputsDir, basename + ".mat");
    outputfile     = fullfile(outputsDir, basename + "_filtered.mat");

 
    if ~exist(outputfile,'file')
        [wf] = dz_getWaveform(datfil, clu, ts, clustersToCheck,parameters);  
        save(outputFileName, 'wf');
        dz_filterWaveform(outputFileName, outputfile,fs);
        wfdata=load(outputfile); 
    else
        wfdata=load(outputfile); 
    end       
    %extracting the waveforms
    clippedWaveforms = wfdata.clippedWaveforms;
    meanWaveforms = wfdata.meanWaveforms; %smooothened wfs
    bestWaveformsChannel = wfdata.bestWaveformsChannel;
    bestWaveforms = wfdata.bestWaveforms;
    channelCorrelations = wfdata.channelCorrelations;
    wfs=length(clippedWaveforms);
     
    %%Noise check evaluation
    disp("Entering into Noise check");
    assessWaveformQuality(wfs, clippedWaveforms,meanWaveforms,bestWaveformsChannel, channelCorrelations,minHW,minAmp,minSlope)

    %% MUA CHECK
    good1=find(good);%%getting rest of the clusters
    clustersToCheck=good1;
    disp("entering into MUA check");
    
  
    acgfile = fullfile(outputsDir, basename + "acg.mat");
    assessAcgQuality(clustersToCheck,acgEvaluationMode,acgallthreshold,acgmaxthreshold,acgfile) %%---will give final decision of good clusters
     
    %% 
    clustersToRecluster=[];
    function evaluateLowFiringRates(recordingDuration,firingThreshold)   
        
        for ix = 1:nclu
            % actualclusterid=20;
            % tst = ts(clu == actualclusterid);
            i = uclu(ix);
            tst = ts(clu == i);

            totalSpikes=length(tst);
            Averagefiringrate=totalSpikes/recordingDuration;

            if(Averagefiringrate<firingThreshold)
                %disp("low");
                %disp(ix);
                noiseReasons{ix} = 'lowRate';
                Reasons{ix} = sprintf('Low firing rate: %.3f Hz < %.3f Hz (Threshold)', Averagefiringrate, firingThreshold);
                good(ix)=false;
                continue;
            else
                good(ix)=true;
            end    
        end 
    end
   
   %%  
    function assessWaveformQuality(wfs, clippedWaveforms,meanWaveforms,bestWaveformsChannel, channelCorrelations,minHW,minAmp,minSlope)

        for ix = 1:wfs
            
            spikingType = [];% Reset spikingType to an empty string
            %current=uclu(i); %%actual cluster id
            if isempty(clippedWaveforms{ix})    
                continue;      
            end 
                            
            actualclusterid=uclu(ix);

            %actualclusterid=100;
            %ix = find(uclu == actualclusterid); 
           %%check for waveform shape to label noise
            waveforms=clippedWaveforms{ix};
            bestCh = bestWaveformsChannel{ix};
            corrValues = channelCorrelations{ix};
            bestwf=waveforms(:,bestCh);
            smoothenedwf=meanWaveforms{ix};
            smoothenedbestwf=smoothenedwf(:,bestCh);
          

            channelCount = size(waveforms, 2);
            variance=cell(channelCount, 1);
            %%CONVERTING SAMPLES IN SECS
            timeVector = (0:length(bestwf)-1) / fs;
            sym=0;
            pass=0;
            

            %%%%Checking juts for best wf to begin with
            [bestFittedValues]=dz_fitPolynomial(smoothenedbestwf,timeVector); %%finding the best polynomial fit for the best wf-just highpassed
            % Find peaks and troughs
            [allLocsSorted, sortIdx] = dz_detectPeaksAndTroughs(bestFittedValues,timeVector); %%smoothened wf
            [allLocsSorted1, sortIdx1] = dz_detectPeaksAndTroughs(bestwf,timeVector);%% from actual wf
            
            % If allLocsSorted1 has fewer indices than allLocsSorted, migrate missing ones: sometimes some points aren't captured 
            while length(allLocsSorted1) ~= length(allLocsSorted)
                % Find missing indices (values present in allLocsSorted but not in allLocsSorted1)
                missingIdx = setdiff(allLocsSorted, allLocsSorted1, 'stable'); 
            
                if isempty(missingIdx)
                    break; % No missing indices, exit loop
                end
            
                % Ensure column format for safe concatenation
                allLocsSorted1 = allLocsSorted1(:);  
                missingIdx = missingIdx(:);
                sortIdx1 = sortIdx1(:);
            
                % Find closest matching indices in timeVector for missingIdx
                closestIdx = arrayfun(@(x) find(timeVector >= x, 1), missingIdx);
            
                % Retrieve amplitudes from bestwf based on the closest index in timeVector
                missingAmplitudes = bestwf(closestIdx);  
            
                % Append the migrated index to allLocsSorted1
                allLocsSorted1 = [allLocsSorted1; missingIdx];  
            
                % Append the corresponding amplitudes found from bestwf (not from sortIdx)
                sortIdx1 = [sortIdx1; missingAmplitudes];  
            end


            % Convert back to row vectors
            allLocsSorted1 = allLocsSorted1.';
            sortIdx1 = sortIdx1.';
            % Re-sort both allLocsSorted1 and sortIdx1 to maintain order
            [allLocsSorted1, newSortIdx] = sort(allLocsSorted1);
            sortIdx1 = sortIdx1(newSortIdx);  % Apply the same sorting to sortIdx1

            %%making sure consicutive indices are atleast 50 uV apart
            i = 1;
            while i < length(sortIdx1)
                % Compute the absolute amplitude difference between consecutive values
                if sign(sortIdx1(i)) == sign(sortIdx1(i+1)) && abs(sortIdx1(i) - sortIdx1(i+1)) < 50
                %if abs(sortIdx1(i) - sortIdx1(i+1)) < 50
                    % Compare absolute amplitudes to keep the stronger signal
                    if abs(sortIdx1(i)) < abs(sortIdx1(i+1))
                        % Remove the smaller one (i)
                        sortIdx1(i) = [];
                        allLocsSorted1(i) = [];
                    else
                        % Remove the smaller one (i+1)
                        sortIdx1(i+1) = [];
                        allLocsSorted1(i+1) = [];
                    end
                    % Do not increment i — recheck new pair at position i
                else
                    i = i + 1;  % Only advance when pair passes threshold
                end
            end

            
            
            %%undestanding if its a positive or regular spiking wf, then selecting wf points: 
            %seleting the wf more towards the left,in positive case
       
            try
                if length(sortIdx1) > 3
                    baselineThresh = 100;  % µV
                    baselinePeak = [];
                    baselineLoc = [];

                    % Step 1: Remove Consecutive Peaks (Keep the Higher One)
                    peakIndices = find(sortIdx1 > 0); % Find all peaks
                    baselineCandidates = peakIndices(abs(sortIdx1(peakIndices)) < baselineThresh);
                    if ~isempty(baselineCandidates)
                        % Find the index with the **minimum amplitude**
                        idx = baselineCandidates(1);
                        baselinePeak = sortIdx1(idx);
                        baselineLoc = allLocsSorted1(idx);
                
                    end

                    i = 1;
                    while i < length(peakIndices)
                        if peakIndices(i+1) - peakIndices(i) == 1  % Consecutive peaks detected
                            if sortIdx1(peakIndices(i)) < sortIdx1(peakIndices(i+1))
                                % Remove the lower peak
                                sortIdx1(peakIndices(i)) = [];
                                allLocsSorted1(peakIndices(i)) = [];
                            else
                                sortIdx1(peakIndices(i+1)) = [];
                                allLocsSorted1(peakIndices(i+1)) = [];
                            end
                            peakIndices = find(sortIdx1 > 0); % Recalculate peaks
                        else
                            i = i + 1;
                        end
                    end
            
                    % Step 2: Remove Consecutive Troughs (Keep the Most Negative One)
                    troughIndices = find(sortIdx1 < 0); % Find all troughs
                    i = 1;
                    while i < length(troughIndices)
                        if troughIndices(i+1) - troughIndices(i) == 1  % Consecutive troughs detected
                            if abs(sortIdx1(troughIndices(i))) < abs(sortIdx1(troughIndices(i+1)))
                                % Remove the less negative trough
                                sortIdx1(troughIndices(i)) = [];
                                allLocsSorted1(troughIndices(i)) = [];
                            else
                                sortIdx1(troughIndices(i+1)) = [];
                                allLocsSorted1(troughIndices(i+1)) = [];
                            end
                            troughIndices = find(sortIdx1 < 0); % Recalculate troughs
                        else
                            i = i + 1;
                        end
                    end
            
                    % Step 3: Check If It Has a Bigger Peak or a Deeper Trough
                    [maxPeak, peakIdx] = max(sortIdx1);  % Find the largest peak
                    [minTrough, troughIdx] = min(sortIdx1);  % Find the deepest trough
            
                    if abs(minTrough) > abs(maxPeak)
                        % More negative trough → Regular Spiking
                        spikingType = 'Regular Spiking';
            
                        % Ensure the deepest trough has peaks on both sides
                        if (troughIdx > 1 && troughIdx < length(sortIdx1))
                            if sortIdx1(troughIdx - 1) > 0 && sortIdx1(troughIdx) < 0 && sortIdx1(troughIdx + 1) > 0
                                selectedIndices = [troughIdx - 1, troughIdx, troughIdx + 1]; % Peak-Trough-Peak
                            end
                        else
                            sym = 1; % Noise – if even after being regular spiking it has less than 3 indices
                            %passReason="Waveform doesn't look physiological";
                        end
                    else
                        % More positive peak → Positive Spiking
                        %%trough before the peak in postive spiking is near to 0 or sometimes above it.
                        spikingType = 'Positive Spiking';
                        % Only do this in Positive Spiking condition before setting selectedIndices
                        if (~isempty(baselinePeak) && peakIdx==1) %%in case because consecutive peak removal, highest peak is at 1, in positive spiking
                            % Prepend baseline peak and loc so it appears before the main peak
                            sortIdx1(peakIdx+2)=[];
                            allLocsSorted1(peakIdx+2)=[];
                            sortIdx1 = [baselinePeak, sortIdx1]; 
                            allLocsSorted1 = [baselineLoc, allLocsSorted1];
                        
                            % Since we inserted a new point at the beginning, we shift peakIdx by +1
                            peakIdx = peakIdx + 1;
                        end
                        
                        % Now build your Trough-Peak-Trough triplet
                        selectedIndices = [peakIdx - 1, peakIdx, peakIdx + 1];
            
                        % Ensure the highest peak has troughs on both sides
                        if (peakIdx > 1 && peakIdx < length(sortIdx1))
                            if (sortIdx1(peakIdx -1) < 0)||(sortIdx1(peakIdx + 1) < 0)
                                selectedIndices = [peakIdx - 1, peakIdx, peakIdx + 1]; % Trough-Peak-Trough
                            end
                        else
                            sym = 1; % Noise – if even after being positive spiking it has less than 3 indices
                        end
                    end
            
                    % Step 4: Extract Only the Relevant Indices
                    if exist('selectedIndices', 'var')
                        allLocsSorted1 = allLocsSorted1(selectedIndices);
                        sortIdx1 = sortIdx1(selectedIndices);
                    end
                end
                %%%if amplitude is >maxAmp then its noise 
                if abs(bestwf)<maxAmp
                     %% find new best channel
                    if length(sortIdx1) == 3
                        [spikingType, quality,passReason1,halfWidth,Slope1,Slope2,amplitude1] = dz_analyzeSpikeType(bestwf, timeVector, allLocsSorted1, actualclusterid, bestCh,minHW,minAmp,minSlope);
                     
                        Half_width=halfWidth;%ms
                        Slope=(min(Slope1,Slope2))/100;% uV/s --->uV/ms
                        Slope=abs(Slope);

                        amplitude=amplitude1;%uV
                    else
                        passReason="Waveform doesn't look physiological";
                        pass=0;
                    end
                else
                    disp('Amplitude exceeds ±2000 µV — likely artifact');
                    passReason="Amplitude exceeds ±2000 µV — likely artifact";
                    pass=0;
                end  

            catch
                passReason="Waveform doesn't look physiological";
                pass=0;
          end


         %%%%%%%%%%%%%%correlation check for channel corr,if more than 50% channels ar highly corr then noise otherwise good
         
        corrValues = channelCorrelations{ix};
        corr=corrValues;
        totalChannels = length(corr); % Total number of remaining channels
        numHighCorr = sum(corr > correlationthreshold);
        percentHighCorr = (numHighCorr / totalChannels) * 100;


        if (percentHighCorr>80) 
                pass=0; %%noise
                %disp(['Cluster ID: ', num2str(actualclusterid),' | Percentage of High Correlation Channels: ', num2str(percentHighCorr), '%']);
                %disp('Noise');
                passReason = sprintf('High Correlation (%.2f%%) on all Channels', percentHighCorr);

        else
                pass=1; %GOOD
                %disp(['Cluster ID: ', num2str(actualclusterid),' | Percentage of High Correlation Channels: ', num2str(percentHighCorr), '%']);
                %disp('Good');
               
        end     
       

        if isempty(Half_width)
            Half_width = '';  % Assign empty string if halfWidth is empty
        end

        if isempty(Slope)
            Slope = '';
        end
        if isempty(percentHighCorr)
            percentHighCorr = '';
        end
        if isempty(amplitude)
            amplitude = '';
        end
    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%FINAL DECISIONS%%%%%%%%%%%%%%%%%%%%%%%%%
        if pass==0  %%for  if best wf doesnt have either peak or trough
            %disp('Noise');
            noiseReasons{ix} = 'Noise';
            Reasons{ix} = passReason;
            halfwidth{ix}=Half_width;
            corrCell{ix}=percentHighCorr;
            slope{ix}=Slope;
            good(ix)=false;
            amplitudes{ix}=amplitude;
    
        else
            %disp('Good');
            noiseReasons{ix} = 'Good';
            good(ix)=true;
            halfwidth{ix}=Half_width;
            slope{ix}=Slope;
            corrCell{ix}=percentHighCorr;
            Reasons{ix} = passReason;
            amplitudes{ix}=amplitude;
        end 
        end
     end 
                   

    %% CHECKING FOR MUA
     function assessAcgQuality(clustersToCheck,acgEvaluationMode, acgallthreshold,acgmaxthreshold,acgfile)
        if exist(acgfile,'file')
            load(acgfile,'acgs');
        end    
        for ix = 1:nclu
           
           if ~ismember(ix, clustersToCheck)
                continue;
           end
       
            
           actualindex=uclu(ix);
           %actualindex=11;
           %ix = find(uclu == actualindex);
           tst = ts(clu == actualindex);
           binSize= 0.001;%1MS

            if exist(acgfile,'file')
               acgData=acgs{ix};
               currentacg = acgData.currentacg;
               binCenters = acgData.binCenters;
            else
                numBins=25;
                [currentacg,binCenters] = dz_autoCorr(tst, tst, binSize, numBins);
                acgs{ix} = struct('currentacg', currentacg, 'binCenters', binCenters);
            end 
               
           
           pred = mean(currentacg(1:10));
           zeroLagIndex = ceil(length(currentacg) / 2); 
          
           zeroLagRange = (zeroLagIndex - 2):(zeroLagIndex + 2);
           cchWithoutZeroLag = currentacg;
           cchWithoutZeroLag(zeroLagRange) = []; % Remove zero lag values
           if(pred==0)
               pred=mean(cchWithoutZeroLag);
           end    
           % visualize acgs
           % % % 
           % figure;
           % stairs(binCenters * 1000, currentacg, 'b', 'LineWidth', 1.5); % Convert to ms
           % hold on;
           % % **Mark Refractory Period (±2ms)**
           % refractoryPeriod = 2; % 2 ms
           % xrp0 = -refractoryPeriod;
           % xrp1 = refractoryPeriod;
           % ylimVal = max(currentacg) * 1.2; % To adjust vertical limits
           % plot([xrp0, xrp0], [0 ylimVal], 'r--', 'LineWidth', 1.2);
           % plot([xrp1, xrp1], [0 ylimVal], 'r--', 'LineWidth', 1.2);
           % 
           % % Labels and Title
           % title(sprintf('Autocorrelogram for Cluster %d', actualindex));
           % xlabel('Lag (ms)');
           % ylabel('Spike Count');
           % legend('ACG', '2 ms Refractory');
           % hold off;

           
           cchzero = currentacg(zeroLagIndex-3:zeroLagIndex+2); % Zero lag INDEXES
           
            %%derivative based peaks in the shoulders
            % Find first peak after central refractory zone
           searchStart = zeroLagIndex + 3;  % skip 3 bins on right side of center
           seg = currentacg(searchStart:end);
          
            %%finding the higheest peak
            k =20;                                 

            if numel(seg) > k && all(seg(1) > seg(2:k+1))
                peakRel  = 1;
                peakVal  = seg(1);
            else
                [pks, locs] = findpeaks(seg, ...
                    'MinPeakProminence', max(5, 0.02*max(seg)), ...
                    'MinPeakDistance', 2);
            
                if isempty(pks)
                    [peakVal, peakRel] = max(seg);          % fallback
                else
                    maxPk   = max(pks);
                    iFirst  = find(pks == maxPk, 1, 'first'); % earliest among tallest
                    peakRel = locs(iFirst);
                    peakVal = pks(iFirst);
                end
            end
            
            % Convert to absolute index 
            peakIdxAbs = searchStart + peakRel - 1;

            %  3  values around max peak  and  avg it:       
            % Get surrounding peak values if available
            peakNeighborhood = [];            
            if peakRel > 1 && seg(peakRel-1) > 0
                peakNeighborhood(end+1) = seg(peakRel-1);
            end          
            peakNeighborhood(end+1) = seg(peakRel);  % always include main peak
            if peakRel < length(seg) && seg(peakRel+1) > 0
                peakNeighborhood(end+1) = seg(peakRel+1);
            end           
            % Final refined peak value
            avgPeak = mean(peakNeighborhood);
            proportion=cchzero/avgPeak;  %proportion of  center  peaks  in  comparison to the shoulder of the  acg
    
            %%block  to  check  for noise using:ACG>>>  if all  center peaks  are >100%  of shoulder peaks
            noiseflag=false;
            if all(proportion>acgallthreshold)   
                noiseflag=true;
                passReason = sprintf('Centerbins proportion greater than 80 percent of the shoulder (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                       proportion(1) * 100, proportion(2) * 100, proportion(4) * 100, proportion(5) * 100, proportion(6) * 100);  
            elseif any(proportion>acgmaxthreshold)
                noiseflag=true;
                passReason = sprintf('Atleast one Centerbin has proportion greater than 110 percent of the shoulder (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                       proportion(1) * 100, proportion(2) * 100, proportion(4) * 100, proportion(5) * 100, proportion(6) * 100);  
            else 
                noiseflag=false;
            end
           if noiseflag==true
                good(ix)=false;
                noiseReasons{ix} = 'Noise';
                Reasons{ix} = passReason;
                continue;  
           end  
     

           flag=false;       
           switch acgEvaluationMode
               case 'lenient'    
                    %±2 ms bins must be < 30
                    if proportion(1) > 0.3 || proportion(6) > 0.3
                        flag = true;
                        disp('MUA: violation at ±2 ms bins (>30%)');
                        passReason = sprintf('High activity in 2ms bins compared to the shoulders — %.2f%%, %.2f%%', proportion(1) * 100, proportion(6) * 100);
        

                    %Allow 0th and ±1ms bins to be < 30% if outer bins are as low as 30
                    elseif (proportion(1) < 0.3 && proportion(6) < 0.3) && all(proportion(2:5) < 0.2)
        
                        flag = false;
                        disp('Single Unit');
                        passReason = sprintf('Passed all the checks — (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                        proportion(1) * 100, proportion(2) * 100, proportion(4) * 100, proportion(5) * 100, proportion(6) * 100);   
    
                    else 
                        flag = true;
                        %disp('MUA');
                        passReason = sprintf('High activity in 0th & 1ms bins compared to the shoulders — (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                        proportion(1) * 100, proportion(2) * 100, proportion(4) * 100, proportion(5) * 100, proportion(6) * 100);

                    end


               case 'strict'    
                    if any((proportion*100)>10)
                       flag = true;
                       disp('MUA: violation at ±2 ms bins');
                       passReason = sprintf('MUA: violation at ±2 ms bins — (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                       proportion(1) * 100, proportion(2) * 100, proportion(4) * 100, proportion(5) * 100, proportion(6) * 100);           
                    else        
                       flag=false;
                       disp('Single Unit');
                       passReason = sprintf('Passed all the checks — (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                        proportion(1) * 100, proportion(2) * 100, proportion(4) * 100, proportion(5) * 100, proportion(6) * 100);     
                    end   
           end    
           % 
        
           if flag==true
                good(ix)=false;
                noiseReasons{ix} = 'MUA ';
                Reasons{ix} = passReason;
                continue;  
           end  
            %%labelling good if it passed all the noise checks
          
            good(ix)=true;
            noiseReasons{ix}='good';
            Reasons{ix} = passReason;
              
    
        end   
    
        if ~exist(acgfile,'file')
            save(acgfile,'acgs', '-v7.3');
        end 
     end    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% WRITING THE OUTPUT IN THE SAME FOLDER AS OF DATA TO BE READ BY PHY
    outputFile = fullfile(a, 'cluster_SpikeCleaner.tsv'); %%inside the animal  file
    
    fid = fopen(fullfile(a,'cluster_SpikeCleaner.tsv'),'wt');
    fprintf(fid,'cluster_id\tSpikeCleaner\r\n');     % header

    for ix = 1:numel(uclu)
        cluster_index = uclu(ix);
        s = noiseReasons{ix};                        % a single string/char
        s = regexprep(s,'\r?\n',' ');               % remove any embedded newlines
        fprintf(fid,'%d\t%s\r\n',cluster_index,s);  % CRLF so Excel/Windows behave
    end
    fclose(fid);
    
    %% WRITING THE OUTPUT IN THE SAME FOLDER AS OF DATA TO BE READ BY PHY
    outputFile = fullfile(a, 'cluster_Spikereasons.tsv'); %%inside the animal  file
    fid = fopen(fullfile(a,'cluster_Spikereasons.tsv'),'wt');
    fprintf(fid,'cluster_id\tReasons\r\n');

    for ix = 1:numel(uclu)
        cluster_index = uclu(ix);
        s = Reasons{ix};
        if isempty(s), s = ' '; end
        s = regexprep(s,'\r?\n',' ');
        fprintf(fid,'%d\t%s\r\n',cluster_index,s);  % two columns, new row each time
    end
    fclose(fid);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%WRITING HALFWIDTHS 
    halfwidthfile=fullfile(a,'halfwidths.tsv' );
    % % Open the file for writing halfwidth
    fid1 = fopen(halfwidthfile, 'wt');
    if fid1 == -1
        error('Cannot open file: %s', halfwidthfile);
    end 
    fprintf(fid1, 'cluster_id\tHalf Widths\n');

    for ix = 1:length(uclu)  % Loop over the actual cluster IDs (since uclu contains the cluster IDs)
        cluster_index=uclu(ix);
        halfWidth=halfwidth{ix};
        if isempty(halfWidth)
            halfWidth = '';  % Default to " " if no specific reason provided
        end
        fprintf(fid1, '%d\t%d\n', cluster_index,halfWidth);

    end
    fclose(fid1);
    disp('Finished processing halfwidth');
        
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%WRITING SLOPES 
    slopes=fullfile(a,'slopes.tsv' );
    % 
    % % Open the file for writing halfwidth
    fid1 = fopen(slopes, 'wt');
    if fid1 == -1
        error('Cannot open file: %s', slopes);
    end 
    fprintf(fid1, 'cluster_id\t slopes \n');

    for ix = 1:length(uclu)  % Loop over the actual cluster IDs (since uclu contains the cluster IDs)
        cluster_index=uclu(ix);
        Slope=slope{ix};
        if isempty(Slope)
            Slope = '';  % Default to " " if no specific reason provided
        end
        fprintf(fid1, '%d\t%d\n', cluster_index,Slope);

    end
    fclose(fid1);
    disp('Finished processing slopes');
    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%WRITING corr 
    correlation=fullfile(a,'correlation.tsv' );
    % 
    % % Open the file for writing halfwidth
    fid1 = fopen(correlation, 'wt');
    if fid1 == -1
        error('Cannot open file: %s', correlation);
    end 
    fprintf(fid1, 'cluster_id\t correlation \n');

    for ix = 1:length(uclu)  % Loop over the actual cluster IDs (since uclu contains the cluster IDs)
        cluster_index=uclu(ix);
        percentHighCorr=corrCell{ix};
        if isempty(percentHighCorr)
            percentHighCorr = '';  % Default to " " if no specific reason provided
        end
        fprintf(fid1, '%d\t%d\n', cluster_index,percentHighCorr);

    end
    fclose(fid1);
    disp('Finished processing corr');
    



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%WRITING amplitude 
    amplitudess=fullfile(a,'amplitudess.tsv' );
    % 
    % % Open the file for writing halfwidth
    fid1 = fopen(amplitudess, 'wt');
    if fid1 == -1
        error('Cannot open file: %s', amplitudess);
    end 
    fprintf(fid1, 'cluster_id\t amplitude \n');

    for ix = 1:length(uclu)  % Loop over the actual cluster IDs (since uclu contains the cluster IDs)
        cluster_index=uclu(ix);
     
        amplitude=amplitudes{ix};
        if isempty(amplitude)
            amplitude = '';  % Default to " " if no specific reason provided
        end
        fprintf(fid1, '%d\t%d\n', cluster_index,amplitude);

    end
    fclose(fid1);
    disp('Finished processing amplitude');
     
end
