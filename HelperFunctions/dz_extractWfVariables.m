%% extract waveform variables
function [amplitudes,halfwidths,slopes,spikeType,wfs,good,noiseReasons,Reasons,chRangeAmp]=dz_extractWfVariables(wfdata,uclu,fs,good,noiseReasons,Reasons,pos)  


    %extracting the waveforms::: in uV
    clippedWaveforms = wfdata.clippedWaveforms;
    meanWaveforms = wfdata.meanWaveforms; %smooothened wfs
    bestWaveformsChannel = wfdata.bestWaveformsChannel;
    bestWaveforms = wfdata.bestWaveforms;
    wfs=length(clippedWaveforms);
    nclu=numel(uclu);
    
    % initializing variables
    halfwidths=cell(nclu, 1);%%gives back all the halfwidths 
    slopes=cell(nclu, 1);
    amplitudes = cell(nclu, 1);%to save all the acgs
    spikeType=cell(nclu, 1);
    chRangeAmp=cell(nclu,1);
    maxDistanceUm=150;%microns
    
    for ix = 1:wfs
            
        spikingType = [];% Reset spikingType to an empty string
        Slope=[];
        peaktopeak=[];
        halfWidth=[]; 
        peaktopeakChRange=[];

        %current=uclu(ix); %%actual cluster id
        if isempty(clippedWaveforms{ix})    
            continue;      
        end 

        actualclusterid=uclu(ix);
% 
%         actualclusterid=21;            
%         ix = find(uclu == actualclusterid);

%             
       %%check for waveform shape to label noise
        waveforms=clippedWaveforms{ix};
        bestCh = bestWaveformsChannel{ix};
        bestwf=waveforms(:,bestCh);
        smoothenedwf=meanWaveforms{ix};
        smoothenedbestwf=smoothenedwf(:,bestCh);

        %%CONVERTING SAMPLES IN SECS
        timeVector = (0:length(bestwf)-1) / fs;
        
        %getting near by channels:
        bestPos = pos(bestCh,:);           
        d = sqrt(sum((pos - bestPos).^2,2));
        chRange = find(d <= maxDistanceUm);   
        chRangeWaveforms=waveforms(:,chRange);

        
 
        %%%%Checking just for best wf to begin with
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
                        % Noise – if even after being regular spiking it has less than 3 indices      
                        noiseReasons{ix} = 'Noise';
                        Reasons{ix} = sprintf('Waveform doesnt look physiological');
                        good(ix)=false;
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
                        % Noise – if even after being positive spiking it has less than 3 indices      
                        noiseReasons{ix} = 'Noise';
                        Reasons{ix} = sprintf('Waveform doesnt look physiological');
                        good(ix)=false;
                    end
                end

                % Step 4: Extract Only the Relevant Indices
                if exist('selectedIndices', 'var')
                    allLocsSorted1 = allLocsSorted1(selectedIndices);
                    sortIdx1 = sortIdx1(selectedIndices);
                end
%             else
%                   noiseReasons{ix} = 'Noise';
%                   Reasons{ix} = sprintf('Waveform doesnt look physiological');
%                   good(ix)=false;
%                 
            end
     
            if ~good(ix) % don't process further if waveform doesn't look physiological
                continue;
            else
               
               peaktopeak=max(bestwf)-min(bestwf);
               peaktopeakChRange=max(chRangeWaveforms)-min(chRangeWaveforms);

               %% for halfwidth and slope:::::::::

               % Convert allLocsSorted (time values) into indices for waveform lookup
              allLocsSorted_idx = arrayfun(@(t) find(timeVector >= t, 1), allLocsSorted1);   
              % Ensure indices stay within valid range
              allLocsSorted_idx = min(max(allLocsSorted_idx, 1), length(bestwf));
              centralPoint = bestwf(allLocsSorted_idx(2));
              centerIndex = allLocsSorted_idx(2); %2nd point is central for analysis


                if centralPoint > 0
                     spikingType = 'Positive Spiking';
                     %%lets determine height of  the peak at the center: amplitude
                     amplitude = centralPoint;
                     amp=bestwf(allLocsSorted_idx);
                     %amp=abs(amp);

                    if  abs(amp(1)> amp(3))
                        amplitude= amp(1)-amp(2);
                        amplitude=abs(amplitude);
                        ypoint=(amp(1)+amp(2))/2;
                        xpoint=(allLocsSorted_idx(1)+allLocsSorted_idx(2))/2;
                        greaterIndex=amp(1);

                    else
                        amplitude= amp(3)-amp(2);
                        amplitude=abs(amplitude);
                        ypoint=(amp(3)+amp(2))/2;
                        xpoint=(allLocsSorted_idx(3)+allLocsSorted_idx(2))/2;
                        greaterIndex=amp(3);

                    end   

                     %%half width = 0.5*amplitude: whats the index there: measure the width of the loop
                     % Calculate half-width at half amplitude
                     halfAmplitude = ypoint;
                     halfAmplitudeIndex=xpoint;


                     [~, leftIdx] = min(abs(bestwf(allLocsSorted_idx(1):allLocsSorted_idx(2)) - halfAmplitude));
                     leftIdx = leftIdx + allLocsSorted_idx(1) - 1;  % Adjust index to global waveform

                     [~, rightIdx] = min(abs(bestwf(allLocsSorted_idx(2):allLocsSorted_idx(3)) - halfAmplitude));
                     rightIdx = rightIdx + allLocsSorted_idx(2) - 1;  % Adjust index to global waveform

                     halfWidth = timeVector(rightIdx) - timeVector(leftIdx);

                     Slope1 =(bestwf(allLocsSorted_idx(2)) - bestwf(leftIdx)) / ...
                        (timeVector(allLocsSorted_idx(2)) - timeVector(leftIdx) + eps);

                     Slope2 =(bestwf(rightIdx) - bestwf(allLocsSorted_idx(2))) / ...
                         (timeVector(rightIdx) - timeVector(allLocsSorted_idx(2)) + eps);



                else
                    spikingType = 'Regular Spiking';

                    %%lets  determine bigger peak to use for amplitude
                    amp=bestwf(allLocsSorted_idx);
                    %amp=abs(amp);

                    if  amp(1)> amp(3)
                        amplitude= abs(amp(1)-amp(2));
                        ypoint=(amp(1)+amp(2))/2;
                        xpoint=(allLocsSorted_idx(1)+allLocsSorted_idx(2))/2;
                        greaterIndex=amp(1);

                    else
                        amplitude= abs(amp(3)-amp(2));
                        ypoint=(amp(3)+amp(2))/2;
                        xpoint=(allLocsSorted_idx(3)+allLocsSorted_idx(2))/2;
                        greaterIndex=amp(3);

                    end    
                    % Half-amplitude for width calculation
                    halfAmplitude = ypoint;
                    halfAmplitudeIndex=xpoint;
                    %find points index in the waveform where y=ypoint

                    [~, leftIdx] = min(abs(bestwf(allLocsSorted_idx(1):allLocsSorted_idx(2)) - halfAmplitude));
                    leftIdx = leftIdx + allLocsSorted_idx(1) - 1;  

                    [~, rightIdx] = min(abs(bestwf(allLocsSorted_idx(2):allLocsSorted_idx(3)) - halfAmplitude));
                    rightIdx = rightIdx + allLocsSorted_idx(2) - 1;  

                    halfWidth = timeVector(rightIdx) - timeVector(leftIdx);

                    Slope1 =(bestwf(allLocsSorted_idx(2)) - bestwf(leftIdx)) / ...
                        (timeVector(allLocsSorted_idx(2)) - timeVector(leftIdx) + eps);

                    Slope2 =(bestwf(rightIdx) - bestwf(allLocsSorted_idx(2))) / ...
                         (timeVector(rightIdx) - timeVector(allLocsSorted_idx(2)) + eps);

                end

                halfWidth=halfWidth*1000;%ms
                halfWidth=round(halfWidth,2);
                Slope1=Slope1/1000;% uV/s --->uV/ms
                Slope2=Slope2/1000;% uV/s --->uV/ms
                if (Slope1 ~= 0) && (Slope2 ~= 0)%%very straight waveform on a side
                   Slope=abs(min(Slope1,Slope2)); 
                   Slope=round(Slope);
                else
                    Slope=abs(max(Slope1,Slope2));
                    Slope=round(Slope,2);
                end    
                
           
           
           
           end
        catch
            noiseReasons{ix} = 'Noise';
            Reasons{ix} = 'Waveform doesnt look physiological';
            good(ix) = false;
            continue;                   
        end    
        
        chRangeAmp{ix}= peaktopeakChRange;
        amplitudes{ix}=peaktopeak;
        halfwidths{ix}=halfWidth; 
        slopes{ix}=Slope;
        spikeType{ix}=spikingType;
        
    end

    
    
end