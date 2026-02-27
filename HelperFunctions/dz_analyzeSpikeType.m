function [spikingType, quality,passReason,halfWidth,Slope1,Slope2,amplitude] = dz_analyzeSpikeType(waveform, timeVector,allLocsSorted,actualClusterID,channelUsed, minHW,minAmp,minSlope)

    quality=false;
  
    % Convert allLocsSorted (time values) into indices for waveform lookup
    allLocsSorted_idx = arrayfun(@(t) find(timeVector >= t, 1), allLocsSorted);   
    % Ensure indices stay within valid range
    allLocsSorted_idx = min(max(allLocsSorted_idx, 1), length(waveform));
    centralPoint = waveform(allLocsSorted_idx(2));
    centerIndex = allLocsSorted_idx(2); %2nd point is central for analysis

    
    if centralPoint > 0
         spikingType = 'Positive Spiking';
         %%lets determine height of  the peak at the center: amplitude
         amplitude = centralPoint;
         amp=waveform(allLocsSorted_idx);
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


         [~, leftIdx] = min(abs(waveform(allLocsSorted_idx(1):allLocsSorted_idx(2)) - halfAmplitude));
         leftIdx = leftIdx + allLocsSorted_idx(1) - 1;  % Adjust index to global waveform
        
         [~, rightIdx] = min(abs(waveform(allLocsSorted_idx(2):allLocsSorted_idx(3)) - halfAmplitude));
         rightIdx = rightIdx + allLocsSorted_idx(2) - 1;  % Adjust index to global waveform

         halfWidth = timeVector(rightIdx) - timeVector(leftIdx);

         Slope1 =(waveform(allLocsSorted_idx(2)) - waveform(leftIdx)) / ...
            (timeVector(allLocsSorted_idx(2)) - timeVector(leftIdx) + eps);

         Slope2 =(waveform(rightIdx) - waveform(allLocsSorted_idx(2))) / ...
             (timeVector(rightIdx) - timeVector(allLocsSorted_idx(2)) + eps);

         

    else
        spikingType = 'Regular Spiking';

        %%lets  determine bigger peak to use for amplitude
        amp=waveform(allLocsSorted_idx);
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
        
        [~, leftIdx] = min(abs(waveform(allLocsSorted_idx(1):allLocsSorted_idx(2)) - halfAmplitude));
        leftIdx = leftIdx + allLocsSorted_idx(1) - 1;  
        
        [~, rightIdx] = min(abs(waveform(allLocsSorted_idx(2):allLocsSorted_idx(3)) - halfAmplitude));
        rightIdx = rightIdx + allLocsSorted_idx(2) - 1;  

        halfWidth = timeVector(rightIdx) - timeVector(leftIdx);

        Slope1 =(waveform(allLocsSorted_idx(2)) - waveform(leftIdx)) / ...
            (timeVector(allLocsSorted_idx(2)) - timeVector(leftIdx) + eps);
    
        Slope2 =(waveform(rightIdx) - waveform(allLocsSorted_idx(2))) / ...
             (timeVector(rightIdx) - timeVector(allLocsSorted_idx(2)) + eps);

    end

    halfWidth=halfWidth*1000;%ms
    Slope1=Slope1/1000;% uV/s --->uV/ms
    Slope2=Slope2/1000;% uV/s --->uV/ms
    if halfWidth<=minHW   
        if amplitude>minAmp
            if (abs(Slope1)>=minSlope || abs(Slope2)>=minSlope) 
                quality =1;  
                passReason="Waveform passes all checks";
            else
                passReason="Slope too low (<100 uV/ms),potential slow waveform"; %write in ms
                quality =0;
            end    
        else
            passReason="Amplitude is too low (<50uV)";
            quality =0;
        end    
    else   
        passReason="Half Width is too large (>0.45ms)";
        quality =0;
    end    
    halfWidth=round(halfWidth,2);
    % disp(halfWidth);
    % disp(Slope1);
    % disp(Slope2);

    %% Plotting    
    xHalfWidth = [timeVector(leftIdx), timeVector(rightIdx)];
    yHalfWidth = [halfAmplitude, halfAmplitude];
    % 
%     figure;
%     plot(timeVector, waveform, 'b', 'LineWidth', 1.5); hold on;
%     plot(timeVector(allLocsSorted_idx), waveform(allLocsSorted_idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Peaks/Troughs
%     plot(timeVector(centerIndex), waveform(centerIndex), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % Center point
%     if ~isinf(halfWidth)
%         plot([timeVector(leftIdx), timeVector(rightIdx)], ...
%      [halfAmplitude, halfAmplitude], 'k--', 'LineWidth', 1.5); % Half width line
%     end
%     
%     title(sprintf('Spike Type: %s | Quality: %d | Cluster ID: %d | Channel: %d', spikingType, quality, actualClusterID,channelUsed));
%     
%     xlabel('Time (s)');
%     ylabel('Amplitude');
%     legend('Waveform', 'Peaks/Troughs', 'Center Point', 'Half Width');
%     hold off;
    % % % 
    % % 

end