%%%By Diksha Zutshi--modified for chenaging tetsing order
function dz_Curate(basename,datfil, clufile,thresholds,pipeline)
        
    %%Algorithm flow:
    %1 Check for firing rate for non-biological clusters-noise
    %2 Check for noise
    %3 Check for cluster to be a MUA cluster
    %4 Assign as Single Unit

    
    disp('Inside assessNoise');
    disp(['Data file: ', datfil]);
    
    %% reading
    [acgEvaluationMode, maxHW, minAmp,maxAmp, minSlope, firingThreshold,acgallthreshold, acgalllabel,correlationthreshold]=thresholds{:};
    
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
    
        
    posfil= [a filesep 'channel_positions.npy'];
    pos = readNPY(posfil);   % Nx2 matrix
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize variables
    uclu = unique(clu);
    nclu=numel(uclu);
    good = true(nclu, 1);
    noiseReasons = cell(nclu, 1);
    Reasons = cell(nclu, 1); %%%this will hold reasons why a cluster isn't Good.
    passReason=[];
    
    

   
    %% LOW FIRING RATE:
    %%%RECORDING SIZE
    % Calculate the recording duration
    fileInfo = dir(datfil); % Get file information
    fileSize = fileInfo.bytes; % Get file size in bytes
    totalSamples1 = fileSize / (2 * nChannels); % Calculate the total number of samples
    recordingDuration = totalSamples1 / fs; % Duration in seconds
    recordingDurationHours = recordingDuration / 3600; % Duration in hours     
   
    
    %% extracting and fiktering wfs from dat
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
    %extracting the waveforms::: in uV
    clippedWaveforms = wfdata.clippedWaveforms;
    wfs=length(clippedWaveforms);    
    channelCorrelations = wfdata.channelCorrelations;
    corrCell=cell(nclu, 1);
    %% Extracting 3 point waveform and waveform features
    [amplitudes,halfwidths,slopes,spikeType,wfs,good,noiseReasons,Reasons,chRangeAmp]=dz_extractWfVariables(wfdata,uclu,fs,good,noiseReasons,Reasons,pos);  
    
    %% Extracting ACG proportions
    acgfile = fullfile(outputsDir, basename + "acg.mat");
    [Proportions,isEmpty]=dz_extractAcgVariables(uclu,acgfile,clu,ts);  
    
     
  
    %% FUNCTION DEFINITIONS:::::::::::::::::::::::>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
    %% firing rate check
      function dz_evaluateLowFiringRates(recordingDuration,firingThreshold)   
        
        for ix = 1:nclu
            % actualclusterid=20;
            % tst = ts(clu == actualclusterid);          
            if ~good(ix)
                continue;
            end           
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
            else
                good(ix)=true;
            end    
        end 
    end
   
    %% correlation
    function dz_analyzeCorrelation(wfs, channelCorrelations,wfdata,pos)
        bestWaveformsChannel=wfdata.bestWaveformsChannel;
        for ix=1:wfs
%             actualindex=99;
%             ix = find(uclu == actualindex);
            corrValues = channelCorrelations{ix};
            thisBestChannel = bestWaveformsChannel{ix};
            if ~good(ix)
                continue;
            end
            if isempty(corrValues)           
                continue;
            end
            
            bestPos = pos(thisBestChannel,:);
            
            d = sqrt(sum((pos - bestPos).^2,2));
            [~, idx] = sort(d);
            chRange = idx(1:11);   % best channel + 10 neighbors      
            localCorrValues = corrValues(chRange);
            localCorrValues=localCorrValues(2:end);
            totalChannels = length(localCorrValues);
            numHighCorr = sum(localCorrValues > correlationthreshold);
            percentHighCorr = (numHighCorr / totalChannels) * 100;
            if (percentHighCorr>80)
                noiseReasons{ix} = 'Noise';
                Reasons{ix} = sprintf('High Correlation (%.2f%%) on all Channels', percentHighCorr);  
                good(ix)=false;  

            else
                noiseReasons{ix} = 'good';
                Reasons{ix} = sprintf('Passed Correlation Check:(%.2f%%) ', percentHighCorr);   
                good(ix)=true;  
            end 
            corrCell{ix}=percentHighCorr;
           
        end
    end    
 
    %% amplitude
    function dz_analyzeAmplitude(wfs,minAmp,maxAmp,amplitudes,chRangeAmp)%max and min
        for ix=1:wfs
%             actualindex=157;
%             ix = find(uclu == actualindex);
            thisAmplitude=amplitudes{ix}; 
            
            if ~good(ix) 
                continue;
            end
            if isempty(thisAmplitude) || ~isscalar(thisAmplitude) || isnan(thisAmplitude)          
                continue;
            end
            thisChRangeAmp=chRangeAmp{ix};  
            Diff=thisAmplitude-thisChRangeAmp;
            Diff(Diff == 0) = [];     
            totalChannels = length(thisChRangeAmp);
            numHighAmpDiff = sum(Diff > maxAmp);
            percentHighAmpDiff = (numHighAmpDiff / totalChannels) * 100;
            
            if (percentHighAmpDiff>70) || thisAmplitude <= minAmp
                noiseReasons{ix} = 'Noise';
                Reasons{ix} = sprintf('Amplitude is too low or all local channel differences exceed %.3f uV', maxAmp);
                good(ix) = false;
            else
                noiseReasons{ix} = 'good';
                Reasons{ix} = sprintf('Passed amplitude check: amplitude = %.3f uV, not all local differences exceed %.3f uV', ...
                    thisAmplitude, maxAmp);
                good(ix) = true;
            end
           
        end   
    end
   %% halfwidth
    function dz_analyzeHalfWidth(wfs, halfwidths,maxHW)
        for ix=1:wfs
            thisHalfWidth=halfwidths{ix};
            if ~good(ix)
                continue;
            end
            if isempty(thisHalfWidth) || ~isscalar(thisHalfWidth) || isnan(thisHalfWidth)          
                continue;
            end
                
            if thisHalfWidth<=maxHW
                noiseReasons{ix} = 'good';
                Reasons{ix} = sprintf('Passed the halfwidth check: %.3f Hz < %.3f Hz (Threshold)', thisHalfWidth, maxHW);
                good(ix)=true;                            
            else
                noiseReasons{ix} = 'Noise';
                Reasons{ix} = sprintf('Halfwidth too wide: %.3f Hz > %.3f Hz (Threshold)', thisHalfWidth, maxHW);
                good(ix)=false;   
            end               
        end    
    end
    %% slope
    function dz_analyzeSlope(wfs,minSlope,slopes)
        
     for ix=1:wfs
        if ~good(ix)
            continue;
        end
%         actualindex=203;
%         ix = find(uclu == actualindex);
        thisSlope=slopes{ix};%uV/ms
        if isempty(thisSlope) || ~isscalar(thisSlope) || isnan(thisSlope)          
            continue;
        end
        if thisSlope>=minSlope             
            noiseReasons{ix} = 'good';
            Reasons{ix} =  sprintf('Passed the slope check: %.3f Hz > %.3f Hz (Threshold)', thisSlope, minSlope);
            good(ix)=true; 
        else
            noiseReasons{ix} = 'Noise';
            Reasons{ix} =  sprintf('Slope too low,potential slow waveform: %.3f Hz < %.3f Hz (Threshold)', thisSlope, minSlope);
            good(ix)=false; 
        end 
     end  
    end

    %% empty acg check
    function dz_acgEmpty(nclu,isEmpty,Proportions)
        for ix = 1:nclu
             if ~good(ix)
                continue;
             else
%                  actualindex=9;
%                  ix = find(uclu == actualindex);
                 thisEmpty=isEmpty(ix);  
                 thisProportion=Proportions{ix};
                if thisEmpty==true
                    good(ix)=false;
                    noiseReasons{ix} = 'Noise';
                    Reasons{ix} = sprintf('Empty ACG: (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                           thisProportion(1) * 100, thisProportion(2) * 100, thisProportion(4) * 100, thisProportion(5) * 100, thisProportion(6) * 100);
                end    
  
             end
        end     
    end
    %% acgall
    function dz_acgAllCheck(nclu,acgallthreshold, acgalllabel,Proportions)   
        for ix = 1:nclu
%             actualindex=53;
%             ix = find(uclu == actualindex);
            if ~good(ix)
                continue;
            else
                thisProportion=Proportions{ix};
                if all(thisProportion>acgallthreshold)   
                    good(ix)=false;
                    noiseReasons{ix} = acgalllabel;
                    Reasons{ix} = sprintf('Centerbins proportion greater than 80 percent of the shoulder (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                           thisProportion(1) * 100, thisProportion(2) * 100, thisProportion(4) * 100, thisProportion(5) * 100, thisProportion(6) * 100);
                else
                    good(ix)=true;
                    noiseReasons{ix} = 'good';
                    Reasons{ix} = sprintf('Passed ACGall check: (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                           thisProportion(1) * 100, thisProportion(2) * 100, thisProportion(4) * 100, thisProportion(5) * 100, thisProportion(6) * 100);
                end
            end    
                   
        end    
    end

    

    %% mua check
    function dz_muaCheck(nclu,acgEvaluationMode,Proportions)   
        for ix = 1:nclu
            if ~good(ix)
                continue;
            else
               thisProportion=Proportions{ix};
               thisProportion=round((thisProportion*100))/100;
               thisProportion=round((thisProportion*100),4)/100;
               flag=false;       
               switch acgEvaluationMode
                   case 'lenient'    
                        %±2 ms bins must be < 30
                        if thisProportion(1) > 0.3 || thisProportion(6) > 0.3
                            flag = true;
                            disp('MUA: violation at ±2 ms bins (>30%)');
                            passReason = sprintf('High activity in 2ms bins compared to the shoulders — %.2f%%, %.2f%%', thisProportion(1) * 100, thisProportion(6) * 100);


                        %Allow 0th and ±1ms bins to be < 30% if outer bins are as low as 30
                        elseif (thisProportion(1) < 0.3 && thisProportion(6) < 0.3) && all(thisProportion(2:5) < 0.2)

                            flag = false;
                            disp('Single Unit');
                            passReason = sprintf('Passed all the checks — (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                            thisProportion(1) * 100, thisProportion(2) * 100, thisProportion(4) * 100, thisProportion(5) * 100, thisProportion(6) * 100);   

                        else 
                            flag = true;
                            %disp('MUA');
                            passReason = sprintf('High activity in 0th & 1ms bins compared to the shoulders — (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                            thisProportion(1) * 100, thisProportion(2) * 100, thisProportion(4) * 100, thisProportion(5) * 100, thisProportion(6) * 100);

                        end


                   case 'strict'    
                        if any((thisProportion*100)>5)
                           flag = true;
                           disp('MUA: violation at ±2 ms bins');
                           passReason = sprintf('MUA: violation at ±2 ms bins — (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                           thisProportion(1) * 100, thisProportion(2) * 100, thisProportion(4) * 100, thisProportion(5) * 100, thisProportion(6) * 100);           
                        else        
                           flag=false;
                           disp('Single Unit');
                           passReason = sprintf('Passed all the checks — (%.2f%%, %.2f%%, %.2f%%, %.2f%%, %.2f%%)', ...
                            thisProportion(1) * 100, thisProportion(2) * 100, thisProportion(4) * 100, thisProportion(5) * 100, thisProportion(6) * 100);     
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
      end    
                     
    end
  %% checking for pipeline order
    survivalCounts = zeros(1, length(pipeline) + 1);
    survivalLabels = cell(1, length(pipeline) + 1);
    goodCounts = zeros(1, length(pipeline) + 1);

    % before any filtering
    survivalCounts(1) = nclu;
    survivalLabels{1} = 'start';
    goodCounts(1) = nclu;   % before any final labeling, none explicitly marked good yet
    radarLabels = pipeline;
    for k = 1:length(pipeline)
        step = pipeline{k};

        switch step
            case 'lowFiring'
                dz_evaluateLowFiringRates(recordingDuration, firingThreshold);
                radarLabels{k} = 'Low Firing';

            case 'correlation'               
                 dz_analyzeCorrelation(wfs, channelCorrelations,wfdata,pos);
                 radarLabels{k} = 'Correlation';

            case 'amplitude'  
                dz_analyzeAmplitude(wfs, minAmp, maxAmp, amplitudes,chRangeAmp);
                radarLabels{k} = 'Amplitude';

            case 'halfWidth'
                dz_analyzeHalfWidth(wfs, halfwidths, maxHW);
                radarLabels{k} = 'Half Width';

            case 'slope'
                dz_analyzeSlope(wfs, minSlope, slopes);
                radarLabels{k} = 'Slope';
            
            case 'acgEmpty'
                dz_acgEmpty(nclu,isEmpty,Proportions);                
                radarLabels{k} = 'ACG Empty';
                
            case 'acgAll'
                dz_acgAllCheck(nclu, acgallthreshold, acgalllabel, Proportions);
                 radarLabels{k} = 'ACG All';

            case 'mua'
                dz_muaCheck(nclu, acgEvaluationMode, Proportions);
                radarLabels{k} = 'MUA';

            otherwise
                error('Unknown pipeline step: %s', step);
        end
        labels=strings(size(noiseReasons));
        for i=1:numel(noiseReasons)
            x=noiseReasons{i};
            if isempty(x)
                labels(i)="";         
            else
                labels(i)=string(x);
            end    
        end    
        
        labels=strtrim(labels);
        isNeuronal = labels == "" | labels == "good" | labels == "MUA";
        isGoodOnly = labels == "" |labels == "good";
        survivalCounts(k+1) = sum(isNeuronal); % k+1 for next, k is before pipeline starts
        goodCounts(k+1) = sum(isGoodOnly);
        survivalLabels{k+1} = step;
    end
    
    %% PLOTS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%     %% survival curve        
    survivalPercent = 100 * goodCounts / nclu;
   
    fig = figure;
    plot(0:length(pipeline), survivalPercent, '-o', 'LineWidth', 1.5); 
    xticks(0:length(pipeline));
    xticklabels(survivalLabels);
    xtickangle(45);
    xlabel('Pipeline step');
    ylabel('Clusters remaining (%)');
    title('SpikeCleaner survival across pipeline');
    grid on;
    savePath = fullfile(a, 'SpikeCleaner_survival_curve.png');
    saveas(fig, savePath);
    
    % saving as csv
    survivalTable = table( ...
    (0:length(pipeline))', ...
    goodCounts(:), ...
    survivalPercent(:), ...
    'VariableNames', {'StepIndex','ClustersRemaining','PercentRemaining'});
% 
%     writetable(survivalTable, fullfile(a, 'SpikeCleaner_survival_curve.csv'));
    print(fig, fullfile(a, 'SpikeCleaner_survival_curve'), '-dsvg');
    
    %% radar plot for unit survival
    removedPerStep = zeros(1, length(pipeline));

    for k = 1:length(pipeline)
        removedPerStep(k) = survivalCounts(k) - survivalCounts(k+1);
    end

    removedPercent = 100 * removedPerStep / nclu;
    theta = linspace(0, 2*pi, length(removedPercent)+1);
    rho = [removedPercent removedPercent(1)];
    
    fig = figure;
    polarplot(theta, rho, '-o', 'LineWidth', 1.5);
    hold on;
    % red marker on first step
    text(theta(1), rho(1)+10, 'counterclockwise', ...
    'Color', 'r', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
    thetaticks(rad2deg(theta(1:end-1)));
    thetaticklabels(radarLabels);
    title('Percentage of Noisy Clusters Removed by each SpikeCleaner step');
%     
%     savePath = fullfile(a, 'SpikeCleaner_clusterremovalcurve.png');
%     saveas(fig, savePath);
    
    print(fig, fullfile(a, 'SpikeCleaner_clusterremovalcurve'), '-dsvg');
   
    
       
    %% radar plot for unit survival
    removedPerStep = zeros(1, length(pipeline));

    for k = 1:length(pipeline)
        removedPerStep(k) = goodCounts(k) - goodCounts(k+1);
    end

    removedPercent = 100 * removedPerStep / nclu;
    theta = linspace(0, 2*pi, length(removedPercent)+1);
    rho = [removedPercent removedPercent(1)];
    
    fig1 = figure;
    polarplot(theta, rho, '-o', 'LineWidth', 1.5);
    hold on;
    % red marker on first step
    text(theta(1), rho(1)+10, 'counterclockwise', ...
    'Color', 'r', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
    thetaticks(rad2deg(theta(1:end-1)));
    thetaticklabels(radarLabels);
    title('Percentage of Non-Single Unit Clusters Removed by each SpikeCleaner step');
    
%     savePath = fullfile(a, 'SpikeCleaner_clusterremovalcurve.png');
%     saveas(fig1, savePath);
    
    print(fig1, fullfile(a, 'SpikeCleaner_clusterremovalcurve_SUs'), '-dsvg');
   
    
    %% statistics saved
     goodPercent = 100 * goodCounts / nclu;
     neuronalPercent = 100 * survivalCounts(:) / nclu;
     survivalLabels  = survivalLabels(:);
     survivalCounts  = survivalCounts(:);
     goodCounts      = goodCounts(:);
     neuronalPercent = neuronalPercent(:);
     goodPercent     = goodPercent(:);
  
     StepwiseStats = table( ...
        survivalLabels(:), ...
        survivalCounts(:), ...
        neuronalPercent, ...
        goodCounts(:), ...
        goodPercent, ...
        'VariableNames', {'Step','Good_MUACount','Good_MUAPercent','GoodCount','GoodPercent'});
     %% save stepwise stats table as an image
    figTable = figure('Position',[100 100 1600 650]);
    % main title
    annotation(figTable,'textbox', ...
        [0 0.82 1 0.08], ...
        'String','SpikeCleaner Stepwise Survival Statistics', ...
        'EdgeColor','none', ...
        'HorizontalAlignment','center', ...
        'FontSize',14, ...
        'FontWeight','bold');

    % explanatory subtitle
    annotation(figTable,'textbox', ...
        [0 0.8 1 0.07], ...
        'String','All clusters begin as candidate single-units (SU). Pipeline progressively removes non-neuronal clusters. Neuronal = Good + MUA, SU = Good', ...
        'EdgeColor','none', ...
        'HorizontalAlignment','center', ...
        'FontSize',10);
       t = uitable(figTable, ...
    'Data', table2cell(StepwiseStats), ...
    'ColumnName', StepwiseStats.Properties.VariableNames, ...
    'RowName', [], ...
    'Units','Normalized', ...
    'Position',[0.03 0.08 0.94 0.72], ...
    'FontSize',12);
    t.ColumnWidth = {160 160 160 160 160};
% 
%     saveas(figTable, fullfile(a, 'SpikeCleaner_stepwise_stats_table.png'));
%     writetable(StepwiseStats, ...
%     fullfile(a, 'SpikeCleaner_stepwise_stats_table.csv'));
%     saveas(figTable, fullfile(a, 'SpikeCleaner_stepwise_stats_table.fig'));
    
    print(figTable, fullfile(a, 'SpikeCleaner_stepwise_stats_table'), '-dsvg');
    %%
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
        halfWidth=halfwidths{ix};
        if isempty(halfWidth)
            halfWidth = '';  % Default to " " if no specific reason provided
        end
        fprintf(fid1, '%d\t%d\n', cluster_index,halfWidth);

    end
    fclose(fid1);
    disp('Finished processing halfwidth');
        
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%WRITING SLOPES 
    slopess=fullfile(a,'slopes.tsv' );
    % 
    % % Open the file for writing halfwidth
    fid1 = fopen(slopess, 'wt');
    if fid1 == -1
        error('Cannot open file: %s', slopess);
    end 
    fprintf(fid1, 'cluster_id\t slopes \n');

    for ix = 1:length(uclu)  % Loop over the actual cluster IDs (since uclu contains the cluster IDs)
        cluster_index=uclu(ix);
        Slope=slopes{ix};
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
