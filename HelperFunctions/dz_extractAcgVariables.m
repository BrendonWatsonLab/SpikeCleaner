%% extract acg variable
function [Proportions,isEmpty]=dz_extractAcgVariables(uclu,acgfile,clu,ts)     
    nclu=numel(uclu);
    acgs = cell(nclu, 1);%to save all the acgs
    if exist(acgfile,'file')
        load(acgfile,'acgs');
    end    
    
    % initializing variables
    Proportions=cell(nclu,1);% of center bins vs shoulders
    isEmpty = false(nclu,1);
    for ix = 1:nclu
   
       actualindex=uclu(ix);
%        actualindex=232;
%        ix = find(uclu == actualindex);
       tst = ts(clu == actualindex);
       binSize= 0.001;%1MS
       
       % making acg
       if exist(acgfile,'file')
           acgData=acgs{ix};
           if isempty(acgData)
                numBins=25;
                [currentacg,binCenters] = dz_autoCorr(tst, tst, binSize, numBins);
                acgs{ix} = struct('currentacg', currentacg, 'binCenters', binCenters);     
           else
                currentacg = acgData.currentacg;
                binCenters = acgData.binCenters;
           end    

       else
            numBins=25;
            [currentacg,binCenters] = dz_autoCorr(tst, tst, binSize, numBins);
            acgs{ix} = struct('currentacg', currentacg, 'binCenters', binCenters);
       end
       
       %% check if empty 15 bins on either side
       zeroLagIndex = ceil(length(currentacg) / 2); 
        fracZero = sum(currentacg == 0) / length(currentacg);
        if fracZero > 0.5
            isEmpty(ix) = true;
            proportion= nan(1,6);
            Proportions{ix}=proportion; % saving it in the cell array            
            continue;
        end
        
        %% if a  full acg
       % finding proportion values
       pred = mean(currentacg(1:10));
       zeroLagRange = (zeroLagIndex - 2):(zeroLagIndex + 2);
       cchWithoutZeroLag = currentacg;
       cchWithoutZeroLag(zeroLagRange) = []; % Remove zero lag values
       if(pred==0)
           pred=mean(cchWithoutZeroLag);
       end 
       cchzero = currentacg(zeroLagIndex-3:zeroLagIndex+2); % Zero lag INDEXES
          
        %derivative based peaks in the shoulders
        % Find first peak after central refractory zone
       searchStart = zeroLagIndex + 3;  % skip 3 bins on right side of center
       seg = currentacg(searchStart:end);

        %%finding the highest peak
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
                [sortedPks, idx] = sort(pks, 'descend');
                % Number of peaks to average
                n = min(3, length(sortedPks));
                % Take top n peaks
                topPeaks = sortedPks(1:n);
                % Average them
                peakVal = mean(topPeaks);                         
            end
        end            
        proportion=cchzero/peakVal;  %proportion of  center  peaks  in  comparison to the shoulder of the  acg
   
        Proportions{ix}=proportion; % saving it in the cell array
        % progress update every 10 clusters
        if mod(ix,10) == 0 || ix == nclu
            fprintf('Processed %d / %d clusters (%.1f%%)\n', ix, nclu, 100*ix/nclu);
        end
              
    end
    %% saving the acg file if it didn't exist
    if ~exist(acgfile,'file')
            save(acgfile,'acgs', '-v7.3');
    end

    
end