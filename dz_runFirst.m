function dz_runFirst(curdir)
%%By Diksha Zutshi
%run first: to chnage the shame of the data folder to make it compatible.
rootdir=fileparts(mfilename('fullpath'));
addpath(genpath(rootdir));

if nargin < 1 || isempty(curdir)
    curdir=pwd;
end    


requiredFiles={
'amplitudes.npy'
'channel_map.npy'
'channel_positions.npy'
'params.py'
'spike_clusters.npy'
'pc_features.npy'
'spike_times.npy'
'templates.npy'
'spike_templates.npy'
'whitening_mat.npy'
'whitening_mat_inv.npy'
'similar_templates.npy'
};
%%create SpikeCleaner folder
outputDir=fullfile(curdir,'SpikeCleaner');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

%%check for requiredFiles in curdir and copy to SpikeCleaner folder
missing = {};

for i = 1:numel(requiredFiles)
    if ~isfile(fullfile(curdir, requiredFiles{i}))
        missing{end+1} = requiredFiles{i}; 
    end
end   
if isempty(missing)
    for i=1:numel(requiredFiles)
        copyfile(fullfile(curdir,requiredFiles{i}), fullfile(outputDir, requiredFiles{i}));
    end    
    disp('All files have been copied to SpikeCleaner');
end    
%%if not present look for them in all present subfolders
for i=1:numel(requiredFiles)
    fname=requiredFiles{i};
    hits=dir(fullfile(curdir,'**',fname)); %%recursively looking inside subfolders
    if isempty(hits)
        warning('Couldnt find these files under current directory:%s',fname);
    end    
    %%if multiple matches are found then keep the latest
    [~,idx]=max([hits.datenum]);
    src=fullfile(hits(idx).folder,hits(idx).name);
    
    %%copy these files to SpikeCleaner folder
    copyfile(src,fullfile(outputDir,fname));
    fprintf('Copied %s from %s /n',fname,hits(idx).folder);
end    

%%copy the dat file as well SpikeCleaner:: just get the path
datfile = dir(fullfile(curdir, '*.dat')); 
[~,basename]=fileparts(curdir);

matchIdx=find(strcmp({datfile.name}, [basename '.dat']),1);
fullname=[basename '.dat'];


if isempty(matchIdx)
    matchIdx = find(strcmp({datfile.name}, 'temp_wh.dat'), 1);
end

if isempty(matchIdx)
    warning('No matching .dat file found.');
    %check subfolders
    hits=dir(fullfile(curdir,'**',fullname)); %%recursively looking inside subfolders
    if isempty(hits)
        warning('Couldnt find these files under current directory:%s',fullname);
        
    else
         %%copy file to SpikeCleaner folder
        src=fullfile(hits.folder,hits.name);
        %copyfile(src,fullfile(outputDir,fullname));
        %fprintf('Copied %s from %s /n',fullname,hits.folder);
        
        
    end 
    
else
%     copyfile( ...
%         fullfile(curdir, datfile(matchIdx).name), ...
%         fullfile(outputDir, datfile(matchIdx).name) );
    src=fullfile(curdir, datfile(matchIdx).name);
end

disp("All files have been copied");



%%make a parameters.mat file similar to a rez.mat
%needs sampling rate,num of channels,name of the animal from basename
%
parameters=struct();
parameterFile=fullfile(outputDir,'parameters.mat');

params=fullfile(outputDir,'params.py');
contents=fileread(params);

%Renaming the dat file 
newPath=src;
lines = regexp(contents, '\r?\n', 'split')';

for i=1:length(lines)
    if contains(lines{i},'dat_path')
        lines{i}=sprintf("dat_path= '%s' ",newPath);
        lines{i} = sprintf("dat_path = '%s'", newPath);
    end    
end    
    

%parsing animal name
parameters.animal=basename;

%parsing channels
chanMatch=regexp(contents,'n_channels_dat\s*=\s*(\d+)','tokens','once');
parameters.nChannels=str2double(chanMatch{1});

%parsing sampling rate
fsmatch=regexp(contents,'sample_rate\s*=\s*([0-9\.]+)','tokens','once');
parameters.samplingRate=str2double(fsmatch{1});

%%writing back to params.py
fid=fopen(params,'w');
for i = 1:length(lines)
    fprintf(fid, '%s\n', lines{i});
end
fclose(fid);

save(parameterFile,'parameters');%saving the structure

disp('parameters.mat created in SpikeCleaner directory in the animal folder');
end



