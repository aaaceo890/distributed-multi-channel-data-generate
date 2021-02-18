clear
clear global
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Junqi Chen
% Copyright 2021 Shanzheng Guan
% Copyright 2021 Ziye Yang
% Center for Intelligent Acoustics and Immersive Communications and School of Marine Science and Technology
% Northwestern Polytechnical University
% If you have any questions, please contact:
% jqchen@mail.nwpu.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% workdir
mpath = mfilename('fullpath');
i = strfind(mpath,'/');
mpath=mpath(1:i(end));
cd(mpath)

% librispeech config
librispeech_dir = '';
% setname={'dev-clean','test-clean','dev-other','test-other',... 
setname = {'train-clean-100'}; 
debug_mode = true;

filetype = '.flac';
target_type = '.wav';

% simulation config
mic = {'ad-hoc',16}; 

option.value = [1,1,1,1]; % generate 1)noise_speech 2)direct_speech 3)bkg_noise 4)info
option.key = {'noise_speech','direct_speech','bkg_noise','info'};
out_dir = struct(option.key{1},'',option.key{2},'',option.key{3},'',option.key{4},'');
% noise config
noise_dir = '';
diffuse_noise_dir = '';
point_noise_dir = '';
max_pnum = 4;

% change target_dir to fit different tasks
target_dir = '';

% debug config
if debug_mode
    librispeech_dir = '../';
    setname = {'test-minidata'};
    noise_dir = '';
    diffuse_noise_dir = '../noise_example'; 
    point_noise_dir = '../noise_example';
    target_dir = '../test_generate_simu_data';
end

% generate all sets
for s = 1:length(setname)
    % mode: train, dev or test
    mode = strsplit(setname{s},'-');
    mode = mode{1};
    %% Creation of the output folder
    if ~exist(target_dir,'dir')
        mkdir(target_dir); 
    end
    
    if ~exist(fullfile(target_dir,setname{s}),'dir')
        mkdir(fullfile(target_dir,setname{s}));
    end

    % Generation of a set folder with the same structure of the original  folder
    fprintf('Folders creation...\n');
    % Create -> orginal_folder, target_folder
    for op = 1:length(option.value)
        if option.value(op) == 0
            fprintf('Not create "%s", skip...\n',option.key{op});
            continue
        end
        if strcmp(option.key{op},'info')
            if exist(fullfile(target_dir,setname{s},option.key{op}),'dir')
                rmdir(fullfile(target_dir,setname{s},option.key{op}),'s')
            end
            mkdir(fullfile(target_dir,setname{s},option.key{op}))
            continue
        end
        if ~exist(fullfile(target_dir,setname{s},option.key{op},'done.txt'),'file')     
        create_folder_str(fullfile(librispeech_dir,setname{s}),fullfile(target_dir,setname{s},option.key{op}));
        fid=fopen(fullfile(target_dir,setname{s},option.key{op},'done.txt'),'a+');
        fclose(fid);
        else
            fprintf('Folders for "%s" exist, skip...\n',option.key{op}); 
        end
    end
           
    % list of all the original set files
    list=find_files(fullfile(librispeech_dir,setname{s}),filetype);


    %% Make noise lists
    global noise_list diffuse_noise_list point_noise_list

    % largescale noise
    if ~strcmp(mode,'test') && isempty(noise_list)
        noise_path = dir(fullfile(noise_dir,'*.wav'));
        noise_list = cell(length(noise_path),1);
        % for i = 1:1
        for i = 1:length(noise_path)
            [noise,fs] = audioread(fullfile(noise_dir,noise_path(i).name));
            noise_list(i) = {noise'};
        end
    elseif strcmp(mode,'test')
        fprintf('mode "%s" not in train or dev, skip...\n', mode)
    else
        disp('noise_list exist, skip...')
    end
    
    % diffuse noise
    if strcmp(mode,'test') && isempty(diffuse_noise_list)
        noise_path = dir(fullfile(diffuse_noise_dir,'*.wav'));
        diffuse_noise_list = cell(length(noise_path),1);
        % for i = 1:1
        for i = 1:length(noise_path)
            [noise,fs] = audioread(fullfile(diffuse_noise_dir,noise_path(i).name));
            diffuse_noise_list(i) = {noise'};
        end
    elseif ~strcmp(mode,'test')
        fprintf('mode "%s" not in test, skip...\n', mode)
    else
        disp('diffuse_noise_list exist, skip...')
    end
    
    % point noise
    if isempty(point_noise_list)
        noise_path = dir(fullfile(point_noise_dir,'*.wav'));
        point_noise_list = cell(length(noise_path),1);
    % for i = 1:1
        for i = 1:length(noise_path)
            [noise,fs] = audioread(fullfile(point_noise_dir,noise_path(i).name));
            point_noise_list(i) = {noise'};
        end
    else
        disp('point_noise_list exist, skip...') 
    end

    disp("noise list accomplished");

    %% Add noise
    fprintf("Adding noise on %s\n", setname{s});
    text = [];
    % for i = 501:650
    for i = 1:length(list)
        text = ['Processing... ',num2str(i),'/',num2str(length(list))]; 
        disp(text)
        
        % TODO
        for op = 1:length(option.value)
            if strcmp(option.key{op},'info')
                out_dir = setfield(out_dir,option.key{op},fullfile(target_dir,setname{s},option.key{op}));
                continue
            end
            out_dir = setfield(out_dir,option.key{op},strrep(list{i},fullfile(librispeech_dir,setname{s}),fullfile(target_dir,setname{s},option.key{op})));
        end
        % name_wav=strrep(list{i},fullfile(librispeech_dir,setname{s}),fullfile(target_dir,setname{s}));
        add_noise(mode,list{i},max_pnum,filetype,target_type,mic,option,out_dir);   
    end
    fprintf('\n');
    fprintf("Successful generate %s\n", setname{s});
    
end

fprintf('\n');
disp("Done");

