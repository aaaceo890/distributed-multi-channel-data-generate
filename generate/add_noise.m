function add_noise(speech_type, speech, max_pnum, org_type, target_type, mic, option, varargin)
    global noise_list diffuse_noise_list point_noise_list
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Junqi Chen
% Copyright 2021 Shanzheng Guan
% Copyright 2021 Ziye Yang
% Center for Intelligent Acoustics and Immersive Communications and School of Marine Science and Technology
% Northwestern Polytechnical University
% If you have any questions, please contact:
% jqchen@mail.nwpu.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    if length(varargin) == 2
        out_dir = cell2mat(varargin(1));
        snr = cell2mat(varargin(2));
        snr_ctr = true;
    elseif length(varargin) == 1
        out_dir = cell2mat(varargin(1));
        snr_ctr = false;
    end
    
   % multi channel 
    if strcmp(mic{1},'single') == 0
        mic_type = mic{1};
        channel_nums = mic{2};
    else
        mic_type = 'single';
        channel_nums = 1; 
    end
   
    [source_speech,fs] = audioread(speech);
    
    % Add samples to account for the reverberation tail
    add_samples = 2500;
    source_speech=[source_speech; zeros(1,add_samples)'];
    
    l = length(source_speech);
    source_speech = source_speech';
    
    if strcmp(speech_type,'train') || strcmp(speech_type,'dev')
        % train, dev mode
        noise_select = unidrnd(length(noise_list));
        startpoint = unidrnd(length(noise_list{noise_select}) - l + 1);
        endpoint = startpoint + l - 1;
        diffuse_noise = noise_list{noise_select}(startpoint:endpoint);

        p_num = unidrnd(max_pnum);
        point_noise = zeros(p_num,l);
        for j = 1:p_num
            noise_select = unidrnd(length(noise_list));
            startpoint = unidrnd(length(noise_list{noise_select}) - l + 1);
            endpoint = startpoint + l - 1;
            point_noise(j,:) = noise_list{noise_select}(startpoint:endpoint);
        end
    elseif strcmp(speech_type,'test')
        % test mode
        noise_select = 1;
      %  noise_select = unidrnd(length(diffuse_noise_list));
        startpoint = unidrnd(length(diffuse_noise_list{noise_select}) - l + 1);
        endpoint = startpoint + l - 1;
        diffuse_noise = diffuse_noise_list{noise_select}(startpoint:endpoint);

        p_num = unidrnd(max_pnum);
        point_noise = zeros(p_num,l);
        for j = 1:p_num
            noise_select = 1;
          %  noise_select = unidrnd(length(point_noise_list));
            startpoint = unidrnd(length(point_noise_list{noise_select}) - l + 1);
            endpoint = startpoint + l - 1;
            point_noise(j,:) = point_noise_list{noise_select}(startpoint:endpoint);
        end
    end
    
    while 1
        try
            if snr_ctr
                snr_out = num2str(snr);
                [noise_speech,direct_speech,bkg_noise,snr,info] = noise_speech_generate(source_speech,diffuse_noise,mic,point_noise,snr);
            else
                snr_out = '';
                [noise_speech,direct_speech,bkg_noise,snr,info] = noise_speech_generate(source_speech,diffuse_noise,mic,point_noise);
            end
            break
        catch
            % fprintf(dt);
        end
    end
   % saving the output wavfile
   
   for op = 1:length(option.value)
        if option.value == '0'
            continue
        end
        
        switch option.key{op}
            case 'noise_speech'
                name_wav = strsplit(getfield(out_dir,'noise_speech'),org_type);
                for ch = 1:channel_nums
                    full_name = [name_wav{1},['-',mic_type,num2str(ch)],target_type];
                    audiowrite(full_name,noise_speech(ch,:),fs)
                end
            case 'direct_speech'
                name_wav = strsplit(getfield(out_dir,'direct_speech'),org_type);
                for ch = 1:channel_nums
                    full_name = [name_wav{1},['-','direct',mic_type,num2str(ch)],target_type];
                    audiowrite(full_name,direct_speech(ch,:),fs)
                end
            case 'bkg_noise'
                name_wav = strsplit(getfield(out_dir,'bkg_noise'),org_type);
                for ch = 1:channel_nums
                    full_name = [name_wav{1},['-','bkg',mic_type,num2str(ch)],target_type];
                    audiowrite(full_name,bkg_noise(ch,:),fs)
                end
            case 'info'
                name_wav = strsplit(getfield(out_dir,'noise_speech'),org_type);
                info.name = name_wav{1};
                fname = fullfile(getfield(out_dir,'info'),'info.txt');
                write_info(fname,info);
        end
   end