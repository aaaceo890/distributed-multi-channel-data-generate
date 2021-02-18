function [noise_speech,direct_speech,bkg_noise,snr,info] = noise_speech_generate(source_speech,diffuse_noise,mic,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Junqi Chen
% Copyright 2021 Shanzheng Guan
% Copyright 2021 Ziye Yang
% Center for Intelligent Acoustics and Immersive Communications and School of Marine Science and Technology
% Northwestern Polytechnical University
% If you have any questions, please contact:
% jqchen@mail.nwpu.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-[cs   dn   (pn   snr)]-%
% Import RIR-Generator
addpath('../RIR-Generator-master');

if length(varargin) == 2
    % p_nums X samples
    point_noise = cell2mat(varargin(1));
    snr_control = cell2mat(varargin(2));
    snr_ctrl_flag = 1;
    
elseif length(varargin) == 1
    if length(cell2mat(varargin(1))) == 1
        % Meaning that no pn input but snr_ctrl
        snr_control = cell2mat(varargin(1));
        snr_ctrl_flag = 1;
        point_noise = zeros(size(diffuse_noise));
    else
        snr_ctrl_flag = 0;
        point_noise = cell2mat(varargin(1));
    end
else
    snr_ctrl_flag = 0;
    point_noise = zeros(size(diffuse_noise));
end

% multi channel config
if strcmp(mic{1},'single') == 0
   mic_type = mic{1};
   channel_nums = mic{2};
else
   mic_type = 'single';
   channel_nums = 1; 
end

% Fixed Hyper-Parameters
c = 340;                    % Sound velocity (m/s)

% Semi-Fixed Hyper-Parameters
fs = 16000;                 % Sample frequency (samples/s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = 2;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

%#########################################################################%
% Hyper-Parameters
mode = 2; % 1 means plane, 2 means stereo.
snr_range_for_diffuse_noise = [5,25]; % snr = 10*lg(E_cleanspeech/E_diffusenoise), snr pool for dn.
snr_range_for_point_noise = [-5,25]; % snr = 10*lg(E_cleanspeech/E_pointnoise), snr pool for pn.

% snr_range = [0,30]; % snr = 10*lg(E_directspeech/E_noise), macro-control the snr of noise speech.
% snr_sample_mean = snr_range(1) + 5;snr_sample_std = 10;

% speech_source_nums = 1; % the number of speech source.
% max_point_noise_nums = 4; % the number of point noise source.

room_shape.length = [5 25]; % room shape range (m)
room_shape.width = [5 25];
room_shape.height = [2.7 4];

min_distance_source2wall = 0.2; % min distance from sound source to wall.(m)
min_distance_source2mic = 0.3; % min distance from sound source and point noise source to mic. (m)
T60 = [0 0.8]; % T60 range. (s)
T60_sample_mean = 0.65;T60_sample_std = 0.1;
T60_disturb_mean = 0; T60_disturb_std = 0.005;
%#########################################################################%

%#########################################################################%
% 1 X s
source_speech = normalize(source_speech); % normalize the sig, just in case.
diffuse_noise = normalize(diffuse_noise);

if max(max(point_noise)) > 0
    % if it is all zeros in pn, (for normalize).
    % p_nums X s
    point_noise = normalize(point_noise);
end
% adjust their range based on SNRs.
% --- diffuse noise ---%
snr_dn  = random_uniform_sample(snr_range_for_diffuse_noise);
alpha = sqrt(1/(10^(snr_dn/10))); % "1" for E_cs, because of the normalization of cs.
% c X s
diffuse_noise = ones(channel_nums,1)*diffuse_noise;
diffuse_noise_adjusted = alpha*diffuse_noise;
% --- point noise ---%
[point_noise_nums,idx] = min(size(point_noise)); % get the nums of point noise source.
point_noise_adjusted = zeros(size(point_noise));
for i = 1:point_noise_nums
    snr_pn = random_uniform_sample(snr_range_for_point_noise);
    beta = sqrt(1/(10^(snr_pn/10))); % "1" for E_cs, because of the normalization of cs.
    if idx == 1
        point_noise_adjusted(i,:) = beta*point_noise(i,:);
    else
        point_noise_adjusted(:,i) = beta*point_noise(:,i);
    end
end
%#########################################################################%

%#########################################################################%
% Cal the RIR for cs ds and pn based on the Hyper-Parameters.
room_shape = generate_room_shape(room_shape); % get the room shape.
% --insert T60 low bound.
V = room_shape(1)*room_shape(2)*room_shape(3);
S = 2*(room_shape(1)*room_shape(2)+room_shape(3)*room_shape(2)+room_shape(1)*room_shape(3));
T60_LB =24*log(10)*V/(c*S);
T60(1) = T60_LB + 0.0001;

% get the source position (1 + point_noise_nums).
% source_pos -> (source + point_noise) X 3
source_pos = generate_source_pos(room_shape,point_noise_nums+1,min_distance_source2wall); % [point_noise_nums+1,3] = size(pos).

% get the rec mic position.mode 1 means plane, 2 means stereo.
% rec_pos -> channel_nums X 3
max_iteration_nums = 1e6;
for i = 1:max_iteration_nums
    flag = 0;
    switch mic_type
        case {'single','ad-hoc'}
            rec_pos = generate_rec_pos(room_shape,channel_nums,mode); % With constraint�� min distance from all sources to mic��
        case 'linear10cm'
            d = 0.1;
            rec_pos = generate_liner_pos(room_shape,channel_nums,d,mode);
        otherwise 
            fprintf("not implement for mictype:%s",mic_type)
    end
    
    for ch = 1:channel_nums
        dis_temp = ones(point_noise_nums+1,1)*rec_pos(ch,:)-source_pos;
        if min(sqrt(sum((dis_temp.^2)'))) < min_distance_source2mic
            flag = 1;
            break
        end
    end
    if flag == 0
        [~,near_mic_idx] = min( sum((rec_pos-source_pos(1,:)).^2,2) );
        break
    end
end
% cal the rir for direct speech.
speech_pos = source_pos(1,:);
beta_direct_speech = 0;
h_direct_speech = [];
for ch = 1:channel_nums
    h_direct_speech = [h_direct_speech;rir_generator(c, fs, rec_pos(ch,:), speech_pos, room_shape, beta_direct_speech, n, mtype, order, dim, orientation, hp_filter)];
end
% h_direct_speech -> channel_nums X s 


% cal the rir for reverberation speech.
beta_reverberation_speech = random_nonuniform_sample(T60,T60_sample_mean,T60_sample_std); % get the T60 (s)
beta_reverberation_speech = ones(channel_nums,1) * beta_reverberation_speech;
h_reverberation_speech = [];
for ch = 1:channel_nums
    disturb = normrnd(T60_disturb_mean,T60_disturb_std);
    beta_reverberation_speech(ch) = beta_reverberation_speech(ch) + disturb;
    h_reverberation_speech = [h_reverberation_speech;rir_generator(c, fs, rec_pos(ch,:), speech_pos, room_shape, beta_reverberation_speech(ch), n, mtype, order, dim, orientation, hp_filter)];
end
% h_reverberation_speech -> channel_nums X s 

% cal the rir for point noise.
h_reverberation_noise = []; % point_noise_nums*n_samples
for i = 1:point_noise_nums
    beta_reverberation_noise = random_nonuniform_sample(T60,T60_sample_mean,T60_sample_std); % get the T60 (s)
    h_reverberation_noise = [h_reverberation_noise;rir_generator(c, fs, rec_pos, source_pos(i+1,:), room_shape, beta_reverberation_noise, n, mtype, order, dim, orientation, hp_filter)];
end
% h_reverberation_noise -> p_nums X s 

%#########################################################################%

%#########################################################################%
% generate the direct_speech, reveberation speech and so on.
% one convolution operation per time
% # conv for ds (direct speech).
ds = [];
for ch = 1:channel_nums 
    ds = [ds;conv(source_speech,h_direct_speech(ch,:))];
end
% ds -> c X s
% # conv for rs (reveberation speech).
rs = [];
for ch = 1:channel_nums 
    rs = [rs;conv(source_speech,h_reverberation_speech(ch,:))];
end
% rs -> c X s
% # non-conv for dn (diffuse noise).
try
    diffuse_noise_adjusted = [diffuse_noise_adjusted diffuse_noise_adjusted(:,1:n-1)];
catch
    diffuse_noise_adjusted = [diffuse_noise_adjusted;diffuse_noise_adjusted(:,1:n-1)];
end
dn = diffuse_noise_adjusted;
% dn -> c X s


% # conv for pn (point noise).
[~,idx] = min(size(point_noise_adjusted));
if idx ~= 1
    point_noise_adjusted = point_noise_adjusted';
end
pn = zeros(1,max(size(point_noise_adjusted))+size(h_reverberation_noise,2)-1);
for i = 1:point_noise_nums
    pn = pn + conv(point_noise_adjusted(i,:),h_reverberation_noise(i,:));
end
if idx ~= 1
    pn = pn';
end
pn = ones(channel_nums,1)*pn;
% pn -> c X s
%#########################################################################%

%#########################################################################%
E_s = mean(ds'.^2); % If we accept that reverb does not contribute noise,it's (rs).or else (ds)
E_n = mean((dn + pn)'.^2);
if snr_ctrl_flag
    gamma = sqrt((E_s(near_mic_idx)/(10^(snr_control/10)))/E_n((near_mic_idx)));
else
    gamma = 1;
end
noise_speech = rs + (ones(channel_nums,1) * gamma).*(dn+pn);
snr = 10*log10(E_s./(gamma^2*E_n));
% The sig that needs to be saved.
direct_speech = ds;
bkg_noise = noise_speech - ds;
%#########################################################################%
am_adj_ns = amplitude_limit_dtc(noise_speech);
am_adj_ds = amplitude_limit_dtc(direct_speech);
am_adj_bn = amplitude_limit_dtc(bkg_noise);

am_adj = min([am_adj_ns am_adj_ds am_adj_bn]);

noise_speech = am_adj*noise_speech;
direct_speech = am_adj*direct_speech;
bkg_noise = am_adj*bkg_noise;

% get info
info.source = source_pos(1,:);
info.point = source_pos(2:1+point_noise_nums,:);
info.mic = rec_pos;
info.T60 = beta_reverberation_speech;
info.snr = snr;

end

function sig = normalize(sig_org)
eps = 1e-10;
sig = zeros(size(sig_org));
% Normalize the signal.
[flag,idx] = min(size(sig_org));
if flag == 1
    sig_org = sig_org - mean(sig_org);
    sig = sig_org/max(eps, std(sig_org));
else
    if idx == 1
        for i = 1:flag
            sig_org(i,:) = sig_org(i,:) - mean(sig_org(i,:));
            sig(i,:) = sig_org(i,:)/max(eps, std(sig_org(i,:)));
        end
    else
        for i = 1:flag
            sig_org(:,i) = sig_org(:,i) - mean(sig_org(:,i));
            sig(:,i) = sig_org(:,i)/max(eps, std(sig_org(:,i)));
        end
    end
end
end

function rec_pos = generate_rec_pos(L,channel_nums,mode)
% mode = 1 means plane, mode = 2 means stereo.
% L (1*3) Matrix
candidate_rec_pos = rand(channel_nums,3).*(ones(channel_nums,1)*L);

if mode == 1
    candidate_rec_pos(:,3) = 0;
end
rec_pos = candidate_rec_pos;
end

function rec_pos = generate_liner_pos(L,channel_nums,d,mode)
% mode = 1 means plane, mode = 2 means stereo.
% L (1*3) Matrix
x_region = [0,L(1),L(1),0,0];
y_region = [0,0,L(2),L(2),0];
while 1
    % random point in room
    rand_point = rand(1,3).*L;
    % random angle from -pi to pi
    rand_angle = (rand(1)*2-1)*pi;
    len_arr = (channel_nums-1) * d;
    end_point = rand_point + len_arr*[cos(rand_angle),sin(rand_angle),0];
    if inpolygon(end_point(1),end_point(2),x_region,y_region)
        break
    end
end
rec_pos = zeros(channel_nums,3);
for i = 1:channel_nums
    rec_pos(i,:) = rand_point + d*(i-1)*[cos(rand_angle),sin(rand_angle),0];
end
if mode == 1
    candidate_rec_pos(:,3) = 0;
end
rec_pos = candidate_rec_pos;
end


function L = generate_room_shape(room_shape)
% room_shape
L = [room_shape.length(1) + (room_shape.length(2) -room_shape.length(1)).*rand(1),...
    room_shape.width(1) + (room_shape.width(2) -room_shape.width(1)).*rand(1),...
    room_shape.height(1) + (room_shape.height(2) -room_shape.height(1)).*rand(1)];
end

function source_pos = generate_source_pos(L,source_nums,min_distance_source2wall)
% mode = 1 means plane, mode = 2 means stereo.
% L (1*3) Matrix
L = [L(1)-2*min_distance_source2wall,L(2)-2*min_distance_source2wall,L(3)];
candidate_rec_pos = rand(source_nums,3).*(ones(source_nums,1)*L);
candidate_rec_pos = [candidate_rec_pos(:,1:2)+min_distance_source2wall candidate_rec_pos(:,3)];
source_pos = candidate_rec_pos;

end

function x = random_uniform_sample(a)
% random sampling from [a(1),a(2)]
x = rand(1)*(a(2)-a(1)) + a(1);
end

function x = random_nonuniform_sample(a,p_mean,p_std)
% random sampling from [a(1),a(2)]
if nargin <= 1  % prefer-pra
    if a(2) > 2
        p_mean = 15;
        p_std = 12;
    else
        p_mean = 0.05;
        p_std = 0.1;
    end
end
x = 1e4; % start iteration.
while x<a(1)||x>a(2)
    x = p_mean + p_std*randn;
end
end

function am_adj = amplitude_limit_dtc(x)
% if many sample in a piece > 1 or <-1
[~,idx] = max(size(x));
xx = abs(x);
num_over_bound = sum(xx>1,idx);
gate = 0.5e-3; %1e-3;
if max(num_over_bound) >= gate
   am_adj = 1/max(xx);
else
   am_adj = 1;
end

end
