
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code considers the stimulation experiments. For each trial, it
% checks whether the receiver's response is associated with high-(low)-theta
% power of the modulator's electrode. Modulators with high-theta power for  
% receiver response are considered "positive", modulators with low theta 
% power for receiver response are considered "negative". 
%
%    @ Gino Del Ferraro, March 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - PATHS and NAMES --- %
%%%%%%%%%%%%%%%%%%%
addpath('/mnt/pesaranlab/Matlab/monkeys')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
dir_high_low_theta = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_Stim_Theta = strcat(dir_main,sprintf('/%s/Stim_data/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

name_structure = '/modulators_decod_accuracy.mat';


% STIM paths and names %%%%%%%%%%%%%%%%%
subjects = {'maverick','archie'};

%for
iSubject = 1% : length(subjects) % Loop on the animals
%     clearvars -except subjects iSubject
if strcmp(subjects{iSubject},'archie')
    archie_vSUBNETS220_rig3
else
    maverick_vSUBNETS220_rig3
end
PreStimSess = PreStimResponseAll_Database_NetworkEdge;

% positive and negative modulators 
positive = 0;
negative = 0;
tot_m  = 0; % total number of modulators across sessions

for s = 1:size(sess_info{1},1)  % For each session with at least one modulator

    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_structure))
    load(strcat(dir_Sess,'/Data_with_theta_band.mat')); % load stim data for info about hits/misses
    load(strcat(dir_Modulators,'/session_data_lfp.mat')); % load RS data for info about the modulators 
    
    
    % hit and miss trials -- they are the same for each modulator within the same session
    hit = Data.spec.lfp.DetectedIndx{1};
    miss = Data.spec.lfp.notDetectedIndx{1};
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD LFP STIM ...                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    StimTrials = Data.StimTrials(Data.goodTrials_index);
    sys = StimTrials(1).MT;
    bn_Pre = [-1005 -5]; % ms
    
    % load all the channel LFPs for that given session 
    [Lfp_Pre] = trialStimPulseLfp(StimTrials, sys, [], [], 'PulseStarts', bn_Pre); % returns monopolar recording

    % %%%%%%%%%%%%%%%%% --->>> LFP data loaded at this point
    
    
    % -- load list electrodes, sender, receiver
    electrode = sess_data_lfp.RecordPair; % ---- all electrode pairs
    receiver = sess_data_lfp.receiver_pair;  % ---- receiver pair
    sender = sess_data_lfp.sender_pair; % ---- sender pair
    
    % %%%%%% Assign LFP %%%%%%
    lfp_E = sq(Lfp_Pre(:,electrode(:,1),:) - Lfp_Pre(:,electrode(:,2),:)); % all electrodes lfp
    lfp_R = sq(Lfp_Pre(:,receiver(1),:) - Lfp_Pre(:,receiver(2),:)); % receiver lfp
    lfp_S = sq(Lfp_Pre(:,sender(1),:) - Lfp_Pre(:,sender(2),:)); % sender lfp
    
    
    sess_data_stim.sess_idx = sess_data_lfp.sess_idx;
    sess_data_stim.day = sess_data_lfp.day;
    sess_data_stim.rec_STIM = sess_data_lfp.rec_STIM;
    sess_data_stim.rec_RS = sess_data_lfp.rec_RS;
    sess_data_stim.RecordPair = sess_data_lfp.RecordPair;
    
    sess_data_stim.MRIlabels = sess_data_lfp.MRIlabels;
    sess_data_stim.RecordPairMRIlabels = sess_data_lfp.RecordPairMRIlabels;
    sess_data_stim.Spec = sess_data_lfp.Spec;
    sess_data_stim.sender_pair = sess_data_lfp.sender_pair;
    sess_data_stim.sender_area = sess_data_lfp.sender_area;
    sess_data_stim.receiver_pair = sess_data_lfp.receiver_pair;
    sess_data_stim.receiver_idx = sess_data_lfp.receiver_idx;
    sess_data_stim.receiver_area = sess_data_lfp.receiver_area;
    sess_data_stim.mod_idx = sess_data_lfp.mod_idx;
    sess_data_stim.mod_areas = sess_data_lfp.mod_areas;
    sess_data_stim.Decod_Accuracy = mod_accuracy.Decod_Accuracy;
    sess_data_stim.auc = mod_accuracy.auc;
    sess_data_stim.se = mod_accuracy.se;
    
    sess_data_stim.lfp_E = lfp_E;
    
    cnt_m = 1;
    for m = sess_data_lfp.mod_idx
        
        display(['-- Modulator ',num2str(m),' --  ',num2str(cnt_m),', out of tot  ',num2str(length(sess_data_lfp.mod_idx)),' '])

        lfp_m = sq(lfp_E(:,m,:)); % modulator lfp
        % Compute the spectrum for each trial. Format: iTrial x times
        W = 3; % frequency smoothing 
        [spec, f,err] = dmtspec(lfp_m(:,501:end),[500/1e3,W],1e3,200); % spectrum 500 ms before onset 
        
        % Find low and high theta from the spectrum 
        theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
        theta_pow_mean = mean(theta_pow); % get the average theta power
        theta_pow = theta_pow - theta_pow_mean; % rescale the theta power by removing the mean value
        
        [sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
        cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution
        
        % low and high theta power indexes 
        low_theta_pow = sort_theta(1:cut);
        low_idx = trial_idx(1:cut);
        
        high_theta = sort_theta(end-cut+1:end);
        high_idx = trial_idx(end-cut+1:end);
        
        HH = intersect(hit,high_idx); % high power - hit trials, true positive
        ML = intersect(miss,low_idx); % low power - miss trials, true negative
        HL = intersect(hit,low_idx); % low power - miss trials, false negative
        MH = intersect(miss,high_idx); % high power - hit trials, false positive
        
        HH_rate = length(HH)/cut;
        ML_rate = length(ML)/cut;
        HL_rate = length(HL)/cut;
        MH_rate = length(MH)/cut;
        
        confusion = [HH_rate, HL_rate; MH_rate, ML_rate]; % confusion matrix 
        
        % store theta power, idx, and confusion matrix for each modulator 
        sess_data_stim.mod(cnt_m).theta_pow = theta_pow;
        sess_data_stim.mod(cnt_m).low_theta_pow = low_theta_pow;
        sess_data_stim.mod(cnt_m).low_idx = low_idx;
        sess_data_stim.mod(cnt_m).high_theta = high_theta;
        sess_data_stim.mod(cnt_m).high_idx = high_idx;
        sess_data_stim.mod(cnt_m).confusion = confusion;
        
        % count number of positive and negative modulators 
        if HH_rate + ML_rate > HL_rate + MH_rate
            positive = positive + 1;
        else 
            negative = negative + 1;
        end
        
        tot_m = tot_m + 1;
        cnt_m = cnt_m + 1;

    end %  for each modulator within session 
      
    save(strcat(dir_Sess,'/sess_data_stim.mat'),'sess_data_stim');
    
end  % for each session



for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Sess,'/sess_data_stim.mat'));
    load(strcat(dir_Modulators,name_structure))

    
    sess_data_stim.Decod_Accuracy = mod_accuracy.Decod_Accuracy;
    sess_data_stim.auc = mod_accuracy.auc;
    sess_data_stim.se = mod_accuracy.se;
    
    
end 


