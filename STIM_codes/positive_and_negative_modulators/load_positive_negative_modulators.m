
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


for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    
    load(strcat(dir_Sess,'/sess_data_stim.mat'));
    
    cnt_m = 1;
   
    for m = sess_data_stim.mod_idx
        
        theta_pow = sess_data_stim.mod(cnt_m).theta_pow;
        [sort_theta, sort_idx] = sort(theta_pow);
        hit = sess_data_stim.hits;
        miss = sess_data_stim.misses;
        L = length(theta_pow); 
        
        acc_list = [];
        for l = 1:L % varying the threshold for the ROC
            pot_miss = sort_idx(1:l); % everything at the left of the threshold
            pot_hit = sort_idx(l+1:end); % everything at the right of the threshold
            
            
            Tmiss = intersect(pot_miss,sess_data_stim.misses);
            Thit = intersect(pot_hit,sess_data_stim.hits);
            acc = (length(Tmiss) + length(Thit))/L;
            acc_list = [l, acc; acc_list];
        end
        [acc_max, idx_max] = max(acc_list(:,2))
        keyboard 
        cnt_m = cnt_m + 1;

    end
    
    
end
