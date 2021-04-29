
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the THETA
% modulators found by Shaoyu's and both the sender and the receiver
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators
%
% In contrast to other codes that employes the coherence-gram to estimate
% the coherence vs frequency, this code employes directly coherency.m
%
% INPUT: sess_data_lfp.mat
%        structure containing all modulator infos + RS LFP split
%
% OUTPUT: txt files with the values of the coherence MR, MS, SR and
% corresponding figures
%
%    @ Gino Del Ferraro, November 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

name_struct_input = '/sess_data_lfp_coherence_fk_200_W_5_movie.mat';

filename = '_movie.mat'; % -- filename for sess_data_info.mat
recording = 'movie';
save_dir = 'movie_corrected';

freq_band = 'theta_band';
monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
fk = 200; W = 5;


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

excluded_sess = [8,22,30,31];
excluded_idx = [2,5,8,9];
sess_list = 1:size(sess_info{1},1);
sess_list(excluded_idx) = [];

cnt_sr = 1; % counter sender-receiver coherencies
cnt_el = 1; % counter for how many modulators excluding the receivers modulators

for i = sess_list %1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    
    
    load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    
    % -- store coherence values sender-receiver and spectrums
    stim(cnt_sr).c_sr = sess_data_lfp.c_sr; % assign S-R coherence value
    stim(cnt_sr).s_s = sess_data_lfp.s_s; % assign sender spectrum
    stim(cnt_sr).s_r = sess_data_lfp.s_r; % receiver spectrum
    cnt_sr = cnt_sr + 1;    % sender/receiver counter
    
    
    mod_Ch = sess_data_lfp.mod_idx; % -- modulators (not controls!) index
    
    
    cnt_m = 1;
    for Ch = mod_Ch % for all the modulators in the session
        
        if Sess == 19 && Ch == 60   % Exclude bad channel
            % do nothing 
        elseif Sess == 41 && Ch == 8 % Exclude bad channel
            % do nothing 
        else
            % -- structure assignements
            mod(cnt_el).c_ms = sess_data_lfp.mod(cnt_m).c_ms ; % assign M-S coherence value for this modulator
            mod(cnt_el).c_mr = sess_data_lfp.mod(cnt_m).c_mr;  % M-R coherence
            mod(cnt_el).s_m = sess_data_lfp.mod(cnt_m).S_m; % Modulator spectrum
            
            cnt_el = cnt_el + 1; % total modulators counter
        end
        cnt_m = cnt_m + 1; % counter for modulators within this session
        
    end
    
end


keyboard

dir_Mod_results = strcat(dir_RS_Theta,sprintf('/Modulators_Controls_avg_results/%s',save_dir));
if ~exist(dir_Mod_results, 'dir')
    mkdir(dir_Mod_results)
end


% Save coherence and spectrum data in structure format
save(strcat(dir_Mod_results,sprintf('/coh_spec_m_fk_%d_W_%d%s',fk,W,filename)),'mod');
save(strcat(dir_Mod_results,sprintf('/coh_spec_sr_fk_%d_W_%d%s',fk,W,filename)),'stim');




