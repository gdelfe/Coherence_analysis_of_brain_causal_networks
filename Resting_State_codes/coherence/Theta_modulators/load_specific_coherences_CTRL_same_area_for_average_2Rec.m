
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load only specific sessions for the computation of the average coherence.
% This code creates the .mat files for the mr, ms, sr coherence relative to
% specific sessions (e.g. with no artifacts, no bad channels)
% 
% This is the version for the beta band 
%
%    @ Gino Del Ferraro, May 2021, Pesaran lab, NYU
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

name_struct_input_1 = '/sess_data_lfp_coherence_fk_200_W_5_rec001.mat';
name_struct_input_2 = '/sess_data_lfp_coherence_fk_200_W_5_rec002.mat';


filename = '_rec001_002_no_bad_sess.mat'; % -- filename for sess_data_info.mat
recording1 = 'rec001';
recording2 = 'rec002';
save_dir = 'rec001_002_no_bad_sessions';

freq_band = 'theta_band';
monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
fk = 200; W = 5;


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


% % beta band excluded sessions - rec 001/002
% excluded_sess = [14,30,41];
% excluded_idx = [2,8,11];
% sess_list = 1:size(sess_info{1},1);
% sess_list(excluded_idx) = [];


% beta band excluded sessions - rec 001/002
% excluded_sess = [14,16,22,30,41];
excluded_idx = [2,5,8,9];
sess_list = 1:size(sess_info{1},1);
sess_list(excluded_idx) = [];


cnt_sr = 1; % counter sender-receiver coherencies
cnt_el = 1; % counter for how many modulators excluding the receivers modulators

for i = sess_list % 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
%     if sess_info{2}{i} == '180702' % -- if RS is available 
    if any([4,6,7,10,11] == i) % -- if RS is available 
        dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_same_area/%s',Sess,recording2));
        load(strcat(dir_Modulators,name_struct_input_2)); % RS LFP split into 1 sec window and artifacts removed
    elseif any([1,3,12] == i) % -- use first recording 
        dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_same_area/%s',Sess,recording1));
        load(strcat(dir_Modulators,name_struct_input_1)); % RS LFP split into 1 sec window and artifacts removed
    end 
    
%     
    % -- movie
%     dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_same_area/%s',Sess,recording));
%     load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed


    
    % -- store coherence values sender-receiver and spectrums
    stim(cnt_sr).c_sr = sess_control_lfp.c_sr; % assign S-R coherence value
    stim(cnt_sr).s_s = sess_control_lfp.s_s; % assign sender spectrum
    stim(cnt_sr).s_r = sess_control_lfp.s_r; % receiver spectrum
    cnt_sr = cnt_sr + 1;    % sender/receiver counter
    
    
    mod_Ch = sess_control_lfp.ctrl_idx; % -- modulators (not controls!) index
    
    
    cnt_m = 1;
    for Ch = mod_Ch % for all the modulators in the session
        
        if Ch ~= sess_control_lfp.receiver_idx            
            
%              if Sess == 29 && Ch == 31   % Exclude bad channel
%                 % do nothing
%              else
                % -- structure assignements
                mod(cnt_el).c_ms = sess_control_lfp.ctrl(cnt_m).c_ms ; % assign M-S coherence value for this modulator
                mod(cnt_el).c_mr = sess_control_lfp.ctrl(cnt_m).c_mr;  % M-R coherence
                mod(cnt_el).s_m = sess_control_lfp.ctrl(cnt_m).S_m; % Modulator spectrum
                
                cnt_el = cnt_el + 1; % total modulators counter
%             end
        end
        cnt_m = cnt_m + 1; % counter for modulators within this session
    end 
    
end
 

dir_Mod_results = strcat(dir_RS_Theta,sprintf('/Modulators_Controls_avg_results/%s',save_dir));
if ~exist(dir_Mod_results, 'dir')
    mkdir(dir_Mod_results)
end


% Save coherence and spectrum data in structure format
save(strcat(dir_Mod_results,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s',fk,W,filename)),'mod');
save(strcat(dir_Mod_results,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s',fk,W,filename)),'stim');




