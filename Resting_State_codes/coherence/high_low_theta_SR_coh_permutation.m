
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code compute the theta coherence between the modulator-sender and
% modulator-receiver for modulators having high/low theta power.
%
%    @ Gino Del Ferraro, March 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
dir_high_low_theta = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

n_iter = 5000;

% ---- parameters for the coherence-gram
tot_time = 150001;
nt = tot_time;
fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;


for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
        
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    dir_Sess_send_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/send_rec/Data',Sess));

    
    load(strcat(dir_Sess_send_rec_data,'/send_rec_coh_for_session.mat'));

        
%     load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info 
    

    % number of modulators in that session 
    n_mod = size(send_rec.mod_idx,2);
    cnt_m = 1; % counter for numb of modulators checked for the outlier trials remover
    
    send_rec_perm.sess_idx = send_rec.sess_idx;
    send_rec_perm.mod_idx = send_rec.mod_idx; 
    
    cnt_m = 1;
    for m = send_rec.mod_idx % for each modulator 
        
        display(['-- Modulator ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);
        send_rec_perm.mod(cnt_m).m = m;

        % matrices to store high/low coherence for permuted data for each modulator
        coh_perm_sr_high = zeros(n_iter,409);
        coh_perm_sr_low = zeros(n_iter,409);
  
       
        low_idx = send_rec.mod(cnt_m).low_pow_idx;
        high_idx = send_rec.mod(cnt_m).high_pow_idx;
        
        lfp_S = send_rec.mod(cnt_m).lfp_S_clean;
        lfp_R = send_rec.mod(cnt_m).lfp_R_clean;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Coherency Modulator-Sender Permutation test
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        low_high_idx = [low_idx; high_idx]; % indexes for the high/low theta wrt sender 
        W = 5;
        
        for iter = 1:n_iter
            
            perm_idx = low_high_idx(randperm(length(low_high_idx)));
            cut = length(perm_idx)/2;
            perm_idx_low = perm_idx(1:cut);
            perm_idx_high = perm_idx(cut+1:end);
            
            % -- coherence calculation via coherency()
            [c_sr_low_perm,f] = coherency(lfp_S(perm_idx_low,:),lfp_R(perm_idx_low,:),[N W],fs,fk,pad,0.05,1,1);
            [c_sr_high_perm,f] = coherency(lfp_S(perm_idx_high,:),lfp_R(perm_idx_high,:),[N W],fs,fk,pad,0.05,1,1);
            
            coh_perm_sr_low(iter,:) = c_sr_low_perm;  % n iteration x frequency
            coh_perm_sr_high(iter,:) = c_sr_high_perm;
            
        end 
        
        
        send_rec_perm.mod(cnt_m).c_ms_low = coh_perm_sr_low;
        send_rec_perm.mod(cnt_m).c_ms_high = coh_perm_sr_high;
        
           
  
        cnt_m = cnt_m +1;
        
    end % end of for cycle for all the modulators in a given session 
    
      
    dir_Sess_send_rec = strcat(dir_high_low_theta,sprintf('/Sess_%d/send_rec/Data',Sess));
    save(strcat(dir_Sess_send_rec,'/mod_send_perm.mat'),'send_rec_perm');
    
    
    clear send_rec_perm 
  
end % end of for cycle for all the sessions 




