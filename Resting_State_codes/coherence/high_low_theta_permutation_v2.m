
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
    
    
    close all
    clear mod_rec mod_send coh_all sess_data_lfp
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info 

    %     % ---  time parameter
    tot_time = 150001;
    
    % outliers time series in sender and receiver 
    outliers_S = sess_data_lfp.outliers_S;
    outliers_R = sess_data_lfp.outliers_R;
    

    % number of modulators in that session 
    n_mod = size(sess_data_lfp.mod_idx,2);
    cnt_m = 1; % counter for numb of modulators checked for the outlier trials remover
    

    mod_send_perm.mod_idx = sess_data_lfp.mod_idx; 
    mod_rec_perm.mod_idx = sess_data_lfp.mod_idx;
    
    for m = sess_data_lfp.mod_idx
        
        display(['-- Modulator ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);
        
        
        % matrices to store high/low coherence wrt sender and receiver for the
        % permuted data for each modulator
        coh_perm_ms_high = zeros(n_iter,409);
        coh_perm_ms_low = zeros(n_iter,409);
        coh_perm_mr_high = zeros(n_iter,409);
        coh_perm_mr_low = zeros(n_iter,409);
    
    
        % When computing MS and MR coherence we have to remove the outliers
        % for the sender/receiver, and this are different time series for
        % the sender and the receiver 
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%
        % In relation to the sender -- LFPs
        
        lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
        lfp_S = sess_data_lfp.lfp_S;
        
        outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;
        
        outliers_ES = [outliers_S, outliers_E];
        outliers_ES = unique(outliers_ES);  % -- remove repeated entries in outliers
        
        lfp_E(outliers_ES,:) = [];
        lfp_S(outliers_ES,:) = [];
        
        % Compute the spectrum for each trial. Format: iTrial x times
        W = 3;
        [spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200);
        
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
        
        mod_send_perm.mod(cnt_m).m = m;
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Coherency Modulator-Sender Permutation test
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        low_high_idx = [low_idx; high_idx]; % indexes for the high/low theta wrt sender 
        W = 5;
        
        for iter = 1:n_iter
            
            perm_idx = low_high_idx(randperm(length(low_high_idx)));
            perm_idx_low = perm_idx(1:cut);
            perm_idx_high = perm_idx(cut+1:end);
            
            % -- coherence calculation via coherency()
            [c_ms_low_perm,f] = coherency(lfp_S(perm_idx_low,:),lfp_E(perm_idx_low,:),[N W],fs,fk,pad,0.05,1,1);
            [c_ms_high_perm,f] = coherency(lfp_S(perm_idx_high,:),lfp_E(perm_idx_high,:),[N W],fs,fk,pad,0.05,1,1);
            
            coh_perm_ms_low(iter,:) = c_ms_low_perm;
            coh_perm_ms_high(iter,:) = c_ms_high_perm;
            
        end 
        
        
        mod_send_perm.mod(cnt_m).c_ms_low = coh_perm_ms_low;
        mod_send_perm.mod(cnt_m).c_ms_high = coh_perm_ms_high;
        

        dir_Sess_mod_send_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_send/Data',Sess));
        if ~exist(dir_Sess_mod_send_data, 'dir')
            mkdir(dir_Sess_mod_send_data)
        end
               
        save(strcat(dir_Sess_mod_send_data,sprintf('/mod_send_perm.mat',cnt_m)),'mod_send_perm');
        
        
        
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % In relation to the receiver -- LFPs
        
        lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
        lfp_R = sess_data_lfp.lfp_R;
        
        outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;
        
        outliers_ER = [outliers_R, outliers_E];
        outliers_ER = unique(outliers_ER);  % -- remove repeated entries in outliers
        
        lfp_E(outliers_ER,:) = [];
        lfp_R(outliers_ER,:) = [];
        
        
        % Compute the spectrum for each trial. Format: iTrial x times
        W = 3;
        [spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200);
        theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
        theta_pow_mean = mean(theta_pow);
        theta_pow = theta_pow - theta_pow_mean;
        
        [sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
        cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution
        
        low_theta_pow = sort_theta(1:cut);
        low_idx = trial_idx(1:cut);
        
        high_theta = sort_theta(end-cut+1:end);
        high_idx = trial_idx(end-cut+1:end);
      
        mod_rec_perm.mod(cnt_m).m = m;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Coherency Modulator-Receiver Permutation Test 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        low_high_idx = [low_idx; high_idx]; % indexes for the high/low theta wrt sender
        W = 5;
        
        for iter = 1:n_iter
            
            perm_idx = low_high_idx(randperm(length(low_high_idx)));
            perm_idx_low = perm_idx(1:cut);
            perm_idx_high = perm_idx(cut+1:end);
            
            % -- coherence calculation via coherency()
            [c_mr_low_perm,f] = coherency(lfp_R(perm_idx_low,:),lfp_E(perm_idx_low,:),[N W],fs,fk,pad,0.05,1,1);
            [c_mr_high_perm,f] = coherency(lfp_R(perm_idx_high,:),lfp_E(perm_idx_high,:),[N W],fs,fk,pad,0.05,1,1);
            
            coh_perm_mr_low(iter,:) = c_mr_low_perm;
            coh_perm_mr_high(iter,:) = c_mr_high_perm;
            
        end
        
        
        mod_rec_perm.mod(cnt_m).c_mr_low = coh_perm_mr_low;
        mod_rec_perm.mod(cnt_m).c_mr_high = coh_perm_mr_high;
        
        
        dir_Sess_mod_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Data',Sess));
        if ~exist(dir_Sess_mod_rec_data, 'dir')
            mkdir(dir_Sess_mod_rec_data)
        end
        
        save(strcat(dir_Sess_mod_rec_data,sprintf('/mod_rec_perm.mat',cnt_m)),'mod_rec_perm');
        
   
        
        cnt_m = cnt_m +1;
        
    end % end of for cycle for all the modulators in a given session 
  
end % end of for cycle for all the sessions 




