
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the theta coherence between the modulator-sender and
% modulator-receiver for modulators having high/low theta power and
% computes the null distribution of zero-coherence hypothesis by 
% performing permutation test 
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

num_iter = 1000; % number of iteration for permutation test

% ---- parameters for the coherence-gram
% ---  time parameter
tot_time = 150001;
nt = tot_time;
fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;

coh_all_c_ms = [];
coh_all_c_mr = [];

for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    clear mod_rec mod_send coh_all sess_data_lfp
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info 
    
    % outliers time series in sender and receiver 
    outliers_S = sess_data_lfp.outliers_S;
    outliers_R = sess_data_lfp.outliers_R;
    

    % number of modulators in that session 
    n_mod = size(sess_data_lfp.mod_idx,2);
    cnt_m = 1; % counter for numb of modulators checked for the outlier trials remover
    
    % matrices to store high/low coherence wrt sender and receiver 
    coh_ms_sess = zeros(n_mod,num_iter,409);
    coh_mr_sess = zeros(n_mod,num_iter,409);

    
    for m = sess_data_lfp.mod_idx
        
        display(['-- Modulator ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);
        
        
        % When computing MS and MR coherence we have to remove the outliers
        % for the sender/receiver, and this are different time series for
        % the sender and the receiver 
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%
        % In relation to the sender
        
        lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
        lfp_S = sess_data_lfp.lfp_S;
        
        outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;
        
        outliers_ES = [outliers_S, outliers_E];
        outliers_ES = unique(outliers_ES);  % -- remove repeated entries in outliers
        
        lfp_E(outliers_ES,:) = [];
        lfp_S(outliers_ES,:) = [];
        
        n_trial = size(lfp_S,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Coherency Modulator-Sender permutation test 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        display(['-- MS coherence permutation test ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);

        for iter = 1:num_iter
            
            perm_idx = randperm(n_trial);
            % -- coherence calculation via coherency()
            [c_ms,f] = coherency(lfp_S(perm_idx,:),lfp_E,[N W],fs,fk,pad,0.05,1,1);
            
            coh_ms_sess(cnt_m,iter,:) = c_ms;
            
        end
        
        % assign coherence MS values (permuted) to structure 
        coh_sess.ms = coh_ms_sess;
        
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % In relation to the receiver
        
        lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
        lfp_R = sess_data_lfp.lfp_R;
        
        outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;
        
        outliers_ER = [outliers_R, outliers_E];
        outliers_ER = unique(outliers_ER);  % -- remove repeated entries in outliers
        
        lfp_E(outliers_ER,:) = [];
        lfp_R(outliers_ER,:) = [];
        
        n_trial = size(lfp_R,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Coherency Modulator-Receiver permutation test 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        display(['-- MR coherence permutation test ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);
        
        for iter = 1:num_iter
            
            perm_idx = randperm(n_trial);
            % -- coherence calculation via coherency()
            [c_mr,f] = coherency(lfp_R(perm_idx,:),lfp_E,[N W],fs,fk,pad,0.05,1,1);
            
            coh_mr_sess(cnt_m,iter,:) = c_mr;
            
        end

        % assign coherence MR values (permuted) to structure 
        coh_sess.mr = coh_mr_sess;
        
        
        dir_Sess = strcat(dir_high_low_theta,sprintf('/Sess_%d',Sess));
        if ~exist(dir_Sess, 'dir')
            mkdir(dir_Sess)
        end
        
        save(strcat(dir_Sess,'/permuted_coh_all_mod_in_sess'),'coh_sess');
        

        
        cnt_m = cnt_m +1;
        
    end % end of for cycle for all the modulators in a given session
    
    % store MS and MR permuted-coherences for all the session into matrix
    coh_all_c_ms = [coh_all_c_ms; coh_ms_sess];
    coh_all_c_mr = [coh_all_c_mr; coh_mr_sess];
    
    
end % end for all the sessions  
    
save(strcat(dir_high_low_theta,'/coh_all_permuted_c_ms.mat'),'coh_all_c_ms','-v7.3')
save(strcat(dir_high_low_theta,'/coh_all_permuted_c_mr.mat'),'coh_all_c_mr','-v7.3')

% load(strcat(dir_high_low_theta,'/coh_all_permuted_c_ms.mat'));
% load(strcat(dir_high_low_theta,'/coh_all_permuted_c_mr.mat'));


    
