
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs the permutation test for the coherence: it shuffles 
% the trail of one of the two channels when computing the MR, SR, MS
% coherence
%
%    @ Gino Del Ferraro, December 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
step = 110;

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

name_struct_input = '/sess_data_lfp.mat';

iter = 1000; % number of iteration for the permutation test
cnt_sr = 1; % counter sender-receiver coherencies
cnt_el = 1; % counter for how many modulators excluding the receivers modulators
list_sess = 1:19;
list_sess(17) = []; % -- Session 17 and 20 are full of artifacts

for i = list_sess %1:size(sess_info{1},1)-1  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));
    
    load(strcat(dir_Sess,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    %     % ---  time parameter
    tot_time = 150001;
    
    outliers_SR = [sess_data_lfp.outliers_S, sess_data_lfp.outliers_R];
    outliers_SR = unique(outliers_SR);  % -- remove repeated entries in outliers
    
    % %%%%%%%%%%% Sender and Receiver LFP %%%%%%%%%%%%%%%%%%%
    lfp_S = sess_data_lfp.lfp_S;
    lfp_R = sess_data_lfp.lfp_R;
    
    
    % -- remove outliers from sender and receiver
    lfp_S(outliers_SR,:) = [];
    lfp_R(outliers_SR,:) = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- COHERENCE- Sender-Receiver -- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ---- parameters for the coherence-gram
    nt = tot_time;
    fs = 1000;
    fk = 200;
    pad = 2;
    N = 1;
    W = 5;
    % --- coherence
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Coherency - Permutation test for sender and receiver %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    display(['Permutation test sender-receiver coherence...'])

    for j = 1:iter
        
        perm = randperm(size(lfp_R,1));
        % -- coherence calculation via coherency()
        [c_sr,f] = coherency(lfp_S,lfp_R(perm,:),[N W],fs,fk,pad,0.05,1,1);
        coh_sr(cnt_sr).perm(j).c_sr = c_sr;
        if mod(j,100) == 0
            display(['Iteration S-R ---- # ',num2str(j)]);
        end
        
    end
    cnt_sr = cnt_sr + 1;
    
    
%     
%     figure;
%     plot(f,abs(permuted(10).c_sr));
%     hold on
%     plot(f,abs(permuted(100).c_sr))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- MODULATORS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     mod_Ch = sess_data_lfp.mod_idx; % -- modulators (not controls!) index
%     
%     display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
%     
%     
%     % %%%%%%% ALL Electrodes LFP %%%%%%%%%%%%%%%%%%%%%
%     lfp_E_all = sess_data_lfp.lfp_E;
%     
%     cnt_m = 1;
%     for Ch = mod_Ch % for all the modulators in the session
%         
%         close all
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % --- COHERENCE- Modulator - Sender/Receiver -- %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         if Ch ~= sess_data_lfp.receiver_idx % if the electrode is not the receiver itself
%             
%             % -- remove outliers from modulator, sender, and receiver
%             
%             
%             lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
%             outliers_tot = sess_data_lfp.outliers_tot(cnt_m).idx;  % -- get the M,R,S shared outliers
%             
%             % -- Sender and Receiver LFP
%             lfp_S = sess_data_lfp.lfp_S;
%             lfp_R = sess_data_lfp.lfp_R;
%             % -- remove outliers from sender, receiver, and control
%             lfp_S(outliers_tot,:) = [];
%             lfp_R(outliers_tot,:) = [];
%             lfp_E(outliers_tot,:) = [];
%             
%             sess_data_lfp.lfp_E_clean(cnt_m).lfp = lfp_E;   % -- save to structure
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % -- Permutation test for the coherence MR, MR
%             display(['Computing modulator-sender and modulator-receiver coherence...'])
%             for j = 1:iter
%                 
%                 perm = randperm(size(lfp_R,1));
%                 % -- coherence for modulator-sender, modulator-receiver
%                 [c_ms,f] = coherency(lfp_E(perm,:),lfp_S,[N W],fs,fk,pad,0.05,1,1);
%                 [c_mr,f] = coherency(lfp_E(perm,:),lfp_R,[N W],fs,fk,pad,0.05,1,1);
%                 
%                 % -- structure assignements
%                 coh(cnt_el).perm(j).c_ms = c_ms ; % assign M-S coherence value for this modulator
%                 coh(cnt_el).perm(j).c_mr = c_mr;  % M-R coherence
%                 
%                 if mod(j,100) == 0
%                     display(['iteration M-R, M-S ------ # ',num2str(j)]);
%                 end
%                 
%             end
%                 
%       
%             cnt_el = cnt_el + 1; % total modulators counter         
%             cnt_m = cnt_m + 1; % counter for modulators within this session
%             
%         end
%     end
end

keyboard

dir_Perm = strcat(dir_RS,'/Permutation_test2');
if ~exist(dir_Perm, 'dir')
    mkdir(dir_Perm)
end

save(strcat(dir_Perm,sprintf('/coh_spec_m_fk_%d_W_%d.mat',fk,W)),'coh'); % M-S, M-R coherence permuted
save(strcat(dir_Perm,sprintf('/coh_sr_permuted_fk_%d_W_%d.mat',fk,W)),'coh_sr'); % S-R coherence permuted 





