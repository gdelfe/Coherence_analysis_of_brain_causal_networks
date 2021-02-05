
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the
% modulators found by Shaoyu's and the receiver only, for comparison with
% the STIM analogous
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators
%
% NOTE: The n_trials of the lfp in the RS are chosen equal to the n_trials
% available in the STIM exp in order to have analogous statistics
%
% INPUT: sess_data_lfp.mat
%        structure containing all modulator infos + RS LFP split
%
% OUTPUT: txt files with the values of the coherence MR, MS, SR and
% corresponding figures
%
%    @ Gino Del Ferraro, January 2021, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Stim = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';

fk = 200;
W = 5;

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

name_struct_input = '/sess_data_lfp.mat';
load(strcat(dir_Stim,sprintf('/coh_spec_mr_sameRStrails_fk_%d_W_%d.mat',fk,W)));


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
    
    % ---- parameters for coherence
    nt = tot_time;
    fs = 1000;
    fk = 200;
    pad = 2;
    N = 1;
    W = 5;
    % --- coherence

    
    mod_Ch = sess_data_lfp.mod_idx; % -- modulators (not controls!) index
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
    
    dir_Modulators = strcat(dir_Sess,'/Modulators_same_trials_RS_STIM');
    if ~exist(dir_Modulators, 'dir')
        mkdir(dir_Modulators)
    end
    
    % %%%%%%% ALL Electrodes LFP %%%%%%%%%%%%%%%%%%%%%
    lfp_E_all = sess_data_lfp.lfp_E;
    
    cnt_m = 1;
    for Ch = mod_Ch % for all the modulators in the session
        
        close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- COHERENCE- Modulator - Sender/Receiver -- %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if Ch ~= sess_data_lfp.receiver_idx % if the electrode is not the receiver itself
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Select same number of trials for RS and STIM
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % -- remove outliers from modulator and receiver
            outliers_tot = [sess_data_lfp.outliers_E(cnt_m).idx, sess_data_lfp.outliers_R];  % -- get the M,R,S shared outliers
            outliers_tot = unique(outliers_tot);
            
            %%%% CHECK HERE: stim needs to be loaded
            n_trials_STIM = stim(cnt_el).n_trials_RS_STIM; % tot numbero of trials STIM
            n_trials_RS = (150 - size(outliers_tot,2)); % tot numb trials RS
            
            % --- Select STIM trials to match the number of RS trial for idential statistics
            n_trials = min(n_trials_RS,n_trials_STIM) % select the min numb of trial between RS and STIM 
            % remove outliers from RS, in order not to include them 
            RS_trials = 1:150;
            RS_trials(outliers_tot)=[];
            
            trials = RS_trials(1:n_trials); % 
            
            % -- Modulator and Receiver LFP
            lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
            lfp_R = sess_data_lfp.lfp_R;
            
            % include only the trails without artifacts and as many as in STIM
            lfp_E = lfp_E(trials,:);
            lfp_R = lfp_R(trials,:);
            
            
            sess_data_lfp.lfp_E_clean(cnt_m).lfp = lfp_E;   % -- save lfp_E to structure
            
            % -- coherence for modulator-receiver            
            display(['Computing modulator-receiver coherence...'])
            [c_mr,f,S_m,S_r] = coherency(lfp_E,lfp_R,[N W],fs,fk,pad,0.05,1,1);
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ABS COHERENCE                 %%%%%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % --- FIGURE --------- %%
            % -- Coherence vs frequency --- %
            fig = figure;
            plot(f,abs(c_mr))
            grid on
            title(sprintf('Abs coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
            legend('M-R coherence')
            %         xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Modulators,sprintf('/coherency_MR_vs_freq_ch_%d_fk_%d.jpg',Ch,fk));
            saveas(fig,fname);
            
            % -- structure assignements
            modulator(cnt_el).c_mr = c_mr;  % M-R coherence
            modulator(cnt_el).s_m = S_m; % Modulator spectrum
          
            
            cnt_el = cnt_el + 1; % total modulators counter                     
        end
        sess_data_lfp.lfp_R_clean.lfp = lfp_R;   % -- save lfp_E to structure
        cnt_m = cnt_m + 1; % counter for modulators within this session

    end
end

dir_coherence = strcat(dir_RS,'/Coherence_STIM_RS_same_trials');
if ~exist(dir_coherence, 'dir')
    mkdir(dir_coherence)
end
    
    
keyboard
% Save coherence and spectrum data in structure format
save(strcat(dir_coherence,sprintf('/coh_spec_MR_only_for_STIM_comparison_fk_%d_W_%d.mat',fk,W)),'mod');

% -- load structure files
fk = 200;
load(strcat(dir_coherence,sprintf('/coh_spec_m_fk_%d_W_%d.mat',fk,W)))

% -- structures to matrices
mod_mat = cell2mat(struct2cell(mod)); % transform struct to mat for modulators
stim_mat = cell2mat(struct2cell(stim)); % transform struct to mat for sender-receiver


% -- assign fields to matrices
coh_ms = sq(mod_mat(1,:,:))'; % 1st field, c_ms
coh_mr = sq(mod_mat(2,:,:))'; %  2nd field, c_mr
spec_m = sq(mod_mat(3,:,:))'; %  3rd field, spec_m

coh_sr = sq(stim_mat(1,:,:))'; % 1st field, c_sr
spec_s = sq(stim_mat(2,:,:))'; %  2nd field, spec_s
spec_r = sq(stim_mat(3,:,:))'; %  3rd field, spec_r



% --- mean coherences
mean_cho_ms = mean(abs(coh_ms));  % modulator - sender
mean_cho_mr = mean(abs(coh_mr));  % modulator - receiver
mean_cho_sr = mean(abs(coh_sr));  % sender - receiver

% --- std coherences
std_cho_ms = std(abs(coh_ms));  % modulator - sender
std_cho_mr = std(abs(coh_mr)); % modulator - receiver
std_cho_sr = std(abs(coh_sr));  % modulator - receiver

% --- Error bars
M = size(coh_ms,1);
S = size(coh_sr,1);
err_ms = std_cho_ms/sqrt(M);
err_mr = std_cho_mr/sqrt(M);
err_sr = std_cho_sr/sqrt(S);



set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, SR, MS
fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all

% shadedErrorBar(f,mean_cho_ms,err_ms,'lineProps','b','patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[255, 83, 26]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_sr,err_sr,'lineprops',{'color',[230, 184 , 0]/255},'patchSaturation',0.4); hold on


shadedErrorBar(f,mean_cho_ms,err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_sr,err_sr,'lineprops',{'color',[0 204 204]/255},'patchSaturation',0.4); hold on

grid on
title('Abs coherency of MS, SR, SR - Resting State, all Sessions','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('M-S abs coherency','M-R abs coherency','S-R abs coherency','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% xlim([0 60])
set(gcf, 'Position',  [100, 600, 1000, 600])

fname = strcat(dir_RS,sprintf('/coherency_mean_split-data_MS_MR_SR_W_%d_fk_%d-all-Sess.png',W,fk));
saveas(fig,fname)
fname = strcat(dir_RS,sprintf('/coherency_mean_split-data_MS_MR_SR_W_%d_fk_%d-all-Sess.fig',W,fk));
saveas(fig,fname)


