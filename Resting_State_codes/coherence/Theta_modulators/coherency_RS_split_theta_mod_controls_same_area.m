
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the CONTROL
% THETA modulators found by Shaoyu's and both the sender and the receiver
% for the control electrodes (for all the non-modulators electrodes)
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators
%
% In contrast to other codes that employes the coherence-gram to estimate
% the coherence vs frequency, this code employes directly coherency.m
%
% INPUT: sess_control_lfp.mat
%        structure containing all controls and modulator infos + RS LFP split
%
% OUTPUT: txt files with the values of the coherence MR, MS, SR and
% corresponding figures
%
%    @ Gino Del Ferraro, November 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

name_struct_input = '/session_controls_same_area_lfp.mat'; % -- name file to load 
filename = '.mat'; % -- write out filename for sess_data_info.mat 
recording = 'last_recording'; % -- write out directory

freq_band = 'beta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% -- exclude bad sessions 
% excluded_sess = [8,22,30,31];
% excluded_idx = [2,5,8,9];
% sess_list = 1:size(sess_info{1},1);
% sess_list(excluded_idx) = [];


cnt_sr = 1; % counter sender-receiver coherencies
cnt_el = 1; % counter for how many modulators excluding the receivers modulators


for i = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess_Ctrl = strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_same_area',Sess));
    
    load(strcat(dir_Sess_Ctrl,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    dir_Sess_Ctrl = strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_same_area/%s',Sess,recording));
    if ~exist(dir_Sess_Ctrl, 'dir')
        mkdir(dir_Sess_Ctrl)
    end
    
    % -- load list electrodes, sender, receiver
    %     electrode = dataG.RecordPair; % ---- all electrode pairs
    %     receiver = dataG.receiver;  % ---- receiver pair
    %     sender = dataG.sender; % ---- sender pair
    %
    %     % ---  time parameter
    tot_time = 150001;
    %     % ---  freq parameter for the masking
    %     fmin = 10;
    %     fmax = 40;
    
    outliers_SR = [sess_control_lfp.outliers_S, sess_control_lfp.outliers_R];
    outliers_SR = unique(outliers_SR);  % -- remove repeated entries in outliers
    
    % %%%%%%%%%%% Sender and Receiver LFP %%%%%%%%%%%%%%%%%%%
    lfp_S = sess_control_lfp.lfp_S;
    lfp_R = sess_control_lfp.lfp_R;
    
    
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
    % % Coherency for SPLIT LFP data   %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    display(['Computing sender-receiver coherence...'])
    % -- coherence calculation via coherency()
    [c_sr,f,S_s,S_r] = coherency(lfp_S,lfp_R,[N W],fs,fk,pad,0.05,1,1);
    
    
    %     [c_sr_ns,f_ns,S_s_ns,S_r_ns] = coherency(lfp_S_ns,lfp_R_ns,[floor(tot_time/1000 W],fs,fk,pad,0.05,1,1);
    %
    %     figure;
    %     plot(f,abs(c_sr));
    %     hold on
    %     plot(f_ns,abs(c_sr_ns));
    
    
    % -- store coherence values sender-receiver and spectrums
    stim(cnt_sr).c_sr = c_sr; % assign S-R coherence value
    stim(cnt_sr).s_s = S_s; % assign sender spectrum
    stim(cnt_sr).s_r = S_r; % receiver spectrum
    cnt_sr = cnt_sr + 1;    % sender/receiver counter
    
     % -- store coherence results into lfp structure 
    sess_control_lfp.c_sr = c_sr;
    sess_control_lfp.s_s = S_s;
    sess_control_lfp.s_r = S_r;
    
    
    ctrl_Ch = sess_control_lfp.ctrl_idx % control channel indexes
    mod_Ch = sess_control_lfp.mod_idx; % modulators channel indexes
    
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch),'  -- contols mod Ch: ',num2str(ctrl_Ch)])
    

    % %%%%%%% ALL Electrodes LFP %%%%%%%%%%%%%%%%%%%%%
    lfp_E_all = sess_control_lfp.lfp_E;
    
    cnt_m = 1;
    for Ch = ctrl_Ch % for all the modulators in the session
        
        close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- COHERENCE- Modulator - Sender/Receiver -- %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if Ch ~= sess_control_lfp.receiver_idx % if the electrode is not the receiver itself
            
            
            % -- MODULATOR - SENDER
            % -- modulaor lfp
            lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
            outliers_E = sess_control_lfp.outliers_E(cnt_m).idx;  % -- get the modulator's outliers
            % -- Sender  LFP
            lfp_S = sess_control_lfp.lfp_S;
            outliers_S = sess_control_lfp.outliers_S;
            % -- outliers
            outliers_ES = [outliers_E, outliers_S];
            outliers_ES = unique(outliers_ES);
            % -- remove outliers from sender and modulator
            lfp_S(outliers_ES,:) = [];
            lfp_E(outliers_ES,:) = [];
            
            sess_control_lfp.lfp_E_clean(cnt_m).lfp = lfp_E;   % -- save to structure
            
            % -- coherence for modulator-sender, modulator-receiver
            display(['Computing modulator-sender coherence...'])
            [c_ms,f,S_m,S_s] = coherency(lfp_E,lfp_S,[N W],fs,fk,pad,0.05,1,1);
           
            
            
            
            % -- MODULATOR - RECEIVER
            % -- modulaor lfp
            lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
            outliers_E = sess_control_lfp.outliers_E(cnt_m).idx;  % -- get the modulator's outliers
            % -- Receiver  LFP
            lfp_R = sess_control_lfp.lfp_R;
            outliers_R = sess_control_lfp.outliers_R;
            % -- outliers
            outliers_ER = [outliers_E, outliers_R];
            outliers_ER = unique(outliers_ER);
            % -- remove outliers from sender and modulator
            lfp_R(outliers_ER,:) = [];
            lfp_E(outliers_ER,:) = [];
            
            
            display(['Computing modulator-receiver coherence...'])
            [c_mr,f,S_m,S_r] = coherency(lfp_E,lfp_R,[N W],fs,fk,pad,0.05,1,1);
         
 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ABS COHERENCE                 %%%%%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % --- FIGURE --------- %%
            % -- Coherence vs frequency --- %
            fig = figure;
            plot(f,abs(c_sr))
            hold on
            plot(f,abs(c_ms))
            hold on
            plot(f,abs(c_mr))
            grid on
            title(sprintf('Abs coherence vs frequency CONTROLS, ch = %d, causal mod',Ch),'FontSize',10);
            legend('S-R coherence','M-S coherence','M-R coherence')
            %         xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Sess_Ctrl,sprintf('/coherency_vs_freq_ch_%d_fk_%d.jpg',Ch,fk));
            saveas(fig,fname);
            
            % -- structure assignements
            mod(cnt_el).c_ms = c_ms ; % assign M-S coherence value for this modulator
            mod(cnt_el).c_mr = c_mr;  % M-R coherence
            mod(cnt_el).s_m = S_m; % Modulator spectrum
            
            
            sess_control_lfp.ctrl(cnt_m).c_ms = c_ms;
            sess_control_lfp.ctrl(cnt_m).c_mr = c_mr;
            sess_control_lfp.ctrl(cnt_m).S_m = S_m;
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   FIGURES     %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % --- full length with artifacts
%             lfp_S_rshape = reshape(sess_control_lfp.lfp_S',[],1)';
%             lfp_R_rshape = reshape(sess_control_lfp.lfp_R',[],1)';
%             lfp_E_rshape = reshape(sq(sess_control_lfp.lfp_E(Ch,:,:))',[],1)';
%             
%             fig = figure;
%             plot(lfp_S_rshape)
%             hold on
% %             plot(lfp_R_rshape)
% %             hold on
%             plot(lfp_E_rshape)
%             grid on
%             title(sprintf('full length - Controls modulator %d',Ch),'FontSize',11)
%             legend('Sender','Receiver','Modulator')
%             set(gcf, 'Position',  [100, 600, 1000, 500])
            
            %                 fig_name = strcat(dir_Ctrl_all,sprintf('/LFP_Controls_S-R-M_full_length_mod_%d.fig',Ch));
            %                 saveas(fig,fig_name);
            %                 fig_name = strcat(dir_Ctrl_all,sprintf('/LFP_Controls_S-R-M_full_length_mod_%d.png',Ch));
            %                 saveas(fig,fig_name);
            
            % -- full length without artifacts
%             lfp_S_rshape = reshape(lfp_S',[],1)';
%             lfp_R_rshape = reshape(lfp_R',[],1)';
%             lfp_E_rshape = reshape(lfp_E',[],1)';
%             
%             fig = figure;
%             plot(lfp_S_rshape)
%             hold on
%             plot(lfp_R_rshape)
%             hold on
%             plot(lfp_E_rshape)
%             grid on
%             title('Cleaned version LFP Controls ','FontSize',11)
%             legend('Sender','Receiver','Modulator')
%             set(gcf, 'Position',  [100, 600, 1000, 500])
            
            %                 fig_name = strcat(dir_Ctrl_all,sprintf('/LFP_Controls_S-R-M_cleaned_version_no-artifacts_%d.fig',Ch));
            %                 saveas(fig,fig_name);
            %                 fig_name = strcat(dir_Ctrl_all,sprintf('/LFP_Controls_S-R-M_cleaned_version_no-artifacts_%d.png',Ch));
            %                 saveas(fig,fig_name);
            
            cnt_el = cnt_el + 1; % total control counter            
        end
        cnt_m = cnt_m + 1; % counter for control within this session

    end    
        save(strcat(dir_Sess_Ctrl,sprintf('/sess_data_lfp_coherence_fk_%d_W_%d%s',fk,W,filename)),'sess_control_lfp');

end


% dir_Ctrl = strcat(dir_RS_Theta,sprintf('/Modulators_Controls_avg_results/%s',recording));
% if ~exist(dir_Ctrl, 'dir')
%     mkdir(dir_Ctrl)
% end   
% 
% % Save coherence and spectrum data in structure format
% save(strcat(dir_Ctrl,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s',fk,W,filename)),'mod');
% save(strcat(dir_Ctrl,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s',fk,W,filename)),'stim');


% 
% keyboard
% 
% 
% % -- load structure files
% fk = 200; W = 5;
% load(strcat(dir_RS,sprintf('/coh_spec_m_all_Controls_same_area_fk_%d_W_%d.mat',fk,W)))
% load(strcat(dir_RS,sprintf('/coh_spec_sr_all_Controls_same_area_fk_%d_W_%d.mat',fk,W)))
% 
% % -- structures to matrices
% mod_mat = cell2mat(struct2cell(mod)); % transform struct to mat for modulators
% stim_mat = cell2mat(struct2cell(stim)); % transform struct to mat for sender-receiver
% 
% 
% % -- assign fields to matrices
% coh_ms = sq(mod_mat(1,:,:))'; % 1st field, c_ms
% coh_mr = sq(mod_mat(2,:,:))'; %  2nd field, c_mr
% spec_m = sq(mod_mat(3,:,:))'; %  3rd field, spec_m
% 
% coh_sr = sq(stim_mat(1,:,:))'; % 1st field, c_sr
% spec_s = sq(stim_mat(2,:,:))'; %  2nd field, spec_s
% spec_r = sq(stim_mat(3,:,:))'; %  3rd field, spec_r
% 
% 
% 
% % --- mean coherences
% mean_cho_ms = mean(abs(coh_ms));  % modulator - sender
% mean_cho_mr = mean(abs(coh_mr));  % modulator - receiver
% mean_cho_sr = mean(abs(coh_sr));  % sender - receiver
% 
% % --- std coherences
% std_cho_ms = std(abs(coh_ms));  % modulator - sender
% std_cho_mr = std(abs(coh_mr)); % modulator - receiver
% std_cho_sr = std(abs(coh_sr));  % modulator - receiver
% 
% % --- Error bars
% M = size(coh_ms,1);
% S = size(coh_sr,1);
% err_ms = std_cho_ms/sqrt(M);
% err_mr = std_cho_mr/sqrt(M);
% err_sr = std_cho_sr/sqrt(S);
% 
% 
% 
% set(0,'DefaultFigureVisible','on')
% % -- FIGURE: Plot average coherence across sessions for MR, SR, MS
% fig = figure;
% % hAx=axes;
% % hAx.XScale='linear'
% % hAx.YScale='log'
% hold all
% 
% % shadedErrorBar(f,mean_cho_ms,err_ms,'lineProps','b','patchSaturation',0.4); hold on
% % shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[255, 83, 26]/255},'patchSaturation',0.4); hold on
% % shadedErrorBar(f,mean_cho_sr,err_sr,'lineprops',{'color',[230, 184 , 0]/255},'patchSaturation',0.4); hold on
% 
% 
% shadedErrorBar(f,mean_cho_ms,err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_sr,err_sr,'lineprops',{'color',[0 204 204]/255},'patchSaturation',0.4); hold on
% 
% grid on
% title('Abs coherency CONTROLS of MS, SR, SR - Resting State, all Sessions','FontSize',11);
% xlabel('freq (Hz)');
% ylabel('coherence');
% % legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% % legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% legend('M-S abs coherency','M-R abs coherency','S-R abs coherency','FontSize',10)
% % legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% % xlim([0 60])
% set(gcf, 'Position',  [100, 600, 1000, 600])
% 
% 
% fname = strcat(dir_RS,sprintf('/coherency_mean_all_Controls_same_area_split-data_MS_MR_SR_W_%d_fk_%d-all-Sess.png',W,fk));
% saveas(fig,fname)
% fname = strcat(dir_RS,sprintf('/coherency_mean_all_Controls_same_area_split-data_MS_MR_SR_W_%d_fk_%d-all-Sess.fig',W,fk));
% saveas(fig,fname)
% 
% 
