
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the THETA
% modulators and both the sender and the receiver
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

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

% addpath('T:/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat 
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

cnt_sr = 1; % counter sender-receiver coherencies
cnt_el = 1; % counter for how many modulators excluding the receivers modulators
% sess_list = [1,3,4,5]

for i = 16 %1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    dir_Mod_recording = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    if ~exist(dir_Mod_recording, 'dir')
        mkdir(dir_Mod_recording)
    end
    %     % ---  time parameter
    tot_time = 150001;
    
    outliers_SR = [sess_data_lfp.outliers_S, sess_data_lfp.outliers_R];
    outliers_SR = unique(outliers_SR)  % -- remove repeated entries in outliers
    
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
    % % Coherency for SPLIT LFP        %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    display(['Computing sender-receiver coherence...'])
    % -- coherence calculation via coherency()
    [c_sr,f,S_s,S_r] = coherency(lfp_S,lfp_R,[N W],fs,fk,pad,0.05,1,1);
    
    
    %     [c_sr_ns,f_ns,S_s_ns,S_r_ns] = coherency(lfp_S_ns,lfp_R_ns,[floor(tot_time/1000 W],fs,fk,pad,0.05,1,1);
    %
        figure;
        plot(f,abs(c_sr));
%         hold on
%         plot(f_ns,abs(c_sr_ns));
    
    
    % -- store coherence values sender-receiver and spectrums
    stim(cnt_sr).c_sr = c_sr; % assign S-R coherence value
    stim(cnt_sr).s_s = S_s; % assign sender spectrum
    stim(cnt_sr).s_r = S_r; % receiver spectrum
    cnt_sr = cnt_sr + 1;    % sender/receiver counter
    
    % -- store coherence results into lfp structure 
    sess_data_lfp.c_sr = c_sr;
    sess_data_lfp.s_s = S_s;
    sess_data_lfp.s_r = S_r;
    
    
    mod_Ch = sess_data_lfp.mod_idx; % -- modulators (not controls!) index
    
    % --- Remove channels with artifacts
    if Sess == 19
        mod_Ch(mod_Ch == 60) = []; % remove channel 
    end
    if Sess == 41
        mod_Ch(mod_Ch == 8) = []; % remove channel 
    end
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
    
    % %%%%%%% ALL Electrodes LFP %%%%%%%%%%%%%%%%%%%%%
    lfp_E_all = sess_data_lfp.lfp_E;
    
    cnt_m = 1;
    for Ch = mod_Ch % for all the modulators in the session
        
        close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- COHERENCE- Modulator - Sender/Receiver -- %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if Ch ~= sess_data_lfp.receiver_idx % if the electrode is not the receiver itself
            
            % -- remove outliers from modulator, sender, and receiver
            
            % -- MODULATOR - SENDER 
            % -- modulaor lfp
            lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
            outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;  % -- get the modulator's outliers
            % -- Sender  LFP
            lfp_S = sess_data_lfp.lfp_S;
            outliers_S = sess_data_lfp.outliers_S;
            % -- outliers
            outliers_ES = [outliers_E, outliers_S];
            outliers_ES = unique(outliers_ES);
            % -- remove outliers from sender and modulator 
            lfp_S(outliers_ES,:) = [];
            lfp_E(outliers_ES,:) = [];
            
            sess_data_lfp.lfp_E_clean(cnt_m).lfp = lfp_E;   % -- save to structure
            
            % -- coherence for modulator-sender, modulator-receiver
            display(['Computing modulator-sender coherence...'])
            [c_ms,f,S_m,S_s] = coherency(lfp_E,lfp_S,[N W],fs,fk,pad,0.05,1,1);
                        
            
            % -- MODULATOR - RECEIVER
            % -- modulaor lfp
            lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
            outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;  % -- get the modulator's outliers
            % -- Receiver  LFP
            lfp_R = sess_data_lfp.lfp_R;
            outliers_R = sess_data_lfp.outliers_R;
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
            title(sprintf('Abs coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
            legend('S-R coherence','M-S coherence','M-R coherence')
            %         xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Mod_recording,sprintf('/coherency_vs_freq_ch_%d_fk_%d.jpg',Ch,fk));
            saveas(fig,fname);
            
            % -- structure assignements
            mod(cnt_el).c_ms = c_ms ; % assign M-S coherence value for this modulator
            mod(cnt_el).c_mr = c_mr;  % M-R coherence
            mod(cnt_el).s_m = S_m; % Modulator spectrum
            
            sess_data_lfp.mod(cnt_m).c_ms = c_ms;
            sess_data_lfp.mod(cnt_m).c_mr = c_mr;
            sess_data_lfp.mod(cnt_m).S_m = S_m;

            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   FIGURES     %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % --- full length with artifacts
            lfp_S_rshape = reshape(sess_data_lfp.lfp_S',[],1)';
            lfp_R_rshape = reshape(sess_data_lfp.lfp_R',[],1)';
            lfp_E_rshape = reshape(sq(sess_data_lfp.lfp_E(Ch,:,:))',[],1)';
            
            fig = figure;
            plot(lfp_S_rshape)
            hold on
            plot(lfp_R_rshape)
            hold on
            plot(lfp_E_rshape)
            grid on
            title(sprintf('full length - modulator %d',Ch),'FontSize',11)
            legend('Sender','Receiver','Modulator')
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fig_name = strcat(dir_Mod_recording,sprintf('/LFP_S-R-M_full_length_mod_%d.fig',Ch));
            saveas(fig,fig_name);
            fig_name = strcat(dir_Mod_recording,sprintf('/LFP_S-R-M_full_length_mod_%d.png',Ch));
            saveas(fig,fig_name);
            
              % -- full length without artifacts
            lfp_S_rshape = reshape(lfp_S',[],1)';
            lfp_R_rshape = reshape(lfp_R',[],1)';
            lfp_E_rshape = reshape(lfp_E',[],1)';
            
            fig = figure;
            plot(lfp_S_rshape)
            hold on
            plot(lfp_R_rshape)
            hold on
            plot(lfp_E_rshape)
            grid on
            title('Cleaned version LFP ','FontSize',11)
            legend('Sender','Receiver','Modulator')
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fig_name = strcat(dir_Mod_recording,sprintf('/LFP_S-R-M_cleaned_version_no-artifacts_%d.fig',Ch));
            saveas(fig,fig_name);
            fig_name = strcat(dir_Mod_recording,sprintf('/LFP_S-R-M_cleaned_version_no-artifacts_%d.png',Ch));
            saveas(fig,fig_name);
            
            cnt_el = cnt_el + 1; % total modulators counter                     
        end
        cnt_m = cnt_m + 1; % counter for modulators within this session
    end
    save(strcat(dir_Mod_recording,sprintf('/sess_data_lfp_coherence_fk_%d_W_%d%s',fk,W,filename)),'sess_data_lfp');
end

keyboard 

% dir_Mod_results = strcat(dir_RS_Theta,sprintf('/Modulators_Controls_avg_results/%s',recording));
% if ~exist(dir_Mod_results, 'dir')
%     mkdir(dir_Mod_results)
% end    
% 
% 
% % Save coherence and spectrum data in structure format
% save(strcat(dir_Mod_results,sprintf('/coh_spec_m_fk_%d_W_%d%s',fk,W,filename)),'mod');
% save(strcat(dir_Mod_results,sprintf('/coh_spec_sr_fk_%d_W_%d%s',fk,W,filename)),'stim');
% 
% keyboard
% % 
% % -- load structure files
% fk = 200;
% load(strcat(dir_Mod_results,sprintf('/coh_spec_m_fk_%d_W_%d%s',fk,W,filename)))
% load(strcat(dir_Mod_results,sprintf('/coh_spec_sr_fk_%d_W_%d%s',fk,W,filename)))
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
% % --- mean coherences
% mean_cho_ms = mean(abs(coh_ms));  % modulator - sender
% mean_cho_mr = mean(abs(coh_mr));  % modulator - receiver
% mean_cho_sr = mean(abs(coh_sr));  % sender - receiver
% 
% % --- std coherences
% std_cho_ms = std(abs(coh_ms));  % modulator - sender
% std_cho_mr = std(abs(coh_mr)); % modulator - receiver
% std_cho_sr = std(abs(coh_sr));  % modulator - receiver
% % 
% % --- Error bars
% M = size(coh_ms,1);
% S = size(coh_sr,1);
% err_ms = std_cho_ms/sqrt(M);
% err_mr = std_cho_mr/sqrt(M);
% err_sr = std_cho_sr/sqrt(S);
% % 
% % 
% % 
% set(0,'DefaultFigureVisible','on')
% % -- FIGURE: Plot average coherence across sessions for MR, SR, MS
% fig = figure;
% % hAx=axes;
% % hAx.XScale='linear'
% % hAx.YScale='log'
% hold all
% 
% 
% shadedErrorBar(f,mean_cho_ms,err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_sr,err_sr,'lineprops',{'color',[0 204 204]/255},'patchSaturation',0.4); hold on
% 
% grid on
% title('Mean abs coherency Theta Modulators - MS, SR, SR - Resting State','FontSize',11);
% xlabel('freq (Hz)');
% ylabel('coherence');
% % legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% % legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% legend('M-S coherency','M-R coherency','S-R coherency','FontSize',10)
% % legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% % xlim([0 60])
% set(gcf, 'Position',  [100, 600, 1000, 600])
% % 
% % fname = strcat(dir_Modulators,sprintf('/coherency_theta_mod_MS_MR_SR_W_%d_fk_%d.png',W,fk));
% % saveas(fig,fname)
% % fname = strcat(dir_Modulators,sprintf('/coherency_theta_mod_MS_MR_SR_W_%d_fk_%d.fig',W,fk));
% % saveas(fig,fname)
% % 
% % 
