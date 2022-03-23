
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the
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
    
    outliers_S = sess_data_lfp.outliers_S;
    outliers_R = sess_data_lfp.outliers_R;
    
    outliers_SR = [sess_data_lfp.outliers_S, sess_data_lfp.outliers_R];
    outliers_SR = unique(outliers_SR);  % -- remove repeated entries in outliers
    
    % %%%%%%%%%%% Sender and Receiver LFP %%%%%%%%%%%%%%%%%%%
    lfp_S = sess_data_lfp.lfp_S;
    lfp_R = sess_data_lfp.lfp_R;
    
    
    % -- remove outliers from sender and receiver
    lfp_S(outliers_S,:) = [];
    lfp_R(outliers_R,:) = [];
    
    
      % ---- parameters for the coherence-gram
    nt = tot_time;
    fs = 1000;
    fk = 200;
    pad = 2;
    N = 1;
    W = 5;
  
    
    % %%%%%%%%%%% Compute 'normalized' spectrum %%%%%%%%%%%%%%%%%%%

    Xs = fft(lfp_S')';
    norm_Fs = sum(sum(Xs.*conj(Xs),2)); % what is delta f here ?
    
    norm_Ts = sum(sum(lfp_S.*lfp_S,2))*0.001;
    norm_Tr = sum(sum(lfp_R.*lfp_R,2))*0.001;
    
    lfp_S = lfp_S/norm_Ts; 
    lfp_R = lfp_R/norm_Tr;
    
    S_s = double(dmtspec(lfp_S,[1 3],fs,fk, pad, 0.05, 1));
    S_r = double(dmtspec(lfp_R,[1 3],fs,fk, pad, 0.05, 1));
    
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Coherency for SPLIT LFP        %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % %%%%%%%%%%% Sender and Receiver LFP %%%%%%%%%%%%%%%%%%%
    lfp_S = sess_data_lfp.lfp_S;
    lfp_R = sess_data_lfp.lfp_R;
    
    
    % -- remove outliers from sender and receiver
    lfp_S(outliers_SR,:) = [];
    lfp_R(outliers_SR,:) = [];
    
    display(['Computing sender-receiver coherence...'])
    
    % -- coherence Sender-Receiver 
    [c_sr,f] = coherency(lfp_S,lfp_R,[N W],fs,fk,pad,0.05,1,1);

        
   
   
    
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
    
    
    mod_Ch = sess_data_lfp.mod_idx; % -- modulators (not controls!) index
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
    
    dir_Modulators = strcat(dir_Sess,'/Modulators');
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
            
            % -- remove outliers from modulator, sender, and receiver
            
            
            lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
            outliers_tot = sess_data_lfp.outliers_tot(cnt_m).idx;  % -- get the M,R,S shared outliers
            
            % -- Sender and Receiver LFP
            lfp_S = sess_data_lfp.lfp_S;
            lfp_R = sess_data_lfp.lfp_R;
            % -- remove outliers from sender, receiver, and control
            lfp_S(outliers_tot,:) = [];
            lfp_R(outliers_tot,:) = [];
            lfp_E(outliers_tot,:) = [];
            
             norm_Ts = sum(sum(lfp_S.*lfp_S,2))*0.001;
             norm_Tr = sum(sum(lfp_R.*lfp_R,2))*0.001;
             norm_Te = sum(sum(lfp_E.*lfp_E,2))*0.001;
             
             lfp_S = lfp_S/norm_Ts;
             lfp_R = lfp_R/norm_Tr;
             lfp_E = lfp_R/norm_Te;
             
             S_s = double(dmtspec(lfp_S,[1 3],fs,fk, pad, 0.05, 1));
             S_r = double(dmtspec(lfp_R,[1 3],fs,fk, pad, 0.05, 1));
             S_m = double(dmtspec(lfp_E,[1 3],fs,fk, pad, 0.05, 1));

            
            sess_data_lfp.lfp_E_clean(cnt_m).lfp = lfp_E;   % -- save to structure
            
            % -- coherence for modulator-sender, modulator-receiver
            display(['Computing modulator-sender coherence...'])
            [c_ms,f] = coherency(lfp_E,lfp_S,[N W],fs,fk,pad,0.05,1,1);
            
            
            display(['Computing modulator-receiver coherence...'])
            [c_mr,f] = coherency(lfp_E,lfp_R,[N W],fs,fk,pad,0.05,1,1);
            
                 % -- structure assignements
            mod(cnt_el).c_ms = c_ms ; % assign M-S coherence value for this modulator
            mod(cnt_el).c_mr = c_mr;  % M-R coherence
            mod(cnt_el).s_m = S_m; % Modulator spectrum
            
                 
            % --- FIGURE --------- %%
%             % -- Coherence vs frequency --- %
%             fig = figure;
%             plot(f,abs(c_sr))
%             hold on
%             plot(f,abs(c_ms))
%             hold on
%             plot(f,abs(c_mr))
%             grid on
%             title(sprintf('Abs coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
%             legend('S-R coherence','M-S coherence','M-R coherence')
%             %         xlim([0 60])
%             set(gcf, 'Position',  [100, 600, 1000, 500])
%             
%             fname = strcat(dir_Modulators,sprintf('/coherency_vs_freq_ch_%d_fk_%d.jpg',Ch,fk));
%             saveas(fig,fname);
            
       
            
%             % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %   FIGURES     %%%%%%%%%%%%%%%%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             % --- full length with artifacts
%             lfp_S_rshape = reshape(sess_data_lfp.lfp_S',[],1)';
%             lfp_R_rshape = reshape(sess_data_lfp.lfp_R',[],1)';
%             lfp_E_rshape = reshape(sq(sess_data_lfp.lfp_E(Ch,:,:))',[],1)';
%             
%             fig = figure;
%             plot(lfp_S_rshape)
%             hold on
%             plot(lfp_R_rshape)
%             hold on
%             plot(lfp_E_rshape)
%             grid on
%             title(sprintf('full length - modulator %d',Ch),'FontSize',11)
%             legend('Sender','Receiver','Modulator')
%             set(gcf, 'Position',  [100, 600, 1000, 500])
%             
%             fig_name = strcat(dir_Modulators,sprintf('/LFP_S-R-M_full_length_mod_%d.fig',Ch));
%             saveas(fig,fig_name);
%             fig_name = strcat(dir_Modulators,sprintf('/LFP_S-R-M_full_length_mod_%d.png',Ch));
%             saveas(fig,fig_name);
%             
%               % -- full length without artifacts
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
%             title('Cleaned version LFP ','FontSize',11)
%             legend('Sender','Receiver','Modulator')
%             set(gcf, 'Position',  [100, 600, 1000, 500])
%             
%             fig_name = strcat(dir_Modulators,sprintf('/LFP_S-R-M_cleaned_version_no-artifacts_%d.fig',Ch));
%             saveas(fig,fig_name);
%             fig_name = strcat(dir_Modulators,sprintf('/LFP_S-R-M_cleaned_version_no-artifacts_%d.png',Ch));
%             saveas(fig,fig_name);
%             
            cnt_el = cnt_el + 1; % total modulators counter                     
        end
        cnt_m = cnt_m + 1; % counter for modulators within this session

    end
end


keyboard

dir_SpecNorm = strcat(dir_RS,'/Modulators_SpecNorm');
if ~exist(dir_SpecNorm, 'dir')
    mkdir(dir_SpecNorm)
end

% Save coherence and spectrum data in structure format
save(strcat(dir_SpecNorm,sprintf('/coh_specNorm_m_fk_%d_W_%d.mat',fk,W)),'mod');
save(strcat(dir_SpecNorm,sprintf('/coh_specNorm_sr_fk_%d_W_%d.mat',fk,W)),'stim');

% -- load structure files
 



