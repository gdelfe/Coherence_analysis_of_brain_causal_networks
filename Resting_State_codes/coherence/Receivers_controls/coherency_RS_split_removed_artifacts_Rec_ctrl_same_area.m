
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the CONTROL
% modulators found by Shaoyu's and both the sender and the receiver
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
name_struct_input = '/sess_Rec_ctrl_other_areas_lfp.mat';


cnt_sr = 1; % counter sender-receiver_ctrl coherencies
list_sess = 1:19;
list_sess(17) = []; % -- Session 17 and 20 are full of artifacts

for i = list_sess %1:size(sess_info{1},1)-1  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Receiver = strcat(dir_RS,sprintf('/Sess_%d/Receiver_controls_other_areas',Sess));
    
    load(strcat(dir_Receiver,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    % ---  time parameter
    tot_time = 150001;
      % ---- parameters for the coherence-gram
    nt = tot_time;
    fs = 1000;
    fk = 200;
    pad = 2;
    N = 1;
    W = 5;
    % --- coherence

    ctrl_Ch = sess_control_lfp.ctrl_idx % control channel indexes for the receiver controls 
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(ctrl_Ch),'  -- contols mod Ch: ',num2str(ctrl_Ch)])
    
    
    % %%%%%%% ALL Electrodes LFP (Receiver controls) %%%%%%%%%%%%%%%%%%%%%
    lfp_R_all = sess_control_lfp.lfp_E;
    
    cnt_m = 1;
    for Ch = ctrl_Ch % for all the receiver controls
        
        close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- COHERENCE - Sender/ Control Receiver -- %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if Ch ~= sess_control_lfp.receiver_idx % if the electrode is not the receiver itself
            
            
            % -- LFPs
            lfp_R = sq(lfp_R_all(Ch,:,:));          % -- get lfp for only that channel -- control receiver
            lfp_S = sess_control_lfp.lfp_S;     % -- sender LFP 
            
            % -- outliers 
            outliers_R = sess_control_lfp.outliers_E(cnt_m).idx;  % -- control receiver outliers 
            outliers_S = sess_control_lfp.outliers_S;  % -- sender outliers
            outliers_SR = [outliers_R, outliers_S]; % 
            
            % -- remove outliers from LFPs
            lfp_R(outliers_SR,:) = [];
            lfp_S(outliers_SR,:) = [];
           
            
            sess_control_lfp.lfp_R_clean(cnt_m).lfp = lfp_R;   % -- save to structure
            
            
            display(['Computing Sender - Control Receiver coherence...'])
            % -- coherence calculation via coherency()
            [c_sr,f,S_s,S_r] = coherency(lfp_S,lfp_R,[N W],fs,fk,pad,0.05,1,1);
                       
            
            % -- store coherence values sender-receiver and spectrums
            stim_ctrl(cnt_sr).c_sr = c_sr; % assign S-R coherence value
            stim_ctrl(cnt_sr).s_s = S_s; % assign sender spectrum
            stim_ctrl(cnt_sr).s_r = S_r; % receiver spectrum
            cnt_sr = cnt_sr + 1;    % sender/receiver counter
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ABS COHERENCE                 %%%%%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % --- FIGURE --------- %%
            % -- Coherence vs frequency --- %
            fig = figure;
            plot(f,abs(c_sr))
            grid on
            title(sprintf('Abs coherence Sender - Receiver-ctrl, receiver ch = %d,',Ch),'FontSize',10);
            legend('S-R ctrl coherence')
            %         xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Receiver,sprintf('/coherency_SR_Rec_ctrl_other_areas_ch_%d_fk_%d.jpg',Ch,fk));
            saveas(fig,fname);
            
        end
        cnt_m = cnt_m + 1; % counter for control within this session

    end    
    
end


dir_R_ctrl = strcat(dir_RS,'/Receiver_controls');
if ~exist(dir_R_ctrl, 'dir')
    mkdir(dir_R_ctrl)
end
% Save coherence and spectrum data in structure format
save(strcat(dir_R_ctrl,sprintf('/coh_spec_SR_Rec_ctrl_other_areas_fk_%d_W_%d.mat',fk,W)),'stim_ctrl');







