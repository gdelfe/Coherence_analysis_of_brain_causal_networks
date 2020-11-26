
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the causal 
% modulators found by Shaoyu's and both the sender and the receiver
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators
%
% In contrast to other codes that employes the coherence-gram to estimate
% the coherence vs frequency, this code employes directly coherency.m
%
% INPUT: file with session modulator info
%        .mat file with structure AM and MA information 
%
% OUTPUT: txt files with the values of the coherence MR, MS, SR and
% corresponding figures
%
%    @ Gino Del Ferraro, August 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_base = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
step = 110;

fid = fopen(strcat(dir_base,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

% -- load structure files
load(strcat(dir_base,'/session_AM.mat'))
load(strcat(dir_base,'/session_MA.mat'))


% -- print structures on stdout
%format short
for s=1:size(sess_info{1},1)
    session_AM(s)
    session_MA(s)
end


cnt_m = 1; % counter for the modulator-receiver/sender coherencies 
cnt_sr = 1; % counter sender-receiver coherencies 

for i=1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    % addpath('/vol/sas8/Maverick_RecStim_vSUBNETS220/160125/004')
    addpath(sprintf('/vol/sas8/Maverick_RecStim_vSUBNETS220/%s/%s/',sess_info{2}{i},sess_info{3}{i})) % add path of the specific RS session
    
    % file = 'rec004.Frontal.lfp.dat'
    file = sprintf('rec%s.Frontal.lfp.dat',sess_info{3}{i})
    fid = fopen(file);
    format = 'float=>single';
    
    CH = 220; % tot number of channels
    FS = 1000; % sampling
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    data = fread(fid,[CH,inf],format); % load the RS data
    % h = fread(fid,[CH,diff(bn)*FS./1e3],format);
    % ---- bipolar referencing, pairs of electrodes
    dir_Sess = strcat(dir_base,sprintf('/Sess_%d',Sess));
    if ~exist(dir_Sess, 'dir')
        mkdir(dir_Sess)
    end
    
    % -- load list electrodes, sender, receiver
    electrode = importdata(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairs
    receiver = importdata(strcat(dir_Sess,sprintf('/receiver_Sess_%d.txt',Sess)));  % ---- receiver pair
    sender = importdata(strcat(dir_Sess,sprintf('/sender_Sess_%d.txt',Sess))); % ---- sender pair
    
    % ---  time parameter
    tot_time = 150001;
    % ---  freq parameter for the masking
    fmin = 10;
    fmax = 40;
    
    % ---- Lfp of the resting state for that specific pair of electrodes
    lfp_E_temp = data(electrode(:,1),:) - data(electrode(:,2),:); % all the electrodes
    lfp_S_temp = data(sender(1),:) - data(sender(2),:); % sender
    lfp_R_temp = data(receiver(1),:) - data(receiver(2),:); % receiver
    
    % include signal up to time where signal is not corrupted
    lfp_E_temp = lfp_E_temp(:,1:tot_time);
    lfp_S_temp = lfp_S_temp(:,1:tot_time);
    lfp_R_temp = lfp_R_temp(:,1:tot_time);
    
    
    % create matrices to store the spitted data: trial x time
    lfp_S = zeros(floor(tot_time/FS),1000);
    lfp_R = zeros(floor(tot_time/FS),1000);
    lfp_E = zeros(size(lfp_E_temp,1),floor(tot_time/FS),1000); % channel x trial x time 
    
   
    % -- sanity check LFP
    figure;
    plot(lfp_S_temp)
    hold on
    plot(lfp_R_temp)
    hold on 
    plot(lfp_E_temp(6,:,:)-lfp_R_temp)
    legend('sender','receiver','electrode')
    
    % split the Lenghty RS time series into 1000 ms windows
    % format: channel x win_indx xtime. For R and S size_channel = 1  
    delta = 1000;
    cnt = 1;
    for j = 0:delta:(tot_time - delta)
        lfp_S(cnt,:) = lfp_S_temp(1,j+1:j+delta);
        lfp_R(cnt,:) = lfp_R_temp(1,j+1:j+delta);
        
        lfp_E(:,cnt,:) =  lfp_E_temp(:,j+1:j+delta,:);
        cnt = cnt + 1;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- COHERENCE- Sender-Receiver -- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ---- parameters for the coherence-gram
    nt = tot_time;
    fs = 1000;
    fk = 200;
    pad = 2;    
    % --- coherence
    
    mod_Ch = session_AM(i).mod_idx; % causal modulator channel
  
    display(['Computing sender-receiver coherence...'])
    % -- coherence calculation via coherency()                
    [c_sr,f,S_s,S_r] = coherency(lfp_S,lfp_R,[1 5],fs,fk,pad,0.05,1,1);
  
    % -- store coherence values sender-receiver and spectrums 
    stim(cnt_sr).c_sr = c_sr; % assign S-R coherence value
    stim(cnt_sr).s_s = S_s; % assign sender spectrum 
    stim(cnt_sr).s_r = S_r; % receiver spectrum 
    cnt_sr = cnt_sr + 1;    % sender/receiver counter 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- COHERENCE- Modulator - Sender/Receiver -- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for Ch = mod_Ch % for all the modulators in the session
        
        close all        
        
        % -- coherence for modulator-sender, modulator-receiver 
        display(['Computing modulator-sender coherence...'])
        [c_ms,f,S_m,S_s] = coherency(sq(lfp_E(Ch,:,:)),lfp_S,[1 5],fs,fk,pad,0.05,1,1);
       
        
        display(['Computing modulator-receiver coherence...'])
        [c_mr,f,S_m,S_r] = coherency(sq(lfp_E(Ch,:,:)),lfp_R,[1 5],fs,fk,pad,0.05,1,1);
        
 
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
    
        fname = strcat(dir_Sess,sprintf('/coherency_vs_freq_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig,fname);
        
        % -- structure assignements 
        mod(cnt_m).c_ms = c_ms ; % assign M-S coherence value for this modulator
        mod(cnt_m).c_mr = c_mr;  % M-R coherence 
        mod(cnt_m).s_m = S_m; % Modulator spectrum
       
        cnt_m = cnt_m + 1; % modulators counter 
        
    end
    
    
end

keyboard 

% Save coherence and spectrum data in structure format
save(strcat(dir_base,sprintf('/coh_spec_m_fk_%d.mat',fk)),'mod');
save(strcat(dir_base,sprintf('/coh_spec_sr_fk_%d.mat',fk)),'stim');

% -- load structure files
fk = 200;
load(strcat(dir_base,sprintf('/coh_spec_m_fk_%d.mat',fk)))
load(strcat(dir_base,sprintf('/coh_spec_sr_fk_%d.mat',fk)))

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
mean_cho_ms = mean(coh_ms);  % modulator - sender
mean_cho_mr = mean(coh_mr);  % modulator - receiver
mean_cho_sr = mean(coh_sr);  % sender - receiver 

% --- std coherences
std_cho_ms = std(coh_ms);  % modulator - sender
std_cho_mr = std(coh_mr); % modulator - receiver
std_cho_sr = std(coh_sr);  % modulator - receiver

% --- Error bars
err_ms = std_cho_ms/sqrt(48);
err_mr = std_cho_mr/sqrt(48);
err_sr = std_cho_sr/sqrt(20);



set(0,'DefaultFigureVisible','on')

% -- FIGURE: Plot average coherence across sessions for MR, SR, MS
fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all

shadedErrorBar(f,mean_cho_ms,err_ms,'lineProps','b','patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[255, 83, 26]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_sr,err_sr,'lineprops',{'color',[230, 184 , 0]/255},'patchSaturation',0.4); hold on


% shadedErrorBar(f,mean_cho_ms,err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_sr,err_sr,'lineprops',{'color',[0 204 204]/255},'patchSaturation',0.4); hold on

grid on
title('Mean Abs coherence of MS, SR, SR','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('M-S mean abs','M-R mean abs','S-R mean abs','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% xlim([0 60])
set(gcf, 'Position',  [100, 600, 1000, 600])


fname = strcat(dir_base,sprintf('/coherency_mean_split-data_MS_MR_SR_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)


