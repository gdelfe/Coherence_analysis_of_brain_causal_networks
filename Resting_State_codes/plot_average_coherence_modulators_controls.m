

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';

f = importdata(strcat(dir_RS,'/frequency_range_fk_200.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           MODULATORS   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fk = 200; W = 5;
load(strcat(dir_RS,sprintf('/coh_spec_m_fk_%d_W_%d.mat',fk,W)));
load(strcat(dir_RS,sprintf('/coh_spec_sr_fk_%d_W_%d.mat',fk,W)));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           CONTROLS   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_m_all_Controls_same_area_fk_%d_W_%d.mat',fk,W)));
load(strcat(dir_RS,sprintf('/coh_spec_sr_all_Controls_same_area_fk_%d_W_%d.mat',fk,W)));

% -- structures to matrices
mod_mat = cell2mat(struct2cell(mod)); % transform struct to mat for modulators
stim_mat = cell2mat(struct2cell(stim)); % transform struct to mat for sender-receiver


% -- assign fields to matrices
coh_ms_ctrl = sq(mod_mat(1,:,:))'; % 1st field, c_ms
coh_mr_ctrl = sq(mod_mat(2,:,:))'; %  2nd field, c_mr
spec_m_ctrl = sq(mod_mat(3,:,:))'; %  3rd field, spec_m

coh_sr_ctrl = sq(stim_mat(1,:,:))'; % 1st field, c_sr
spec_s_ctrl = sq(stim_mat(2,:,:))'; %  2nd field, spec_s
spec_r_ctrl = sq(stim_mat(3,:,:))'; %  3rd field, spec_r

% --- mean coherences
mean_cho_ms_ctrl = mean(abs(coh_ms_ctrl));  % modulator - sender
mean_cho_mr_ctrl = mean(abs(coh_mr_ctrl));  % modulator - receiver
mean_cho_sr_ctrl = mean(abs(coh_sr_ctrl));  % sender - receiver

% --- std coherences
std_cho_ms_ctrl = std(abs(coh_ms_ctrl));  % modulator - sender
std_cho_mr_ctrl = std(abs(coh_mr_ctrl)); % modulator - receiver
std_cho_sr_ctrl = std(abs(coh_sr_ctrl));  % modulator - receiver

% --- Error bars
M = size(coh_ms_ctrl,1);
S = size(coh_sr_ctrl,1);
err_ms_ctrl = std_cho_ms_ctrl/sqrt(M);
err_mr_ctrl = std_cho_mr_ctrl/sqrt(M);
err_sr_ctrl = std_cho_sr_ctrl/sqrt(S);

keyboard




set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, SR, MS
fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all


shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_mr_ctrl,err_mr_ctrl,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on

grid on
title('Abs coherency CONTROLS of MR, CR-control - Resting State','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('M-R abs coherency','C-R abs coherency','S-R abs coherency','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])



fname = strcat(dir_RS,sprintf('/coherency_mean_all_Controls_same_area_split-data_MS_MR_SR_W_%d_fk_%d-all-Sess.png',W,fk));
saveas(fig,fname)
fname = strcat(dir_RS,sprintf('/coherency_mean_all_Controls_same_area_split-data_MS_MR_SR_W_%d_fk_%d-all-Sess.fig',W,fk));
saveas(fig,fname)