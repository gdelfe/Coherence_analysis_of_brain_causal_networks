

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)



addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Receiver = strcat(dir_RS,'/Receiver_controls');
dir_Sender = strcat(dir_RS,'/Sender_controls');

fk = 200; W = 5;
% %%%%%%%%% Sender - Receivers ctrl same area  %%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_sr_fk_%d_W_%d.mat',fk,W)),'stim'); % -- stim structure  
data_SR = mean_SR_coh_and_spec_RestState(stim);

% %%%%%%%%% Sender - Receivers ctrl same area  %%%%%%
load(strcat(dir_Receiver,sprintf('/coh_spec_SR_Rec_ctrl_same_area_fk_%d_W_%d.mat',fk,W))); % structure mod
ctrl_SA = mean_SR_coh_and_spec_RestState(stim_ctrl);


% %%%%%%%%% Sender - Receivers ctrl other areas  %%%%%%
load(strcat(dir_Receiver,sprintf('/coh_spec_SR_Rec_ctrl_other_areas_fk_%d_W_%d.mat',fk,W))); % structure mod
ctrl_OA = mean_SR_coh_and_spec_RestState(stim_ctrl);


f = linspace(1,fk,size(data_SR.mean_coh_sr,2)); % frequency values (range)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES Sender vs Receiver controls  %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all


shadedErrorBar(f,data_SR.mean_coh_sr,data_SR.err_sr,'lineprops',{'color',[0, 51, 153]/255},'patchSaturation',0.5); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_sr,ctrl_SA.err_sr,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.5); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_sr,ctrl_OA.err_sr,'lineprops',{'color',[0, 255, 255]/255},'patchSaturation',0.5); hold on

grid on
title('Abs SR coherence - Sender vs Receiver Controls - Resting State','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Sender-Receivers','Sender-Receivers same area','Sender-Receiver  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Receiver,sprintf('/coherency_SR_Sender_vs_Receiver_ctrl_W_%d_fk_%d-all-Sess.png',W,fk));
saveas(fig,fname)

