

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the average coherence for the Receiver - Sender_controls
% vs the Receiver - Sender average coherence

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)



addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Sender = strcat(dir_RS,'/Sender_controls');

fk = 200; W = 5;
% %%%%%%%%% Sender - Receiver coherence data  %%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_sr_fk_%d_W_%d.mat',fk,W)),'stim'); % -- stim structure  
data_SR = mean_SR_coh_and_spec_RestState(stim);
clear stim

% %%%%%%%%% Sender ctrl same area - Receivers %%%%%%
load(strcat(dir_Sender,sprintf('/coh_spec_SR_Send_ctrl_same_area_fk_%d_W_%d.mat',fk,W))); % structure stim_ctrl
ctrl_SA = mean_SR_coh_and_spec_RestState(stim_ctrl);
clear stim

% %%%%%%%%% Sender ctrl other areas - Receivers %%%%%%
load(strcat(dir_Sender,sprintf('/coh_spec_SR_Send_ctrl_other_areas_fk_%d_W_%d.mat',fk,W))); % structure stim_ctrl
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
title('Abs SR coherence - Sender controls vs Receiver - Resting State','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Sender-Receivers','Sender same area - Receiver','Sender other areas - Receiver','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Sender,sprintf('/coherency_SR_Sender_ctrl_vs_Receiver_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)
fname = strcat(dir_Sender,sprintf('/coherency_SR_Sender_ctrl_vs_Receiver_W_%d_fk_%d.fig',W,fk));
saveas(fig,fname)

