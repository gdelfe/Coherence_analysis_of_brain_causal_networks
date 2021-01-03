

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_Stim = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           MODULATORS   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fk = 200; W = 5;
load(strcat(dir_Stim,sprintf('/coh_spec_mr_fk_%d_W_%d.mat',fk,W)));

stim_mat = cell2mat(struct2cell(stim)); % transform struct to mat for sender-receiver

coh_mr = sq(stim_mat(1,:,:))'; % 1st field, c_mr
spec_m = sq(stim_mat(2,:,:))'; %  2nd field, spec_m
spec_r = sq(stim_mat(3,:,:))'; %  3rd field, spec_r

coh_mr_H = sq(stim_mat(4,:,:))'; % 1st field, c_mr Hits
spec_m_H = sq(stim_mat(5,:,:))'; %  2nd field, spec_m Hits
spec_r_H = sq(stim_mat(6,:,:))'; %  3rd field, spec_r Hits

coh_mr_M = sq(stim_mat(7,:,:))'; % 1st field, c_mr Hits
spec_m_M = sq(stim_mat(8,:,:))'; %  2nd field, spec_m Hits
spec_r_M = sq(stim_mat(9,:,:))'; %  3rd field, spec_r Hits


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS - MEAN COHERENCES vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- mean coherences
mean_cho = mean(abs(coh_mr));  % modulator - sender
mean_cho_H = mean(abs(coh_mr_H));  % modulator - receiver
mean_cho_M = mean(abs(coh_mr_M));  % sender - receiver

% --- std coherences
std_mr = std(abs(coh_mr));  % modulator - sender
std_mr_H = std(abs(coh_mr_H)); % modulator - receiver
std_mr_M = std(abs(coh_mr_M));  % modulator - receiver

% --- Error bars choerencies
M = size(stim,2); % -- number of electrodes
err_mr = std_mr/sqrt(M);
err_mr_H = std_mr_H/sqrt(M);
err_mr_M = std_mr_M/sqrt(M);

% --- mean spectrum
mean_spec_m = mean(abs(spec_m));  % mean spectrum modulator 
mean_spec_r = mean(abs(spec_r));  % mean spectrum receiver
mean_spec_m_H = mean(abs(spec_m_H));  % mean spectrum modulator Hits
mean_spec_r_H = mean(abs(spec_r_H));  % mean spectrum receiver Hits
mean_spec_m_M = mean(abs(spec_m_M));  % mean spectrum modulator Misses
mean_spec_r_M = mean(abs(spec_r_M));  % mean spectrum receiver Misses

% --- std coherences
std_spec_m = std(abs(spec_m));  % mean spectrum modulator 
std_spec_r = std(abs(spec_r));  % mean spectrum receiver
std_spec_m_H = std(abs(spec_m_H));  % mean spectrum modulator Hits
std_spec_r_H = std(abs(spec_r_H));  % mean spectrum receiver Hits
std_spec_m_M = std(abs(spec_m_M));  % mean spectrum modulator Misses
std_spec_r_M = std(abs(spec_r_M));  % mean spectrum receiver Misses

% --- Error bars choerencies
M = size(stim,2); % -- number of electrodes
err_S_m = std_spec_m/sqrt(M); % -- all
err_S_r = std_spec_r/sqrt(M);
err_S_m_H = std_spec_m_H/sqrt(M); % -- Hits
err_S_r_H = std_spec_r_H/sqrt(M);
err_S_m_M = std_spec_m_M/sqrt(M); % -- misses
err_S_r_M = std_spec_r_M/sqrt(M);

set(0,'DefaultLineLineWidth',2)



set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, SR, MS
fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all


shadedErrorBar(f,mean_cho_mr,err_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_mr_ctrl,err_mr_ctrl,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_mr_ctrl_OA,err_mr_ctrl_OA,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); hold on

grid on
title('Abs MR coherence MODULATORS vs CONTROLS - Resting State','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Receivers','Controls-Receivers  same area','Controls-Receiver  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_RS,sprintf('/coherency_MR_Modulators_vs_all_Controls_both_versions_W_%d_fk_%d-all-Sess.png',W,fk));
saveas(fig,fname)




fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all
errorbar(f,mean_Sm_H,err_Sm_H); hold on
errorbar(f,mean_Sm_M,err_Sm_M); hold on
% errorbar(f,mean_cho_ms_AM,err_ms_AM,'Color',[0.4940, 0.1840, 0.5560]); hold on
% errorbar(f,mean_cho_mr_AM,err_mr_AM,'Color',[102/255, 153/255 0]); hold on
% errorbar(f,mean_cho_sr_AM,err_sr_AM,'color',[26/255 198/255 1]);
grid on
title('STIM: Hits/Misses Spectrum caus mod','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('HITS Mod Spectrum','MISSES Mod Spectrum','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 700, 500])

fname = strcat(dir_Stim,'/mean_spectrum_causal_modulators.png');
saveas(fig,fname)
