

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';



fk = 200; W = 5;
% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_m_fk_%d_W_%d.mat',fk,W))); % structure mod
load(strcat(dir_RS,sprintf('/coh_spec_sr_fk_%d_W_%d.mat',fk,W))); % structure stim
mod = mean_coh_and_spec_RS(mod,stim);

%%%%%%%%% CONTROLS SAME AREA %%%%%%%%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_m_all_Controls_same_area_fk_%d_W_%d.mat',fk,W)));
load(strcat(dir_RS,sprintf('/coh_spec_sr_all_Controls_same_area_fk_%d_W_%d.mat',fk,W)));
ctrl_SA = mean_coh_and_spec_RS(mod,stim);


%%%%%%%%% CONTROLS OTHER AREAS %%%%%%%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_m_all_Controls_other_areas_fk_%d_W_%d.mat',fk,W)));
load(strcat(dir_RS,sprintf('/coh_spec_sr_all_Controls_other_areas_fk_%d_W_%d.mat',fk,W)));
ctrl_OA = mean_coh_and_spec_RS(mod,stim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,fk,size(coh_mr,2)); % frequency values (range)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES MODULATORS vs CONTROLS %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
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

% --- ELECTRODE-SENDER coherence   -------%

fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all


shadedErrorBar(f,mean_cho_ms,err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_ms_ctrl,err_ms_ctrl,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_ms_ctrl_OA,err_ms_ctrl_OA,'lineprops',{'color',[255, 128, 128]/255},'patchSaturation',0.4); hold on

grid on
title('Abs MS coherence MODULATORS vs CONTROLS - Resting State','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Senders','Controls-Senders  same area','Controls-Senders  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_RS,sprintf('/coherency_MS_Modulators_vs_all_Controls_both_versions_W_%d_fk_%d-all-Sess.png',W,fk));
saveas(fig,fname)









