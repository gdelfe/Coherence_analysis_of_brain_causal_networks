

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_Stim = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';



fk = 200; W = 5;
% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_Stim,sprintf('/coh_spec_mr_fk_%d_W_%d.mat',fk,W)));
mod = mean_coh_and_spec_STIM(stim);

%%%%%%%%% CONTROLS SAME AREA %%%%%%%%%%%%
load(strcat(dir_Stim,sprintf('/coh_spec_mr_controls_same_area_fk_%d_W_%d.mat',fk,W)));
ctrl_SA = mean_coh_and_spec_STIM(stim);

%%%%%%%%% CONTROLS OTHER AREAS %%%%%%%%%%%%
load(strcat(dir_Stim,sprintf('/coh_spec_mr_controls_other_areas_fk_%d_W_%d.mat',fk,W)));
ctrl_OA = mean_coh_and_spec_STIM(stim);


f = linspace(1,fk,size(mod.mean_coh_mr,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% COHERENCES MODULATORS vs CONTROLS %%%%%%%%%%%%

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all

shadedErrorBar(f,mod.mean_coh_mr,mod.err_coh_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_mr,ctrl_SA.err_coh_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_mr,ctrl_OA.err_coh_mr,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); hold on

grid on
title('STIM: Abs MR coherence MODULATORS vs CONTROLS ','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Receivers','Controls-Receivers  same area','Controls-Receiver  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Stim,sprintf('/coherency_MR_Modulators_vs_Controls_STIM_W_%d_fk_%d-all-Sess.png',W,fk));
saveas(fig,fname)

%%%%%%%%% COHERENCES HITS vs MISSES %%%%%%%%%%%%

% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all

shadedErrorBar(f,mod.mean_coh_mr_H,mod.err_coh_mr_H,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mod.mean_coh_mr_M,mod.err_coh_mr_M,'lineprops',{'color',[255, 102, 0]/255},'patchSaturation',0.4); hold on

grid on
title('STIM: Abs MR coherence MODULATORS: HITS vs MISSES ','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Mod-Receivers HITS','Mod-Receivers MISSES','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Stim,sprintf('/coherency_MR_Modulators_HITS_vs_MISSES_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)


%%%%%%%%% SPECTRUMS MODULATORS vs RECEIVERS %%%%%%%%%%%%

fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all
shadedErrorBar(f,mod.mean_spec_m,mod.err_S_m,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mod.mean_spec_r,mod.err_S_r,'lineprops',{'color',[255, 102, 0]/255},'patchSaturation',0.4); hold on

grid on
title('STIM: Modulators vs Receivers Spectrum ','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('Modulators Spectrum','Receivers Spectrum','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 800, 500])

fname = strcat(dir_Stim,'/mean_spectrum_modulators_vs_receivers.png');
saveas(fig,fname)


%%%%%%%%% SPECTRUMS MODULATORS vs RECEIVERS %%%%%%%%%%%%

fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all
shadedErrorBar(f,mod.mean_spec_m_H,mod.err_S_m_H,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mod.mean_spec_m_M,mod.err_S_m_M,'lineprops',{'color',[0, 0, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mod.mean_spec_r_H,mod.err_S_r_H,'lineprops',{'color',[255, 102, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mod.mean_spec_r_M,mod.err_S_r_M,'lineprops',{'color',[204, 0, 0]/255},'patchSaturation',0.4); hold on

grid on
title('STIM: Modulators vs Receivers Spectrum HITS and MISSES ','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('Mod Spec HITS','Mod Spec MISSES','Rec Spec HITS','Rec Spec MISSES','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 800, 500])

fname = strcat(dir_Stim,'/mean_spectrum_modulators_vs_receivers_Hits_Misses.png');
saveas(fig,fname)


%%%%%%%%% SPECTRUMS MODULATORS vs CONTROLS vs RECEIVERS %%%%%%%%%%%%

fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all
shadedErrorBar(f,mod.mean_spec_m,mod.err_S_m,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_spec_r,ctrl_SA.err_S_r,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_spec_m,ctrl_OA.err_S_m,'lineprops',{'color',[0, 204, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mod.mean_spec_r,mod.err_S_r,'lineprops',{'color',[255, 102, 0]/255},'patchSaturation',0.4); hold on

grid on
title('STIM: Modulators vs Controls vs Receivers Spectrum ','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('Modulators','CTLR same area','CTLR other areas','Receiver','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 800, 500])

fname = strcat(dir_Stim,'/mean_spectrum_modulators_vs_controls_vs_receivers.png');
saveas(fig,fname)


%%%%%%%%% SPECTRUMS MODULATORS vs CONTROLS %%%%%%%%%%%%

fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all
shadedErrorBar(f,mod.mean_spec_m,mod.err_S_m,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_spec_r,ctrl_SA.err_S_r,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_spec_m,ctrl_OA.err_S_m,'lineprops',{'color',[0, 204, 153]/255},'patchSaturation',0.4); hold on

grid on
title('STIM: Modulators vs Controls vs Receivers Spectrum ','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('Modulators','CTRL same area','CTRL other areas','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 800, 500])

fname = strcat(dir_Stim,'/mean_spectrum_modulators_vs_controls.png');
saveas(fig,fname)



