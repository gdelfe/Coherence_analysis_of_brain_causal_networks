


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/STIM_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Stim = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';


fk = 200; W = 5;
% %%%%%%%%% RESTING STATE  %%%%%%

% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_m_fk_%d_W_%d.mat',fk,W))); % structure mod
load(strcat(dir_RS,sprintf('/coh_spec_sr_fk_%d_W_%d.mat',fk,W))); % structure stim
stim_mod_RS = stim;
mod_mod_RS = mod;
modulators_RS = mean_coh_and_spec_RS(mod,stim);

%%%%%%%%% CONTROLS SAME AREA %%%%%%%%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_m_all_Controls_same_area_fk_%d_W_%d.mat',fk,W)));
load(strcat(dir_RS,sprintf('/coh_spec_sr_all_Controls_same_area_fk_%d_W_%d.mat',fk,W)));
stim_ctrl_SA_RS = stim;
mod_ctrl_SA_RS = mod;
ctrl_SA_RS = mean_coh_and_spec_RS(mod,stim);

%%%%%%%%% CONTROLS OTHER AREAS %%%%%%%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_m_all_Controls_other_areas_fk_%d_W_%d.mat',fk,W)));
load(strcat(dir_RS,sprintf('/coh_spec_sr_all_Controls_other_areas_fk_%d_W_%d.mat',fk,W)));
stim_ctrl_OA_RS = stim;
mod_ctrl_OA_RS = mod;
ctrl_OA_RS = mean_coh_and_spec_RS(mod,stim);

% %%%%%%%%% STIMuLATION  %%%%%%

% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_Stim,sprintf('/coh_spec_mr_fk_%d_W_%d.mat',fk,W)));
stim_mod_STIM = stim;
modulators_STIM = mean_coh_and_spec_STIM(stim);


%%%%%%%%% CONTROLS SAME AREA %%%%%%%%%%%%
load(strcat(dir_Stim,sprintf('/coh_spec_mr_controls_same_area_fk_%d_W_%d.mat',fk,W)));
stim_ctrl_SA_STIM = stim;
ctrl_SA_STIM = mean_coh_and_spec_STIM(stim);

%%%%%%%%% CONTROLS OTHER AREAS %%%%%%%%%%%%

load(strcat(dir_Stim,sprintf('/coh_spec_mr_controls_other_areas_fk_%d_W_%d.mat',fk,W)));
stim_ctrl_OA_STIM = stim;
ctrl_OA_STIM = mean_coh_and_spec_STIM(stim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,fk,size(modulators_RS.mean_coh_ms,2)); % frequency values (range)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES MR RS vs MR STIM %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all

shadedErrorBar(f,modulators_RS.mean_coh_mr,modulators_RS.err_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators_STIM.mean_coh_mr,modulators_STIM.err_coh_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on

grid on
title('Abs coherence Mod-Rec Resting State vs STIM','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Resting State','STIM','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Stim,sprintf('/coherency_MR_RS_vs_STIM_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)

% --- CONTROL-RECEIVER coherence  RESTING STATE vs STIMULATION  -------%

fig = figure;
hold all

shadedErrorBar(f,ctrl_SA_RS.mean_coh_mr,ctrl_SA_RS.err_mr,'lineprops',{'color',[179, 0, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA_STIM.mean_coh_mr,ctrl_SA_STIM.err_coh_mr,'lineprops',{'color',[255, 92, 51]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA_RS.mean_coh_mr,ctrl_OA_RS.err_mr,'lineprops',{'color',[0, 0, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA_STIM.mean_coh_mr,ctrl_OA_STIM.err_coh_mr,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on

grid on
title('Abs coherence Controls-Receivers - Resting State vs STIM','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('CR same area - Resting State','CR same area - STIM','CR other areas - Resting State','CR other areas - STIM','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Stim,sprintf('/coherency_CR_RS_vs_STIM_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)


% --- MODULATOR-RECEIVER and CONTROL-RECEIVER coherence  RESTING STATE vs STIMULATION  -------%

fig = figure;
hold all

shadedErrorBar(f,modulators_RS.mean_coh_mr,modulators_RS.err_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators_STIM.mean_coh_mr,modulators_STIM.err_coh_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA_RS.mean_coh_mr,ctrl_SA_RS.err_mr,'lineprops',{'color',[179, 0, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA_STIM.mean_coh_mr,ctrl_SA_STIM.err_coh_mr,'lineprops',{'color',[255, 92, 51]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA_RS.mean_coh_mr,ctrl_OA_RS.err_mr,'lineprops',{'color',[0, 0, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA_STIM.mean_coh_mr,ctrl_OA_STIM.err_coh_mr,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on

grid on
title('Abs MR and CR coherence - Resting State vs STIM','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('MR Resting State','MR STIM','CR same area - Resting State','CR same area - STIM','CR other areas - Resting State','CR other areas - STIM','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Stim,sprintf('/coherency_MR_and_CR_RS_vs_STIM_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FIGURES  SPECTRUMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% SPECTRUMS MODULATORS vs RECEIVERS vs CONTROLS %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all


shadedErrorBar(f,modulators_RS.mean_spec_m,modulators_RS.err_S_m,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators_STIM.mean_spec_m,modulators_STIM.err_S_m,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA_RS.mean_spec_m,ctrl_SA_RS.err_S_m,'lineprops',{'color',[179, 0, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA_STIM.mean_spec_m,ctrl_SA_STIM.err_S_m,'lineprops',{'color',[255, 92, 51]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA_RS.mean_spec_m,ctrl_OA_RS.err_S_m,'lineprops',{'color',[0, 0, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA_STIM.mean_spec_m,ctrl_OA_STIM.err_S_m,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on



grid on
title('Spectrum Modulators and Controls - RS vs STIM','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
legend('Modulators RS','Modulators STIM','Ctrl same area - RS','Ctrl same area - STIM','Ctrl other areas - RS','Ctrl other areas - STIM','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Stim,sprintf('/spectrum_Modulators_and_Controls_RS_vs_STMI_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)


figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all
for i=1:18
    
plot(f,abs(stim_mod(i).s_s))
hold on
end
hold on 
shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
grid on


set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all


shadedErrorBar(f,modulators.mean_spec_m,modulators.err_S_m,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators.mean_spec_r,modulators.err_S_r,'lineprops',{'color',[0255, 102, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on

grid on
title('RS: Modulators, Senders, Receivers spectrum','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
legend('Modulators','Receiver','Sender','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_RS,sprintf('/spectrum_RS_modulators_sender_receiver_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)


