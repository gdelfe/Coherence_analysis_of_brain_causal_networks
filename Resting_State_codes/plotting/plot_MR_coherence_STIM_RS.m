

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Stim = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';
dir_coh_RS = strcat(dir_RS,'/Coherence_STIM_RS_same_trials');
fk = 200; W = 5;

% %%%%%%%%% MODULATORS  STIM %%%%%%
load(strcat(dir_Stim,sprintf('/coh_spec_mr_sameRStrails_fk_%d_W_%d.mat',fk,W))); % name structure: stim 
% %%%%%%%%% MODULATORS  RS %%%%%%
load(strcat(dir_coh_RS,sprintf('/coh_spec_MR_only_for_STIM_comparison_fk_%d_W_%d.mat',fk,W))); % name structure: resting

% -- compute coh and spec for RS and STIM with error bars
data = mean_coh_and_spec_RS_STIM(resting,stim);

f = linspace(1,fk,size(data.RS_mean_coh_mr,2)); % frequency values (range)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES STIM vs RESTING %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all


shadedErrorBar(f,data.RS_mean_coh_mr,data.RS_err_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,data.STIM_mean_coh_mr,data.STIM_err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on

title('Abs MR coherence RS vs STIM -- same # of trials','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('RS coherence MR','STIM coherence MR','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_coh_RS,sprintf('/coherency_MR_RS_vs_STIM_W_%d_fk_%d-all-Sess.png',W,fk));
saveas(fig,fname)
fname = strcat(dir_coh_RS,sprintf('/coherency_MR_RS_vs_STIM_W_%d_fk_%d-all-Sess.fig',W,fk));
saveas(fig,fname)



keyboard
%

%
% SPECTRUMS NEEDS TO BE RECALCULATED 
%
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


shadedErrorBar(f,modulators.mean_spec_m,modulators.err_S_m,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators.mean_spec_r,modulators.err_S_r,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); hold on

grid on
title('RS: Modulators, Senders, Receivers spectrum','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
legend('Modulators','Receiver','Sender','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_RS,sprintf('/spectrum_RS_modulators_sender_receiver_W_%d_fk_%d.png',W,fk));
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

