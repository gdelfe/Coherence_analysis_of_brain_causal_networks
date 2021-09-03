

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

filename = ''; % -- filename for sess_data_info.mat 
figname = '';
recording = 'last_recording';
title_caption = 'last rec';

freq_band = 'theta_band';
monkey = 'Maverick';

dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_avg = strcat(dir_RS_Theta,sprintf('/Modulators_Controls_avg_results/%s',recording));



fk = 200; W = 5;
% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_avg,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,filename))); % structure mod
load(strcat(dir_avg,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,filename))); % structure stim
stim_mod = stim;
mod_mod = mod;
modulators = mean_coh_and_spec_RS(mod,stim);

%%%%%%%%% CONTROLS SAME AREA %%%%%%%%%%%%
load(strcat(dir_avg,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s.mat',fk,W,filename)));
load(strcat(dir_avg,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s.mat',fk,W,filename)));
stim_ctrl_SA = stim;
mod_ctrl_SA = mod;
ctrl_SA = mean_coh_and_spec_RS(mod,stim);


%%%%%%%%% CONTROLS OTHER AREAS %%%%%%%%%%%
load(strcat(dir_avg,sprintf('/coh_spec_m_Controls_other_areas_fk_%d_W_%d%s.mat',fk,W,filename)));
load(strcat(dir_avg,sprintf('/coh_spec_sr_Controls_other_areas_fk_%d_W_%d%s.mat',fk,W,filename)));
stim_ctrl_OA = stim;
mod_ctrl_OA = mod;
ctrl_OA = mean_coh_and_spec_RS(mod,stim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,fk,size(modulators.mean_coh_ms,2)); % frequency values (range)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES MODULATORS vs CONTROLS %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_mr,ctrl_SA.err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_mr,ctrl_OA.err_mr,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('Mean Abs MR coherence THETA MODULATORS vs CONTROLS - %s',title_caption),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Receivers','Controls-Receivers same area','Controls-Receiver other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_avg,sprintf('/coherency_MR_Theta_Modulators_vs_Controls_W_%d_fk_%d%s.png',W,fk,figname));
saveas(fig,fname)
fname = strcat(dir_avg,sprintf('/coherency_MR_Theta_Modulators_vs_Controls_W_%d_fk_%d%s.fig',W,fk,figname));
saveas(fig,fname)

% --- ELECTRODE-SENDER coherence   -------%

fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_ms,ctrl_SA.err_ms,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_ms,ctrl_OA.err_ms,'lineprops',{'color',[255, 128, 128]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('Mean Abs MS coherence Theta MODULATORS vs CONTROLS - %s',title_caption),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Senders','Controls-Senders  same area','Controls-Senders  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_avg,sprintf('/coherency_MS_Theta_Modulators_vs_Controls_W_%d_fk_%d%s.png',W,fk,figname));
saveas(fig,fname)
fname = strcat(dir_avg,sprintf('/coherency_MS_Theta_Modulators_vs_Controls_W_%d_fk_%d%s.fig',W,fk,figname));
saveas(fig,fname)


keyboard 

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

fname = strcat(dir_Controls,sprintf('/spectrum_theta_modulators_sender_receiver_W_%d_fk_%d.png',W,fk));
saveas(fig,fname)
fname = strcat(dir_Controls,sprintf('/spectrum_theta_modulators_sender_receiver_W_%d_fk_%d.fig',W,fk));
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

