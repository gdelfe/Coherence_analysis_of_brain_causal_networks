

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/vol/bd5/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/vol/bd5/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
% dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
dir_Maverick = strcat(dir_main,sprintf('Maverick/Resting_state/%s',freq_band));
dir_Archie = strcat(dir_main,sprintf('Archie/Resting_state/%s',freq_band));

recording_mav = 'last_recording';
recording_arc = 'rec001_002_all_sessions';

fn_fig = '_lastRec_rec001-002_all';
fn_arc = '_rec001_002_all_sess'

dir_Ctrl_Maverick = strcat(dir_Maverick,sprintf('/Modulators_Controls_avg_results/%s',recording_mav));
dir_Ctrl_Archie = strcat(dir_Archie,sprintf('/Modulators_Controls_avg_results/%s',recording_arc));


fk = 200; W = 5;
% %%%%%%%%% MODULATORS Maverick %%%%%%
load(strcat(dir_Ctrl_Maverick,sprintf('/coh_spec_m_fk_%d_W_%d.mat',fk,W))); % structure mod
load(strcat(dir_Ctrl_Maverick,sprintf('/coh_spec_sr_fk_%d_W_%d.mat',fk,W))); % structure stim
struct_mod_mav = mod;
struct_stim_mav = stim;

% %%%%%%%%% MODULATORS Archie %%%%%%
load(strcat(dir_Ctrl_Archie,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,fn_arc))); % structure mod
load(strcat(dir_Ctrl_Archie,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,fn_arc))); % structure stim
struct_mod_arc = mod;
struct_stim_arc = stim;

modulators = mean_coh_and_spec_across_monkeys(struct_mod_mav,struct_stim_mav,struct_mod_arc,struct_stim_arc);



%%%%%%%%% CONTROLS SAME AREA Maverick %%%%%%%%%%%%
load(strcat(dir_Ctrl_Maverick,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d.mat',fk,W)));
load(strcat(dir_Ctrl_Maverick,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d.mat',fk,W)));
mod_ctrl_SA_mav = mod;
stim_ctrl_SA_mav = stim;

%%%%%%%%% CONTROLS SAME AREA Archie %%%%%%%%%%%%
load(strcat(dir_Ctrl_Archie,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s.mat',fk,W,fn_arc)));
load(strcat(dir_Ctrl_Archie,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s.mat',fk,W,fn_arc)));
mod_ctrl_SA_arc = mod;
stim_ctrl_SA_arc = stim;


ctrl_SA = mean_coh_and_spec_across_monkeys(mod_ctrl_SA_mav,stim_ctrl_SA_mav,mod_ctrl_SA_arc,stim_ctrl_SA_arc);


%%%%%%%%% CONTROLS OTHER AREA Maverick %%%%%%%%%%%%
load(strcat(dir_Ctrl_Maverick,sprintf('/coh_spec_m_Controls_other_areas_fk_%d_W_%d.mat',fk,W)));
load(strcat(dir_Ctrl_Maverick,sprintf('/coh_spec_sr_Controls_other_areas_fk_%d_W_%d.mat',fk,W)));
mod_ctrl_OA_mav = mod;
stim_ctrl_OA_mav = stim;

%%%%%%%%% CONTROLS OTHER AREA Archie %%%%%%%%%%%%
load(strcat(dir_Ctrl_Archie,sprintf('/coh_spec_m_Controls_other_areas_fk_%d_W_%d%s.mat',fk,W,fn_arc)));
load(strcat(dir_Ctrl_Archie,sprintf('/coh_spec_sr_Controls_other_areas_fk_%d_W_%d%s.mat',fk,W,fn_arc)));
mod_ctrl_OA_arc = mod;
stim_ctrl_OA_arc = stim;


ctrl_OA = mean_coh_and_spec_across_monkeys(mod_ctrl_OA_mav,stim_ctrl_OA_mav,mod_ctrl_OA_arc,stim_ctrl_OA_arc);



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
title('Both animals: Abs MR coherence MODULATORS vs CONTROLS, rec001 - Resting State','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Receivers','Controls-Receivers  same area','Controls-Receiver  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on
xlim([0 90])
xlim([0 95])

fname = strcat(dir_Ctrl_Maverick,sprintf('/coherency_MR_Modulators_vs_Controls_both_monkeys_W_%d_fk_%d%s.png',W,fk,fn_fig));
saveas(fig,fname)
fname = strcat(dir_Ctrl_Maverick,sprintf('/coherency_MR_Modulators_vs_Controls_both_monkeys_W_%d_fk_%d%s.fig',W,fk,fn_fig));
saveas(fig,fname)

% --- ELECTRODE-SENDER coherence   -------%

fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_ms,ctrl_SA.err_ms,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_ms,ctrl_OA.err_ms,'lineprops',{'color',[255, 128, 128]/255},'patchSaturation',0.4); hold on

grid on
title('Both animals: Abs MS coherence MODULATORS vs CONTROLS, rec001 - Resting State','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Senders','Controls-Senders  same area','Controls-Senders  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Ctrl_Maverick,sprintf('/coherency_MS_Modulators_vs_Controls_both_monkeys_W_%d_fk_%d%s.png',W,fk,fn_fig));
saveas(fig,fname)
fname = strcat(dir_Ctrl_Maverick,sprintf('/coherency_MS_Modulators_vs_Controls_both_monkeys_W_%d_fk_%d%s.fig',W,fk,fn_fig));
saveas(fig,fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FIGURES  SPECTRUMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% SPECTRUMS MODULATORS vs RECEIVERS vs CONTROLS %%%%%%%%%%%%

 
% set(0,'DefaultFigureVisible','on')
% % -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
% fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
% hold all
% 
% 
% shadedErrorBar(f,modulators.mean_spec_m,modulators.err_S_m,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,modulators.mean_spec_r,modulators.err_S_r,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); hold on
% 
% grid on
% title('RS: Modulators, Senders, Receivers spectrum','FontSize',11);
% xlabel('freq (Hz)');
% ylabel('spectrum');
% legend('Modulators','Receiver','Sender','FontSize',10)
% set(gcf, 'Position',  [100, 600, 1000, 600])
% grid on
% 
% fname = strcat(dir_Controls,sprintf('/spectrum_modulators_sender_receiver_W_%d_fk_%d.png',W,fk));
% saveas(fig,fname)
% fname = strcat(dir_Controls,sprintf('/spectrum_modulators_sender_receiver_W_%d_fk_%d.fig',W,fk));
% saveas(fig,fname)
% 
% 
% figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
% hold all
% for i=1:18
%     
% plot(f,abs(stim_mod(i).s_s))
% hold on
% end
% hold on 
% shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% grid on
% 
% 
% set(0,'DefaultFigureVisible','on')
% % -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
% fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
% hold all
% 
% 
% shadedErrorBar(f,modulators.mean_spec_m,modulators.err_S_m,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,modulators.mean_spec_r,modulators.err_S_r,'lineprops',{'color',[0255, 102, 0]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% 
% grid on
% title('RS: Modulators, Senders, Receivers spectrum','FontSize',11);
% xlabel('freq (Hz)');
% ylabel('spectrum');
% legend('Modulators','Receiver','Sender','FontSize',10)
% set(gcf, 'Position',  [100, 600, 1000, 600])
% grid on
% 
% fname = strcat(dir_RS,sprintf('/spectrum_RS_modulators_sender_receiver_W_%d_fk_%d.png',W,fk));
% saveas(fig,fname)
% 
