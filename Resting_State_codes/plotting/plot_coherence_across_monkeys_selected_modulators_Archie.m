

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';


freq_band = 'beta_band';

N = 30 % --- max number of modulators 
titleN = sprintf('%d top modulators - movie',N);
namefile = '_last_rec_rec001_002.mat';
namefig = '_last_rec_rec001_002.fig';
namepng = '_last_rec_rec001_002.png';

recording_save = 'last_rec-rec001_002';
recording_mav = 'last_recording';
recording_arc = 'rec001_002_no_bad_sessions';
score = 'AUC';
N = 20;
mod_fname_mav = 'AUC';
mod_fname_arc = 'AUC_all';


dir_Maverick = strcat(dir_main,sprintf('Maverick/Resting_state/%s',freq_band));
dir_Archie = strcat(dir_main,sprintf('Archie/Resting_state/%s',freq_band));

dir_ctrl_Maverick = strcat(dir_main,sprintf('Maverick/Resting_state/%s/Modulators_controls',freq_band));
dir_ctrl_Archie = strcat(dir_main,sprintf('Archie/Resting_state/%s/Modulators_controls',freq_band));

dir_avg_coh_Maverick = strcat(dir_Maverick,sprintf('/Modulators_Controls_avg_results/%s/%s',recording_mav,score));
dir_avg_coh_Archie = strcat(dir_Archie,sprintf('/Modulators_Controls_avg_results/%s/%s',recording_arc,score));


fk = 200; W = 5;

dir_both_monkeys = strcat(dir_main,sprintf('both_monkeys/%s/modulators_vs_controls/%s/%s',freq_band,recording_save,score));



% %%%%%%%%% MODULATORS  %%%%%%
mav = load(strcat(dir_avg_coh_Maverick,sprintf('/modulators_N_%d_%s.mat',N,mod_fname_mav))); % structure mod
arc = load(strcat(dir_avg_coh_Archie,sprintf('/modulators_N_%d_%s.mat',N,mod_fname_arc))); % structure mod


% %%%%%%%%% Computes modulators for Maverick and Archie %%%%%%
modulators = mean_coh_across_monkeys(mav,arc);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS SAME AREA               %

mav_SA = load(strcat(dir_avg_coh_Maverick,sprintf('/controls_same_area_N_%d_%s.mat',N,mod_fname_mav))); % structure mod
arc_SA = load(strcat(dir_avg_coh_Archie,sprintf('/controls_same_area_N_%d_%s.mat',N,mod_fname_arc))); % structure mod

ctrl_SA = mean_coh_across_monkeys_SA(mav_SA,arc_SA);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS OTHER AREAS               %

mav_OA = load(strcat(dir_avg_coh_Maverick,sprintf('/controls_other_areas_N_%d_%s.mat',N,mod_fname_mav))); % structure mod
arc_OA = load(strcat(dir_avg_coh_Archie,sprintf('/controls_other_areas_N_%d_%s.mat',N,mod_fname_arc))); % structure mod

ctrl_OA = mean_coh_across_monkeys_OA(mav_OA,arc_OA);


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
title(sprintf('Both animals: Abs MR coherence, %s - Resting State',titleN),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Receivers','Controls-Receivers  same area','Controls-Receiver  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_both_monkeys,sprintf('/coherency_MR_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namefig));
saveas(fig,fname)
fname = strcat(dir_both_monkeys,sprintf('/coherency_MR_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namepng));
saveas(fig,fname)

% --- ELECTRODE-SENDER coherence   -------%

fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_ms,ctrl_SA.err_ms,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_ms,ctrl_OA.err_ms,'lineprops',{'color',[255, 128, 128]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Senders','Controls-Senders  same area','Controls-Senders  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_both_monkeys,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namefig));
saveas(fig,fname)
fname = strcat(dir_both_monkeys,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namepng));
saveas(fig,fname)

% keyboard 

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

