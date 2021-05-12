

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
freq_modulators = 'theta band';

% maverick
filename_mav = '_movie';
recording_mav = 'movie';
% archie
filename_arc = '_movie_all';
recording_arc = 'movie_all_Sess'
% both
title_rec = 'movie (all Sess)';
both_fname = '_movie_all';

dir_Maverick = strcat(dir_main,sprintf('Maverick/Resting_state/%s/Modulators_Controls_avg_results/%s',freq_band,recording_mav));
dir_Archie = strcat(dir_main,sprintf('Archie/Resting_state/%s/Modulators_Controls_avg_results/%s',freq_band,recording_arc));


dir_out = strcat(dir_main,sprintf('Maverick/Resting_state/%s/both_monkeys/modulators_sender_receiver',freq_band));


fk = 200; W = 5;
%%%%%%%%% MODULATORS Maverick %%%%%%
load(strcat(dir_Maverick,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,filename_mav))); % structure mod
load(strcat(dir_Maverick,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,filename_mav))); % structure stim
struct_mod_mav = mod;
struct_stim_mav = stim;

% data_last = mean_coh_and_spec_RS(mod,stim);


% %%%%%%%%% MODULATORS Archie %%%%%%
load(strcat(dir_Archie,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,filename_arc))); % structure mod
load(strcat(dir_Archie,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,filename_arc))); % structure stim
struct_mod_arc = mod;
struct_stim_arc = stim;

modulators = mean_coh_and_spec_across_monkeys(struct_mod_mav,struct_stim_mav,struct_mod_arc,struct_stim_arc);

% data_1 = mean_coh_and_spec_RS(mod,stim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,fk,size(modulators.mean_coh_ms,2)); % frequency values (range)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES MR, MS, SR %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all

shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
% % shadedErrorBar(f,modulators.mean_coh_sr,modulators.err_sr,'lineprops',{'color',[0 204 204]/255},'patchSaturation',0.4); hold on

% shadedErrorBar(f,data_last.mean_coh_ms,data_last.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% shadedErrorBar(f,data_1.mean_coh_ms,data_1.err_ms,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on

% shadedErrorBar(f,data_1.mean_coh_ms,data_1.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% shadedErrorBar(f,data_1.mean_coh_mr,data_1.err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('Both animals - %s - Abs coherency of MR, MS, - %s',freq_modulators,title_rec),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
% legend('M-S abs coherency','M-R abs coherency','S-R abs coherency','FontSize',10)
legend('M-S abs coherency','M-R abs coherency','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on


keyboard 

% fname = strcat(dir_out,sprintf('/coherency_MS_MR_SR_W_%d_fk_%d_%s_both_monkeys.png',W,fk,freq_band));
% saveas(fig,fname)
% fname = strcat(dir_out,sprintf('/coherency_MS_MR_SR_W_%d_fk_%d_%s_both_monkeys.fig',W,fk,freq_band));
% saveas(fig,fname)



fname = strcat(dir_out,sprintf('/coherency_MS_MR_W_%d_fk_%d_%s_both_monkeys%s.png',W,fk,freq_band,both_fname));
saveas(fig,fname)
fname = strcat(dir_out,sprintf('/coherency_MS_MR_W_%d_fk_%d_%s_both_monkeys%s.fig',W,fk,freq_band,both_fname));
saveas(fig,fname)












