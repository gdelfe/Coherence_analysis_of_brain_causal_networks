
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;  close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)



addpath('/vol/bd5/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/vol/bd5/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
monkey = 'Archie';
filename = '_movie_all';
recording = 'movie';
title_rec = 'movie - ALL Sessions';
freq_modulators = 'theta modulators';


dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Controls = strcat(dir_RS,sprintf('/Modulators_Controls_avg_results/%s',recording));


fk = 200; W = 5;
% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_Controls,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,filename))); % structure mod
load(strcat(dir_Controls,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,filename))); % structure stim
stim_mod = stim;
mod_mod = mod;

% if monkey == 'Archie' 
%     % -- exclude bad sessions
%     excluded_sess = [8,22,30,31];
%     excluded_idx = [2,5,8,9];
%     bad_idx = importdata(strcat(dir_RS,'/Modulators_Controls/Archie_modulator_idx_bad_data.txt'));
%     mod(bad_idx) = [];
%     stim(excluded_idx) = [];
% end

modulators = mean_coh_and_spec_RS(mod,stim);

f = linspace(1,fk,size(modulators.mean_coh_ms,2)); % frequency values (range)


set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, SR, MS
fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all

shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,modulators.mean_coh_sr,modulators.err_sr,'lineprops',{'color',[0 204 204]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('%s - %s - Abs coherency of MR, MS - %s',monkey,freq_modulators,title_rec),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('M-S abs coherency','M-R abs coherency','S-R abs coherency','FontSize',10)
% % % legend('M-S abs coherency','M-R abs coherency','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% xlim([0 60])
set(gcf, 'Position',  [100, 600, 1000, 600])

keyboard

fname = strcat(dir_Controls,sprintf('/coherency_MS_MR_SR_W_%d_fk_%d_%s%s.png',W,fk,freq_band,filename));
saveas(fig,fname)
fname = strcat(dir_Controls,sprintf('/coherency_MS_MR_SR_W_%d_fk_%d_%s%s.fig',W,fk,freq_band,filename));
saveas(fig,fname)



fname = strcat(dir_Controls,sprintf('/coherency_MS_MR_W_%d_fk_%d_%s%s.png',W,fk,freq_band,filename));
saveas(fig,fname)
fname = strcat(dir_Controls,sprintf('/coherency_MS_MR_W_%d_fk_%d_%s%s.fig',W,fk,freq_band,filename));
saveas(fig,fname)





