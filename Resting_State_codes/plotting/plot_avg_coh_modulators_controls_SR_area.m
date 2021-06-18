% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the average coherence for the modulators and controls for
% only a specific sender-receiver edge
%
% @ Gino Del Ferraro, NYU, June 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

sess_list_file = '/Sessions_with_modulator_info_movie.txt';
freq_band = 'theta_band';
freq = 'theta band';
monkey = 'Maverick';

filename_mod = ''; % -- loading file name for coherence averages ******************
filename_ctrl = ''; % -- loading file name for the list the coherences in sess_data_lfp_coherence
title_caption = 'S:OFC - R:CN'; % -- title caption 
SR_brain_areas = 'OFC_CN'; % -- name of SR brain area for the figures and coherence files 


recording = 'last_recording'; % -- folder where to load coherency files  *************

dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS,sess_list_file)); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

fk = 200; W = 5;
sess_list_idx = [1,2,5,6,9,24];
mod = [];
ctrl_SA = [];
ctrl_OA = [];

for i = sess_list_idx
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess)]);
    
    dir_Modulators = strcat(dir_RS,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    load(strcat(dir_Modulators,sprintf('/sess_data_lfp_coherence_fk_200_W_5%s.mat',filename_mod)));
    
    mod = [mod, sess_data_lfp.mod];
    
    dir_ctrl_SA = strcat(dir_RS,sprintf('/Sess_%d/Controls_same_area/%s',Sess,recording));
    load(strcat(dir_ctrl_SA,sprintf('/sess_data_lfp_coherence_fk_200_W_5%s.mat',filename_ctrl)));
    
    ctrl_SA = [ctrl_SA, sess_control_lfp.ctrl]; % stack up control same area structure for coherence
    
    
    dir_ctrl_OA = strcat(dir_RS,sprintf('/Sess_%d/Controls_other_areas/%s',Sess,recording));
    load(strcat(dir_ctrl_OA,sprintf('/sess_data_lfp_coherence_fk_200_W_5%s.mat',filename_ctrl)));
    
    ctrl_OA = [ctrl_OA, sess_control_lfp.ctrl];  % stack up control other areas structure for coherence
    
    
end

modulators = mean_coh_MR(mod);
ctrl_SA_avg = mean_coh_MR(ctrl_SA);
ctrl_OA_avg = mean_coh_MR(ctrl_OA);

dir_output = strcat(dir_RS,sprintf('/Modulators_Controls_avg_results/%s/by_SR_brain_area/%s',recording,SR_brain_areas));
if ~exist(dir_output, 'dir')
    mkdir(dir_output)
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES MODULATORS vs CONTROLS %%%%%%%%%%%%

f = linspace(1,fk,size(modulators.mean_coh_ms,2)); % frequency values (range)

set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA_avg.mean_coh_mr,ctrl_SA_avg.err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA_avg.mean_coh_mr,ctrl_OA_avg.err_mr,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('%s: %s, MR coherence MODs vs CTRLs, %s',monkey,freq,title_caption),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Receivers','Controls-Receivers  same area','Controls-Receiver  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_output,sprintf('/coherency_MR_mods_vs_ctrl_W_%d_fk_%d_%s.png',W,fk,SR_brain_areas));
saveas(fig,fname)
fname = strcat(dir_output,sprintf('/coherency_MR_mods_vs_ctrl_W_%d_fk_%d_%s.fig',W,fk,SR_brain_areas));
saveas(fig,fname)

% --- ELECTRODE-SENDER coherence   -------%

fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA_avg.mean_coh_ms,ctrl_SA_avg.err_ms,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA_avg.mean_coh_ms,ctrl_OA_avg.err_ms,'lineprops',{'color',[255, 128, 128]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('%s: %s, MS coherence MODs vs CTRLs, %s ',monkey,freq,title_caption),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Senders','Controls-Senders  same area','Controls-Senders  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_output,sprintf('/coherency_MS_mods_vs_ctrl_W_%d_fk_%d_%s.png',W,fk,SR_brain_areas));
saveas(fig,fname)
fname = strcat(dir_output,sprintf('/coherency_MS_mods_vs_ctrl_W_%d_fk_%d_%s.fig',W,fk,SR_brain_areas));
saveas(fig,fname)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------  SAVE STRUCTURES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save(strcat(dir_output,sprintf('/modulators_%s.mat',SR_brain_areas)),'modulators');
save(strcat(dir_output,sprintf('/controls_same_area_%s.mat',SR_brain_areas)),'ctrl_SA');
save(strcat(dir_output,sprintf('/controls_other_areas_%s.mat',SR_brain_areas)),'ctrl_OA');

