% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the average coherence for the modulators and controls for
% only a specified number of modulators sorted according to their decoding
% accuracy
%
% @ Gino Del Ferraro, NYU, March 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
monkey = 'Archie';
filename = '_rec001_002.mat'; % -- loading file name 
recording = 'rec001_002_corrected'; % -- folder where to load file
N = 41; % --- max number of modulators 
figstr = sprintf('%d_theta_modulators_removed_Sess_rec001_002',N);
title_caption = sprintf('rec001 rec002 - %d modulators - AUC',N)

dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Controls = strcat(dir_RS,sprintf('/Modulators_Controls_avg_results/%s',recording));
dir_Controls_auc = strcat(dir_Controls,'/AUC');
dir_mod_ctrl_list = strcat(dir_RS,'/Modulators_controls');

mod_list = importdata(strcat(dir_mod_ctrl_list,'/modulators_sorted_decod_accuracy_removed_Sess_AUC_rec001_002.txt'));
display([sprintf('---- > Total number of modulators for %s is : ',monkey),num2str(size(mod_list,1))])
fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

   
% if monkey == 'Maverick'
%     sess_info{1}(20) = [];
%     sess_info{1}(17) = [];
% end 

% -- select the first N index
mod_idx = mod_list(1:N,4);
sess_numb = unique(mod_list(1:N,1)); % session label with top modulators 

% -- find the session index corresponding to the session with top modulators 
% -- exclude bad sessions 
excluded_idx = [2,5,8,9];
sess_info{1}(excluded_idx) = [];

sess_idx = [];
for i=1:length(sess_numb)
    sess_idx = [sess_idx, find(sess_info{1}==sess_numb(i))];
end

fk = 200; W = 5;
% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_Controls,sprintf('/coh_spec_m_fk_%d_W_%d%s',fk,W,filename))); % structure mod
load(strcat(dir_Controls,sprintf('/coh_spec_sr_fk_%d_W_%d%s',fk,W,filename))); % structure stim
mod_mod = mod;
stim_mod = stim;

mod_mod = mod_mod(mod_idx);
stim_mod = stim_mod(sess_idx);

modulators = mean_coh_and_spec_RS(mod_mod,stim_mod);

%%%%%%%%% CONTROLS SAME AREA %%%%%%%%%%%%
ctrl_list = importdata(strcat(dir_mod_ctrl_list,'/control_list_same_area_removed_Sess_rec001_002.txt')); % session, modulator idx, order index i

load(strcat(dir_Controls,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s',fk,W,filename)));
load(strcat(dir_Controls,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s',fk,W,filename)));
mod_ctrl_SA = mod;
stim_ctrl_SA = stim;

ctrl_idx = [];
for i=1:length(sess_numb)
    ctrl_idx = [ctrl_idx; find(ctrl_list(:,1)==sess_numb(i))];
end

mod_ctrl_SA = mod_ctrl_SA(ctrl_idx);
stim_ctrl_SA = stim_ctrl_SA(sess_idx);

ctrl_SA = mean_coh_and_spec_RS(mod_ctrl_SA,stim_ctrl_SA);


%%%%%%%%% CONTROLS OTHER AREAS %%%%%%%%%%%
ctrl_list = importdata(strcat(dir_mod_ctrl_list,'/control_list_other_areas_removed_Sess_rec001_002.txt')); % session, modulator idx, order index i

load(strcat(dir_Controls,sprintf('/coh_spec_m_Controls_other_areas_fk_%d_W_%d%s',fk,W,filename)));
load(strcat(dir_Controls,sprintf('/coh_spec_sr_Controls_other_areas_fk_%d_W_%d%s',fk,W,filename)));
mod_ctrl_OA = mod;
stim_ctrl_OA = stim;

ctrl_idx = [];
for i=1:length(sess_numb)
    ctrl_idx = [ctrl_idx; find(ctrl_list(:,1)==sess_numb(i))];
end

mod_ctrl_OA = mod_ctrl_OA(ctrl_idx);
stim_ctrl_OA = stim_ctrl_OA(sess_idx);

ctrl_OA = mean_coh_and_spec_RS(mod_ctrl_OA,stim_ctrl_OA);

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
title(sprintf('%s: Abs MR coherence MODs vs CTRLs, %s',monkey,title_caption),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Receivers','Controls-Receivers  same area','Controls-Receiver  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Controls_auc,sprintf('/coherency_MR_mods_vs_ctrl_W_%d_fk_%d_%s.png',W,fk,figstr));
saveas(fig,fname)
fname = strcat(dir_Controls_auc,sprintf('/coherency_MR_mods_vs_ctrl_W_%d_fk_%d_%s.fig',W,fk,figstr));
saveas(fig,fname)

% --- ELECTRODE-SENDER coherence   -------%

fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_ms,ctrl_SA.err_ms,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_ms,ctrl_OA.err_ms,'lineprops',{'color',[255, 128, 128]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('%s: Abs MS coherence MODs vs CTRLs, %s ',monkey,title_caption),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Senders','Controls-Senders  same area','Controls-Senders  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_Controls_auc,sprintf('/coherency_MS_mods_vs_ctrl_W_%d_fk_%d_%s.png',W,fk,figstr));
saveas(fig,fname)
fname = strcat(dir_Controls_auc,sprintf('/coherency_MS_mods_vs_ctrl_W_%d_fk_%d_%s.fig',W,fk,figstr));
saveas(fig,fname)


% keyboard 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %      FIGURES  SPECTRUMS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%% SPECTRUMS MODULATORS vs RECEIVERS vs CONTROLS %%%%%%%%%%%%
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
