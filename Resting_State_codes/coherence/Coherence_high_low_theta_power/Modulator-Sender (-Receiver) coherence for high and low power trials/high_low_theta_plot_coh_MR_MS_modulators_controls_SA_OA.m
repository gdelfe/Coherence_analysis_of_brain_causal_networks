
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code compute the theta coherence between the controls (Other Area) and
% sender/receiver for modulators having high/low theta power.
%
%    @ Gino Del Ferraro, March 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
dir_high_low_theta = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% load modulators MR and MS
load(strcat(dir_high_low_theta,'/coh_all_sess_ms_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_ms_low.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_mr_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_mr_low.mat'))

load(strcat(dir_high_low_theta,'/coh_all_sess_controls_SA_high-low_pow_controls.mat'))
ctrl_coh_SA = ctrl_coh; 
clear ctrl_coh;

load(strcat(dir_high_low_theta,'/coh_all_sess_controls_OA_high-low_pow_controls.mat'))
ctrl_coh_OA = ctrl_coh; 
clear ctrl_coh;

f = linspace(0,200,409);

% modulator sender and receiver 
mean_all_coh_ms_high = mean(abs(coh_all_c_ms_high),1);
mean_all_coh_ms_low = mean(abs(coh_all_c_ms_low),1);
mean_all_coh_mr_high = mean(abs(coh_all_c_mr_high),1);
mean_all_coh_mr_low = mean(abs(coh_all_c_mr_low),1);

err_all_coh_ms_high = std(abs(coh_all_c_ms_high),0,1)/sqrt(size(coh_all_c_ms_high,1));
err_all_coh_ms_low = std(abs(coh_all_c_ms_low),0,1)/sqrt(size(coh_all_c_ms_low,1));
err_all_coh_mr_high = std(abs(coh_all_c_mr_high),0,1)/sqrt(size(coh_all_c_mr_high,1));
err_all_coh_mr_low = std(abs(coh_all_c_mr_low),0,1)/sqrt(size(coh_all_c_mr_low,1));



% controls SA sender and receveir 
coh_cs_high = ctrl_coh_SA.cs_high; 
coh_cs_low = ctrl_coh_SA.cs_low;
coh_cr_high = ctrl_coh_SA.cr_high;
coh_cr_low = ctrl_coh_SA.cr_low;


mean_all_coh_ms_high_SA = mean(abs(coh_cs_high),1);
mean_all_coh_ms_low_SA = mean(abs(coh_cs_low),1);
mean_all_coh_mr_high_SA = mean(abs(coh_cr_high),1);
mean_all_coh_mr_low_SA = mean(abs(coh_cr_low),1);

err_all_coh_ms_high_SA = std(abs(coh_cs_high),0,1)/sqrt(size(coh_cs_high,1));
err_all_coh_ms_low_SA = std(abs(coh_cs_low),0,1)/sqrt(size(coh_cs_low,1));
err_all_coh_mr_high_SA = std(abs(coh_cr_high),0,1)/sqrt(size(coh_cr_high,1));
err_all_coh_mr_low_SA = std(abs(coh_cr_low),0,1)/sqrt(size(coh_cr_low,1));

clear coh_cs_high coh_cs_low coh_cr_high coh_cr_low

% controls OA sender and receveir 
coh_cs_high = ctrl_coh_OA.cs_high; 
coh_cs_low = ctrl_coh_OA.cs_low;
coh_cr_high = ctrl_coh_OA.cr_high;
coh_cr_low = ctrl_coh_OA.cr_low;


mean_all_coh_ms_high_OA = mean(abs(coh_cs_high),1);
mean_all_coh_ms_low_OA = mean(abs(coh_cs_low),1);
mean_all_coh_mr_high_OA = mean(abs(coh_cr_high),1);
mean_all_coh_mr_low_OA = mean(abs(coh_cr_low),1);

err_all_coh_ms_high_OA = std(abs(coh_cs_high),0,1)/sqrt(size(coh_cs_high,1));
err_all_coh_ms_low_OA = std(abs(coh_cs_low),0,1)/sqrt(size(coh_cs_low,1));
err_all_coh_mr_high_OA = std(abs(coh_cr_high),0,1)/sqrt(size(coh_cr_high,1));
err_all_coh_mr_low_OA = std(abs(coh_cr_low),0,1)/sqrt(size(coh_cr_low,1));



%%%%%%%%%%%%%%%%%%%%%%
% FIGURES %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% %%%% RECEIVER %%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureVisible','on')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));

fig = figure;
shadedErrorBar(f,mean_all_coh_mr_high,err_all_coh_mr_high,'lineprops',{'color',[255, 179, 102]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_mr_low,err_all_coh_mr_low,'lineprops',{'color',[0, 255, 255]/255 },'patchSaturation',0.5);
shadedErrorBar(f,mean_all_coh_mr_high_SA,err_all_coh_mr_high_SA,'lineprops',{'color',[255, 128, 0]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_mr_low_SA,err_all_coh_mr_low_SA,'lineprops',{'color',[0, 179, 179]/255 },'patchSaturation',0.5);
shadedErrorBar(f,mean_all_coh_mr_high_OA,err_all_coh_mr_high_OA,'lineprops',{'color',[179, 89, 0]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_mr_low_OA,err_all_coh_mr_low_OA,'lineprops',{'color',[0, 102, 102]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title(["Coherence MR high vs low theta power trial (of the electrode):", "modulators vs controls SReg/OReg"],'FontSize',10)
legend('mod high theta pow',' mod low theta pow','CSReg high theta pow',' CSReg low theta pow', 'COReg high theta pow',' COReg low theta pow', 'FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,'/MR_mod_vs_controls_coherence_mean_high-low_pow_electrode.jpg');
saveas(fig,fname);



% %%%% SENDER %%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureVisible','on')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));

fig = figure;
shadedErrorBar(f,mean_all_coh_ms_high,err_all_coh_ms_high,'lineprops',{'color',[194, 102, 255]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_ms_low,err_all_coh_ms_low,'lineprops',{'color',[255, 102, 102]/255 },'patchSaturation',0.5);
shadedErrorBar(f,mean_all_coh_ms_high_SA,err_all_coh_ms_high_SA,'lineprops',{'color',[153, 0, 255]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_ms_low_SA,err_all_coh_ms_low_SA,'lineprops',{'color',[255, 26, 26]/255 },'patchSaturation',0.5);
shadedErrorBar(f,mean_all_coh_ms_high_OA,err_all_coh_ms_high_OA,'lineprops',{'color',[92, 0, 153]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_ms_low_OA,err_all_coh_ms_low_OA,'lineprops',{'color',[179, 0, 0]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title(["Coherence MS high vs low theta power trial (of the electrode):", "modulators vs controls SReg/OReg"],'FontSize',10)
legend('mod high theta pow',' mod low theta pow','CSReg high theta pow',' CSReg low theta pow', 'COReg high theta pow',' COReg low theta pow', 'FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,'/MS_mod_vs_controls_coherence_mean_high-low_pow_electrode.jpg');
saveas(fig,fname);





