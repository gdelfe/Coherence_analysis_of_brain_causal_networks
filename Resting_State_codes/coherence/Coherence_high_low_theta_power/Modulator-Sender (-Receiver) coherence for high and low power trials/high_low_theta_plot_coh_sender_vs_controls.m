
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


freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


load(strcat(dir_high_low_theta,'/coh_all_sess_sr_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_low.mat'))


load(strcat(dir_high_low_theta,'/coh_all_sess_controls_sender_SA-receiver.mat'))
% control-sender receiver coherence
coh_cr_SA_high = ctrl_send_SA_coh.cr_high;
coh_cr_SA_low = ctrl_send_SA_coh.cr_low;

load(strcat(dir_high_low_theta,'/coh_all_sess_controls_sender_OA-receiver.mat'))
% control-sender receiver coherence
coh_cr_OA_high = ctrl_send_OA_coh.cr_high;
coh_cr_OA_low = ctrl_send_OA_coh.cr_low;

f = linspace(0,200,409);

% sender and receiver 
mean_all_coh_sr_high = mean(abs(coh_all_c_sr_high),1);
mean_all_coh_sr_low = mean(abs(coh_all_c_sr_low),1);
err_all_coh_sr_high = std(abs(coh_all_c_sr_high),0,1)/sqrt(size(coh_all_c_sr_high,1));
err_all_coh_sr_low = std(abs(coh_all_c_sr_low),0,1)/sqrt(size(coh_all_c_sr_low,1));



% controls SA sender and receveir 
mean_all_coh_cr_SA_high = mean(abs(coh_cr_SA_high),1);
mean_all_coh_cr_SA_low = mean(abs(coh_cr_SA_low),1);
err_all_coh_cr_SA_high = std(abs(coh_cr_SA_high),0,1)/sqrt(size(coh_cr_SA_high,1));
err_all_coh_cr_SA_low = std(abs(coh_cr_SA_low),0,1)/sqrt(size(coh_cr_SA_low,1));

% controls OA sender and receveir 

mean_all_coh_cr_OA_high = mean(abs(coh_cr_OA_high),1);
mean_all_coh_cr_OA_low = mean(abs(coh_cr_OA_low),1);
err_all_coh_cr_OA_high = std(abs(coh_cr_OA_high),0,1)/sqrt(size(coh_cr_OA_high,1));
err_all_coh_cr_OA_low = std(abs(coh_cr_OA_low),0,1)/sqrt(size(coh_cr_OA_low,1));


%%%%%%%%%%%%%%%%%%%%%%
% FIGURES %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% %%%% RECEIVER %%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureVisible','on')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));

fig = figure;
shadedErrorBar(f,mean_all_coh_sr_high,err_all_coh_sr_high,'lineprops',{'color',[255, 179, 102]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_sr_low,err_all_coh_sr_low,'lineprops',{'color',[0, 255, 255]/255 },'patchSaturation',0.5);
shadedErrorBar(f,mean_all_coh_cr_SA_high,err_all_coh_cr_SA_high,'lineprops',{'color',[255, 128, 0]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_cr_SA_low,err_all_coh_cr_SA_low,'lineprops',{'color',[0, 179, 179]/255 },'patchSaturation',0.5);
shadedErrorBar(f,mean_all_coh_cr_OA_high,err_all_coh_cr_OA_high,'lineprops',{'color',[179, 89, 0]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_cr_OA_low,err_all_coh_cr_OA_low,'lineprops',{'color',[0, 102, 102]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title('Coherence SR high vs low theta power trial: sender vs controls SReg/OReg','FontSize',10)
legend('mod high theta pow',' mod low theta pow','CSReg high theta pow',' CSReg low theta pow', 'COReg high theta pow',' COReg low theta pow', 'FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,'/sender-receiver_vs_controls-receiver_coherence_mean.jpg');
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
title('Coherence MS high vs low theta power trial: modulators vs controls SReg/OReg','FontSize',10)
legend('mod high theta pow',' mod low theta pow','CSReg high theta pow',' CSReg low theta pow', 'COReg high theta pow',' COReg low theta pow', 'FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,'/MS_mod_vs_controls_coherence_mean.jpg');
saveas(fig,fname);





