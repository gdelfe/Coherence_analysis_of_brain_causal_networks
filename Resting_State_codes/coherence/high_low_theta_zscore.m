
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the theta coherence between the modulator-sender and
% modulator-receiver for modulators having high/low theta power and
% z-scores the results by using the data from the  permutation test used to
% create a null distribution for the zero-coherence hypothesis 
%
%    @ Gino Del Ferraro, November 2020, Pesaran lab, NYU
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


% load MR and MS coherence for low and high theta power 
load(strcat(dir_high_low_theta,'/coh_all_sess_ms_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_ms_low.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_mr_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_mr_low.mat'))

% load MR and MS permuted coherences for all the sessions and all the modulators
load(strcat(dir_high_low_theta,'/coh_all_permuted_c_ms.mat'))
load(strcat(dir_high_low_theta,'/coh_all_permuted_c_mr.mat'))

% mean of the null distributions for MR and MS
% mean_c_ms_null = sq(mean(mean(cc)))';
% mean_c_mr_null = sq(mean(mean(abs(coh_all_c_mr))))';
% 
% coh_ms_reshape = reshape(coh_all_c_ms,[],409);
% coh_mr_reshape = reshape(coh_all_c_mr,[],409);
% 
% mean_c_ms_null = mean(abs(coh_ms_reshape));
% mean_c_mr_null = mean(abs(coh_mr_reshape));
% 
% std_c_ms_null = std(abs(coh_ms_reshape));
% std_c_mr_null = std(abs(coh_mr_reshape));



% mean of the null distribution for MS and MR for each mod separately
mean_c_ms_null = sq(mean(abs(coh_all_c_ms),2)); % 109x409  = #mod x freq
mean_c_mr_null = sq(mean(abs(coh_all_c_mr),2));

% std of the null distribution for MS and MR for each mod separately
std_c_ms_null = sq(std(abs(coh_all_c_ms),0,2));
std_c_mr_null = sq(std(abs(coh_all_c_mr),0,2));

% One unique mean and std for the null distribution 
mean_all_null_c_ms = mean(mean_c_ms_null,1);
mean_all_null_c_mr = mean(mean_c_mr_null,1);
std_all_null_c_ms = mean(std_c_ms_null,1);
std_all_null_c_mr = mean(std_c_mr_null,1);


zscore_c_ms_high = (abs(coh_all_c_ms_high) - mean_all_null_c_ms)./std_all_null_c_ms;
zscore_c_ms_low = (abs(coh_all_c_ms_low) - mean_all_null_c_ms)./std_all_null_c_ms;
zscore_c_mr_high = (abs(coh_all_c_mr_high) - mean_all_null_c_mr)./std_all_null_c_mr;
zscore_c_mr_low = (abs(coh_all_c_mr_low) - mean_all_null_c_mr)./std_all_null_c_mr;

% zscore_c_ms_high = (abs(coh_all_c_ms_high) - mean_c_ms_null)./std_c_ms_null;
% zscore_c_ms_low = (abs(coh_all_c_ms_low) - mean_c_ms_null)./std_c_ms_null;
% zscore_c_mr_high = (abs(coh_all_c_mr_high) - mean_c_mr_null)./std_c_mr_null;
% zscore_c_mr_low = (abs(coh_all_c_mr_low) - mean_c_mr_null)./std_c_mr_null;


err_coh_ms_high = std(zscore_c_ms_high)/sqrt(size(coh_all_c_ms_high,1));
err_coh_ms_low = std(zscore_c_ms_low)/sqrt(size(coh_all_c_ms_low,1));
err_coh_mr_high = std(zscore_c_mr_high)/sqrt(size(coh_all_c_mr_high,1));
err_coh_mr_low = std(zscore_c_mr_low)/sqrt(size(coh_all_c_mr_low,1));

f = linspace(0,200,409);

fig = figure;
plot(f,mean_all_null_c_ms)
hold on
plot(f,mean_all_null_c_mr)
hold on
plot(f,mean(abs(coh_all_c_mr_high)))
hold on
plot(f,mean(abs(coh_all_c_mr_low)))
hold on
plot(f,mean(abs(coh_all_c_ms_high)))
hold on
plot(f,mean(abs(coh_all_c_ms_low)))
legend('null MS', 'null MR','MR high theta','MR low theta','MS high theta','MS low theta');
title('coherence vs frequency - high/low theta and null')
grid on 
set(gcf, 'Position',  [100, 600, 898, 500])

fname = strcat(dir_high_low_theta,'/All_high_low_coherences_and_null.jpg');
saveas(fig,fname);


figure;
plot(f,mean(zscore_c_ms_high))
hold on 
plot(f,mean(zscore_c_ms_low))




%%%%%%%%%%%%%%%%%%%%%%
% FIGURES %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% %%%% RECEIVER %%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureVisible','on')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));

fig = figure;
shadedErrorBar(f,mean(zscore_c_mr_high),err_coh_mr_high,'lineprops',{'color',[28 199 139]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean(zscore_c_mr_low),err_coh_mr_low,'lineprops',{'color',[50 250 93]/255 },'patchSaturation',0.5);


grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('z-scored Coherence','FontName','Arial','FontSize',15);
title('Coherence MR high vs low theta power trial','FontSize',12)
legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,'/MR_zscored_coherence_mean.jpg');
saveas(fig,fname);


% %%%% SENDER %%%%%%%%%%%%%%%%%%%%%%


fig = figure;
hold all

shadedErrorBar(f,mean(zscore_c_ms_high),err_coh_ms_high,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.5); hold on
shadedErrorBar(f,mean(zscore_c_ms_low),err_coh_ms_low,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('z-scored Coherence','FontName','Arial','FontSize',15);
legend('Modulator - Receiver','Modulator - Sender','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.36])
grid on


fname = strcat(dir_high_low_theta,'/MS_zscored_coherence_mean.jpg');
saveas(fig,fname);






