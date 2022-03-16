
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code ..
%
%    @ Gino Del Ferraro, November 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - PATHS and NAMES --- %
%%%%%%%%%%%%%%%%%%%
addpath('/mnt/pesaranlab/Matlab/monkeys')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
dir_high_low_theta = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_Stim_Theta = strcat(dir_main,sprintf('/%s/Stim_data/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


load(strcat(dir_high_low_theta,'/coh_all_sess_ms_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_ms_low.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_mr_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_mr_low.mat'))

f = linspace(0,200,409); % frequency range 


% MS coherence -- low power
fig_histo = figure;
histogram(angle(coh_all_c_ms_low(:,9:19)),20,'FaceAlpha',.6); grid on
hold on
legend('phase')
title('Phase dist in theta-freq band -- MS-coh low-power','FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')

fname = strcat(dir_high_low_theta,'/phase_dist_MS_low_pow_all_mod.jpg');
saveas(fig_histo,fname);


fig_histo = figure;
histogram(mean(angle(coh_all_c_ms_low(:,9:19)),2),20,'FaceColor','r','FaceAlpha',.6); grid on
hold on
legend('phase')
title('Mean-Phase dist in theta-freq band -- MS-coh low-power','FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')

fname = strcat(dir_high_low_theta,'/phase_dist_MS_low_pow_all_mod_mean.jpg');
saveas(fig_histo,fname);


% MS coherence -- high power
fig_histo = figure;
histogram(angle(coh_all_c_ms_high(:,9:19)),20,'FaceAlpha',.6); grid on
hold on
legend('phase')
title('Phase dist in theta-freq band -- MS-coh high-power','FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')

fname = strcat(dir_high_low_theta,'/phase_dist_MS_high_pow_all_mod.jpg');
saveas(fig_histo,fname);


fig_histo = figure;
histogram(mean(angle(coh_all_c_ms_high(:,9:19)),2),20,'FaceColor','r','FaceAlpha',.6); grid on
hold on
legend('phase')
title('Mean-Phase dist in theta-freq band -- MS-coh high-power','FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')

fname = strcat(dir_high_low_theta,'/phase_dist_MS_high_pow_all_mod_mean.jpg');
saveas(fig_histo,fname);


% MR coherence -- low power
fig_histo = figure;
histogram(angle(coh_all_c_mr_low(:,9:19)),20,'FaceAlpha',.6); grid on
hold on
legend('phase')
title('Phase dist in theta-freq band -- MR-coh low-power','FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')

fname = strcat(dir_high_low_theta,'/phase_dist_MR_low_pow_all_mod.jpg');
saveas(fig_histo,fname);


fig_histo = figure;
histogram(mean(angle(coh_all_c_mr_low(:,9:19)),2),20,'FaceColor','r','FaceAlpha',.6); grid on
hold on
legend('phase')
title('Mean-Phase dist in theta-freq band -- MR-coh low-power','FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')

fname = strcat(dir_high_low_theta,'/phase_dist_MR_low_pow_all_mod_mean.jpg');
saveas(fig_histo,fname);


% MR coherence -- high power
fig_histo = figure;
histogram(angle(coh_all_c_mr_high(:,9:19)),20,'FaceAlpha',.6); grid on
hold on
legend('phase')
title('Phase dist in theta-freq band -- MR-coh high-power','FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')

fname = strcat(dir_high_low_theta,'/phase_dist_MR_high_pow_all_mod.jpg');
saveas(fig_histo,fname);


fig_histo = figure;
histogram(mean(angle(coh_all_c_mr_high(:,9:19)),2),20,'FaceColor','r','FaceAlpha',.6); grid on
hold on
legend('phase')
title('Mean-Phase dist in theta-freq band -- MR-coh high-power','FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')

fname = strcat(dir_high_low_theta,'/phase_dist_MR_high_pow_all_mod_mean.jpg');
saveas(fig_histo,fname);



