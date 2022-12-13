
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

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';

% coherence for both monkeys for high/low modulator power
both_sr_coh_high = [];
both_sr_coh_low = [];

% %%%%%%%%%%%%%%
% MAVERICK 
% %%%%%%%%%%%%%

monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_Stim_Theta = strcat(dir_main,sprintf('/%s/Stim_data/%s',monkey,freq_band));
dir_high_low_theta = strcat(dir_main,sprintf('/%s/Resting_state/high_low_theta',monkey));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% Sender-Receiver coherence across all brain for high/low modulator theta power 
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_low.mat'))

f = linspace(0,200,409); % frequency range 


% SR coherence phase -- low power
fig_histo = figure;
histogram(angle(coh_all_c_sr_low(:,9:19)),20,'FaceAlpha',.6); grid on
hold on
legend('phase')
title([sprintf("%s",monkey),"Phase dist in theta-freq band -- SR-coh low-power"],'FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')
xticks([-pi, -pi/2, 0, pi/2, pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

fname = strcat(dir_high_low_theta,'/phase_dist_SR_low_pow_all_mod.jpg');
saveas(fig_histo,fname);


% SR coherence -- high power
fig_histo = figure;
histogram(angle(coh_all_c_sr_high(:,9:19)),20,'FaceAlpha',.6); grid on
hold on
legend('phase')
title([sprintf("%s",monkey),"Phase dist in theta-freq band -- SR-coh high-power"],'FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')
xticks([-pi, -pi/2, 0, pi/2, pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

fname = strcat(dir_high_low_theta,'/phase_dist_SR_high_pow_all_mod.jpg');
saveas(fig_histo,fname);



% store Maverick coherence
both_sr_coh_high = [both_sr_coh_high; coh_all_c_sr_high];
both_sr_coh_low = [both_sr_coh_low; coh_all_c_sr_low];


%%%%%%%%%%%%%%%
% ARCHIE 
%%%%%%%%%%%%%%%

monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_Stim_Theta = strcat(dir_main,sprintf('/%s/Stim_data/%s',monkey,freq_band));
dir_high_low_theta = strcat(dir_main,sprintf('/%s/Resting_state/high_low_theta',monkey));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


load(strcat(dir_high_low_theta,'/coh_all_sess_sr_high.mat'))
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_low.mat'))

both_sr_coh_high = [both_sr_coh_high; coh_all_c_sr_high];
both_sr_coh_low = [both_sr_coh_low; coh_all_c_sr_low];

f = linspace(0,200,409); % frequency range 



% SR coherence -- low power
fig_histo = figure;
histogram(angle(coh_all_c_sr_low(:,9:19)),20,'FaceAlpha',.6); grid on
hold on
legend('phase')
title([sprintf("%s",monkey),"Phase dist in theta-freq band -- SR-coh low-power"],'FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')
xticks([-pi, -pi/2, 0, pi/2, pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

fname = strcat(dir_high_low_theta,'/phase_dist_SR_low_pow_all_mod.jpg');
saveas(fig_histo,fname);


% SR coherence -- high power
fig_histo = figure;
histogram(angle(coh_all_c_sr_high(:,9:19)),20,'FaceAlpha',.6); grid on
hold on
legend('phase')
title([sprintf("%s",monkey),"Phase dist in theta-freq band -- SR-coh high-power"],'FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')
xticks([-pi, -pi/2, 0, pi/2, pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

fname = strcat(dir_high_low_theta,'/phase_dist_SR_high_pow_all_mod.jpg');
saveas(fig_histo,fname);




%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOTH MONKEYS
%%%%%%%%%%%%%%%%%%%%%%%%%

dir_both = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys/high_low_theta';
monkey = 'Both monkeys';

% SR coherence -- low power
fig_histo = figure;
histogram(angle(both_sr_coh_low(:,9:19)),20,'FaceAlpha',.6,'Normalization','probability'); grid on
hold on
histogram(angle(both_sr_coh_high(:,9:19)),20,'FaceAlpha',.6,'Normalization','probability'); grid on
legend('phase')
title([sprintf("%s",monkey),"Phase dist in theta-freq band -- SR-coh low-power"],'FontSize',10)
ylabel('count')
xlabel('Phase value (rad)')
legend('low power','high power')
xticks([-pi, -pi/2, 0, pi/2, pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

fname = strcat(dir_both,'/phase_dist_SR_low_and_high.jpg');
saveas(fig_histo,fname);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABS VALUES FOR THE PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % SR coherence -- low power
% fig_histo = figure;
% histogram(abs(angle(both_sr_coh_low(:,9:19))),20,'FaceAlpha',.6); grid on
% hold on
% legend('phase')
% title([sprintf("%s",monkey),"Abs Phase dist in theta-freq band -- SR-coh low-power"],'FontSize',10)
% ylabel('count')
% xlabel('Phase value (rad)')
% xticks([-pi, -pi/2, 0, pi/2, pi])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% 
% fname = strcat(dir_both,'/abs_phase_dist_SR_low_pow_all_mod.jpg');
% saveas(fig_histo,fname);
% 
% 
% fig_histo = figure;
% histogram(abs(mean(angle(both_sr_coh_low(:,9:19)),2)),20,'FaceColor','r','FaceAlpha',.6); grid on
% hold on
% legend('phase')
% title([sprintf("%s",monkey),"Abs Mean-Phase dist in theta-freq band -- SR-coh low-power"],'FontSize',10)
% ylabel('count')
% xlabel('Phase value (rad)')
% xticks([-pi, -pi/2, 0, pi/2, pi])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% 
% fname = strcat(dir_both,'/abs_phase_dist_SR_low_pow_all_mod_mean.jpg');
% saveas(fig_histo,fname);
% 
% 
% % SR coherence -- high power
% fig_histo = figure;
% histogram(abs(angle(both_sr_coh_high(:,9:19))),20,'FaceAlpha',.6); grid on
% hold on
% legend('phase')
% title([sprintf("%s",monkey),"Abs Phase dist in theta-freq band -- SR-coh high-power"],'FontSize',10)
% ylabel('count')
% xlabel('Phase value (rad)')
% xticks([-pi, -pi/2, 0, pi/2, pi])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% 
% fname = strcat(dir_both,'/abs_phase_dist_SR_high_pow_all_mod.jpg');
% saveas(fig_histo,fname);
% 
% 
% fig_histo = figure;
% histogram(abs(mean(angle(both_sr_coh_high(:,9:19)),2)),20,'FaceColor','r','FaceAlpha',.6); grid on
% hold on
% legend('phase')
% title([sprintf("%s",monkey),"Abs Mean-Phase dist in theta-freq band -- SR-coh high-power"],'FontSize',10)
% ylabel('count')
% xlabel('Phase value (rad)')
% xticks([-pi, -pi/2, 0, pi/2, pi])
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
% 
% fname = strcat(dir_both,'/abs_phase_dist_SR_high_pow_all_mod_mean.jpg');
% saveas(fig_histo,fname);
% 
% 

angle_low = angle(both_sr_coh_low(:,9:19))';
angle_high = angle(both_sr_coh_high(:,9:19))';

angle_low = angle_low(:)';
angle_high = angle_high(:)';


figure;
scatter(angle_low,angle_high,25,'o','filled')
xticks([-pi, -pi/2, 0, pi/2, pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yticks([-pi, -pi/2, 0, pi/2, pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
grid on



angle_low = angle(both_sr_coh_low(:,15))';
angle_high = angle(both_sr_coh_high(:,15))';

angle_low = angle_low(:)';
angle_high = angle_high(:)';

figure;
scatter(angle_low,angle_high,25,'o','filled')
xticks([-pi, -pi/2, 0, pi/2, pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
yticks([-pi, -pi/2, 0, pi/2, pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
grid on
