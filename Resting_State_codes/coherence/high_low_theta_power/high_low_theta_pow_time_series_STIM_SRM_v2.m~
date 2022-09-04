
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the time series (LFP) for low and high theta power of the
% modulator for the STIM experiment. High theta power time series should show oscillatory rhythms
% at the theta frequency.
%
%    @ Gino Del Ferraro, June 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
% close all;

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

name_struct_input = '/session_data_lfp.mat'; % -- name file to load

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_Stim_Theta = strcat(dir_main,sprintf('/%s/Stim_data/%s',monkey,freq_band));
dir_sort = strcat(dir_main,sprintf('/%s/Resting_state/%s/Modulators_controls',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

modulators = importdata(strcat(dir_sort,'/modulators_sorted_decod_accuracy_v2.txt')); % session, modulator idx, decod accuracy, order index i

% Select specific Session to look at
s = 23;

for s = 11
Sess = sess_info{1}(s); % Session number

dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info


% Print modulators brain regions
unique(sess_data_lfp.mod_areas)

% select modulator:
cnt_m = 1; % modulator index
for cnt_m = 1 % 1 : length(sess_data_lfp.mod_idx)

    el = sess_data_lfp.mod_idx(cnt_m); % electrode
    lfp_m = sq(sess_data_lfp.lfp_E(el,:,:));

    % Select high and low theta power trials
    W = 3;
    [spec, f, err] = dmtspec(lfp_m,[1000/1e3,W],1e3,200);

    theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
    theta_pow_mean = mean(theta_pow); % get the average theta power
    theta_pow = theta_pow - theta_pow_mean; % rescale the theta power by removing the mean value

    [sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
    cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution

    % low and high theta power indexes
    low_theta = sort_theta(1:cut);
    low_idx = trial_idx(1:cut);

    high_theta = sort_theta(end-cut+1:end);
    high_idx = trial_idx(end-cut+1:end);
end



lfp_S = sess_data_lfp.lfp_S;
lfp_R = sess_data_lfp.lfp_R;


% Select trials based on threshold 
% [lfp_m_H,idx_H] = amplitude_H(lfp_m(high_idx,:),100);
% [lfp_m_L,idx_L] = amplitude_L(lfp_m(low_idx,:),30);
% 
% Xm = lfpfilter_low(lfp_m_H,1000);
% filter_m = thetafilter(Xm,800);
% Xm_smooth = Xm-smooth_array(Xm,100);
% 
% for i=1:15
% figure;
% plot(Xm_smooth(i,:)); hold on
% plot(filter_m(i,:))
% end


% select high and low power trials and then merge in a unique time series
% for visual investigation 

% High power 

% modulator
lfp_m_H = lfp_m(high_idx,:);
Xm = lfpfilter_low(lfp_m_H,1000);
filter_m = thetafilter(Xm,800);
Xm_smooth = Xm-smooth_array(Xm,100);
lfp_T = Xm_smooth';
Xm_smooth = lfp_T(:)';
lfp_T = filter_m';
filter_m = lfp_T(:)';

% sender
lfp_S_H = lfp_S(high_idx,:);
Xs = lfpfilter_low(lfp_S_H,1000);
filter_s = thetafilter(Xs,800);
Xs_smooth = Xs-smooth_array(Xs,100);
lfp_T = Xs_smooth';
Xs_smooth = lfp_T(:)';
lfp_T = filter_s';
filter_s = lfp_T(:)';

% receiver 
lfp_R_H = lfp_R(high_idx,:);
Xr = lfpfilter_low(lfp_R_H,1000);
filter_r = thetafilter(Xr,800);
Xr_smooth = Xr-smooth_array(Xr,100);
lfp_T = Xr_smooth';
Xr_smooth = lfp_T(:)';
lfp_T = filter_r';
filter_r = lfp_T(:)';

a = 30;
b = 70;
figure;
plot(Xm_smooth); hold on
plot(filter_m); hold on
plot(Xr_smooth + a); grid on
plot(filter_r + a); hold on
plot(Xs_smooth + b); hold on
plot(filter_s + b); hold on
% xlim([xa xb]);
grid on
title(sprintf('Sess = %d high power',s),'Fontsize',12);
legend('modulator','m theta','receiver','r theta','sender','s theta');


% close all 
% xa = 47000;
% xb = 48000;
% ya = -60;
% yb = 65;
% fig_ts = figure;
% ha = tight_subplot(3,1,[.04 .04],[.2 .01],[.1 .1])
% 
% axes(ha(1))
% plot(Xm_smooth); hold on
% plot(filter_m); 
% ylabel('Modulator','FontName','Arial','FontSize',12);
% xlim([xa xb])
% ylim([ya yb])
% 
% axes(ha(2))
% plot(Xr_smooth); hold on
% plot(filter_r); 
% ylabel('Receiver','FontName','Arial','FontSize',12);
% xlim([xa xb])
% ylim([ya yb])
% 
% axes(ha(3))
% plot(Xs_smooth); hold on
% plot(filter_s); 
% ylabel('Sender','FontName','Arial','FontSize',12);
% 
% xlim([xa xb])
% ylim([ya yb])
% set(gcf, 'Position',  [100, 600, 320, 1050]);
% 
% fname = strcat(dir_RS_Theta,sprintf('/time_series/SRM/lfp_SRM_%d_sess_%d_ts_47_48.pdf',cnt_m,Sess))
% saveas(fig_ts,fname);


% Low power

lfp_m_L = lfp_m(low_idx,:);
Xm = lfpfilter_low(lfp_m_L,1000);
filter_m = thetafilter(Xm,800);
Xm_smooth = Xm-smooth_array(Xm,100);
lfp_T = Xm_smooth';
Xm_smooth = lfp_T(:)';
lfp_T = filter_m';
filter_m = lfp_T(:)';

% sender
lfp_S_L = lfp_S(low_idx,:);
Xs = lfpfilter_low(lfp_S_L,1000);
filter_s = thetafilter(Xs,800);
Xs_smooth = Xs-smooth_array(Xs,100);
lfp_T = Xs_smooth';
Xs_smooth = lfp_T(:)';
lfp_T = filter_s';
filter_s = lfp_T(:)';

% receiver 
lfp_R_L = lfp_R(low_idx,:);
Xr = lfpfilter_low(lfp_R_L,1000);
filter_r = thetafilter(Xr,800);
Xr_smooth = Xr-smooth_array(Xr,100);
lfp_T = Xr_smooth';
Xr_smooth = lfp_T(:)';
lfp_T = filter_r';
filter_r = lfp_T(:)';

close all 
xa = 15000;
xb = 16000;
ya = -30;
yb = 30;
fig_ts = figure;
ha = tight_subplot(3,1,[.04 .04],[.2 .01],[.1 .1])

axes(ha(1))
plot(Xm_smooth); hold on
plot(filter_m); 
ylabel('Modulator','FontName','Arial','FontSize',12);
xlim([xa xb])
ylim([ya yb])

axes(ha(2))
plot(Xr_smooth); hold on
plot(filter_r); 
ylabel('Receiver','FontName','Arial','FontSize',12);
xlim([xa xb])
ylim([ya yb])

axes(ha(3))
plot(Xs_smooth); hold on
plot(filter_s); 
ylabel('Sender','FontName','Arial','FontSize',12);

xlim([xa xb])
ylim([ya yb])
set(gcf, 'Position',  [100, 600, 320, 1050]);

fname = strcat(dir_RS_Theta,sprintf('/time_series/SRM/lfp_low_pow_SRM_%d_sess_%d_ts_15_16.pdf',cnt_m,Sess))
saveas(fig_ts,fname);





end 


