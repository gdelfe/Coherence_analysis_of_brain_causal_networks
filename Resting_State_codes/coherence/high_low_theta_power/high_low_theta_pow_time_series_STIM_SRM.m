
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

for s = 10
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
%
%     % Power trials
%     lfp_E_high = lfp_m(high_idx,:); % high theta
%     lfp_E_low = lfp_m(low_idx,:); % low theta
%
%     % plot theta power histogram
%     % fig = figure;
%     % histogram(high_theta,30,'FaceAlpha',0.5)
%     % hold on
%     % histogram(low_theta,30,'FaceAlpha',0.5)
%
%     fig = figure;
%     histogram(theta_pow,30,'FaceAlpha',0.5,'Normalization','probability','EdgeColor','w')
%     grid on
%     hold on; xline(low_theta(end),'b');
%     hold on; xline(high_theta(1),'r');
%     legend('theta power','low theta','high theta')
%
%
%     lfp_S = sess_data_lfp.lfp_S;
%     lfp_R = sess_data_lfp.lfp_R;
%
%     lfp_S_high = lfp_S(high_idx,:);
%     lfp_R_high = lfp_R(high_idx,:);
%
%     lfp_S_low = lfp_S(low_idx,:);
%     lfp_R_low = lfp_R(low_idx,:);
%
%
%     % plot LFP for a few trials
%
%         hm = 47; hs = 47; hr = 47;
%         lm = 3; ls = 3; lr = 3;
%
%         a = 45; b = 5;
%         hm = a; hs = a; hr = a;
%         lm = b; ls = b; lr = b;
%
%
%         fig_ts = figure;
%         [ha] = tight_subplot(6,1,[.02 .02],[.05 .05],[.1 .1]);
%
%         for i = 1:1
%             % %%%%%%%%%%% high power
%             axes(ha(1));
%             plot(lfp_E_high(hm,:),'color',[0 0.4470 0.7410]) % plot trials on the right tail of the distribution
%             grid on
%             axes(ha(2));
%             plot(lfp_S_high(hs,:),'color',[0.4940 0.1840 0.5560]) % plot trials on the right tail of the distribution
%             grid on
%             axes(ha(3));
%             plot(lfp_R_high(hr,:),'color',[0.4660 0.6740 0.1880]) % plot trials on the right tail of the distribution
%             grid on
%             %set(gca,'xticklabel',[]);
%             ylim([-330,330])
%
%             % %%%%%%%%%%% low power
%             axes(ha(end-2));
%             plot(lfp_E_low(lm,:),'color',[0 0.4470 0.7410]) % plot trials on the left tail of the distribution
%             grid on;
%             axes(ha(end-1));
%             plot(lfp_S_low(ls,:),'color',[0.4940 0.1840 0.5560]) % plot trials on the right tail of the distribution
%             grid on;
%             axes(ha(end));
%             plot(lfp_R_low(lr,:),'color',[0.4660 0.6740 0.1880]) % plot trials on the right tail of the distribution
%             grid on
%             %             ylim([-330,330])
%             j = j + 1;
%             % set(ha(1:5),'XTickLabel',''); set(ha,'YTickLabel','')
%
%             %     ylabel('Lfp','FontName','Arial','FontSize',12);
%             set(gcf, 'Position',  [100, 600, 700, 1000]);
%             % %     save the plots
%             dir_ts_rec_fig = strcat(dir_Stim_Theta,sprintf('/high_low_theta',Sess));
%
%         end
%
%         % save the plots
%         dir_ts_rec_fig = strcat(dir_high_low_theta);
%         fname = strcat(dir_ts_rec_fig,sprintf('/time_series_MSR_high_low_power.jpg',cnt_m));
%         saveas(fig_ts,fname);
%         fname = strcat(dir_ts_rec_fig,sprintf('/time_series_MSR_high_low_power.pdf',cnt_m));
%         saveas(fig_ts,fname);
%
%
%
% end

lfp_S = sess_data_lfp.lfp_S;
lfp_R = sess_data_lfp.lfp_R;

i = 1;
m = sess_data_lfp.mod_idx(cnt_m); % electrode
lfp_m = sq(sess_data_lfp.lfp_E(m(i),:,:));
lfp_T = lfp_m';
Xm = lfp_T(:)';
Xm = lfpfilter_low(Xm,1000);

filter_m = thetafilter(Xm,800);

lfp_T = lfp_S';
Xs = lfp_T(:)';
Xs = lfpfilter_low(Xs,1000);

filter_s = thetafilter(Xs,800);


lfp_T = lfp_R';
Xr = lfp_T(:)';
Xr = lfpfilter_low(Xr,1000);

filter_r = thetafilter(Xr,800);

xa = 86*1000;
xb = 89*1000;
a = 50;
b = 100;

figure;
plot(Xm-smooth(Xm,100)); hold on
plot(filter_m); hold on
plot(Xr-smooth(Xr,100) + a); grid on
plot(filter_r + a); hold on
plot(Xs-smooth(Xs,100) + b); hold on
plot(filter_s + b); hold on
xlim([xa xb]);
grid on
title(sprintf('Sess = %d',s),'Fontsize',12);
legend('modulator','m theta','receiver','r theta','sender','s theta');


figure;
a = 5;
b = 10;
plot(filter_m); hold on
plot(filter_r+a); hold on
plot(filter_s + b); hold on
% xlim([61000 68000]);
xlim([xa xb]);
title(sprintf('Sess = %d',s),'Fontsize',12);
legend('m theta','r theta','s theta');
grid on

%%%%%%%%%%%%%%%%%%%%%%%
% LOW TRIALS 
%%%%%%%%%%%%%%%%%%%%%%%

xa = 67*1000;
xb = 69*1000;

fig_ts = figure;
[ha] = tight_subplot(3,1,[.04 .04],[.05 .05],[.1 .1]);

axes(ha(1));
plot(Xm-smooth(Xm,100),'color','k'); hold on
plot(filter_m); xlim([xa xb]);
axes(ha(2));
plot(Xr-smooth(Xr,100),'color',[0, 102, 0]/255); hold on
plot(filter_r,'color',[0, 179, 0]/255);
xlim([xa xb]);
axes(ha(3));
plot(Xs-smooth(Xs,100),'color',[102, 0, 102]/255); hold on
plot(filter_s,'color',[204, 0, 204]/255); 
xlim([xa xb]);

set(gcf, 'Position',  [100, 600, 700, 400]);

dir_ts_rec_fig = strcat(dir_high_low_theta,'/time_series');
fname = strcat(dir_ts_rec_fig,sprintf('/time_series_MSR_high_power.fig',cnt_m));
saveas(fig_ts,fname);
fname = strcat(dir_ts_rec_fig,sprintf('/time_series_MSR_high_power.pdf',cnt_m));
saveas(fig_ts,fname);




low = sort(low_idx)'
high = sort(high_idx)'

end 

% M1  = find(strcmp(mod_rec_stim.RecordPairMRIlabels(:,1),'M1'));
% CN  = find(strcmp(mod_rec_stim.RecordPairMRIlabels(:,1),'CN'));

