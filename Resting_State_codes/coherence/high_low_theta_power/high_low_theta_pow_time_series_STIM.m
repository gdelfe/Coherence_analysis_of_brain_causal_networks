
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the time series (LFP) for low and high theta power of the
% modulator for the STIM experiment. High theta power time series should show oscillatory rhythms
% at the theta frequency.
%
%    @ Gino Del Ferraro, June 2022, Pesaran lab, NYU
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

in_name = 'mod_rec_stim.mat'
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

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
s = 16;
Sess = sess_info{1}(s); % Session number

% Load STIM data for that session 
dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d/',Sess));
load(strcat(dir_Stim_Sess,in_name));

% Print modulators brain regions  
unique(mod_rec_stim.mod_areas)

% select modulator:
cnt_m = 14;
lfp_m = mod_rec_stim.lfp_E(:,cnt_m,:);

% Select high and low theta power trials 
W = 3;
[spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200);

theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
theta_pow_mean = mean(theta_pow); % get the average theta power
theta_pow = theta_pow - theta_pow_mean; % rescale the theta power by removing the mean value

[sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution

% low and high theta power indexes 
low_theta_pow = sort_theta(1:cut);
low_idx = trial_idx(1:cut);

high_theta = sort_theta(end-cut+1:end);
high_idx = trial_idx(end-cut+1:end);



fig_ts = figure;
[ha, pos] = tight_subplot(10,1,[.01 .03],[.2 .01],[.1 .1])

% plot LFP for a few trials 
for i = 0:4
    
    axes(ha(i+1))
    plot(mod_rec.lfp_E_high(end-i,:),'color','b')
    %set(gca,'xticklabel',[]);
    grid on
    
    axes(ha(10-i))
    plot(mod_rec.lfp_E_low(i+1,:),'color','r')
    grid on
    
    % set(ha(1:5),'XTickLabel',''); set(ha,'YTickLabel','')
    
    ylabel('Lfp','FontName','Arial','FontSize',12);
    set(gcf, 'Position',  [100, 600, 1500, 1000]);
    
    
    
    
    
    % save the plots
    dir_ts_rec_fig = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Figures/Time_Series',Sess));
    if ~exist(dir_ts_rec_fig, 'dir')
        mkdir(dir_ts_rec_fig)
    end
    
    fname = strcat(dir_ts_rec_fig,sprintf('/modulator_%d_lfp_high_low_theta.jpg',cnt_m));
    saveas(fig_ts,fname);
    
    
    
end




