
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
s = 15;
Sess = sess_info{1}(s); % Session number

% Load STIM data for that session
dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d/',Sess));
load(strcat(dir_Stim_Sess,in_name));

% Print modulators brain regions
unique(mod_rec_stim.mod_areas)

% select modulator:
% cnt_m = 1; % modulator index
for cnt_m = 1 %1 : length(mod_rec_stim.mod_idx)
    
    % % Modulators
%     el = mod_rec_stim.mod_idx(cnt_m); % electrode
%     lfp_m = sq(mod_rec_stim.lfp_E(:,el,:));
    
    % General electrodes 
    el = 5;
    lfp_m = sq(mod_rec_stim.lfp_E(:,el,:));
    
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
    
    % Power trials
    lfp_E_high = lfp_m(high_idx,:); % high theta
    lfp_E_low = lfp_m(low_idx,:); % low theta
    
    % plot theta power histogram
    % fig = figure;
    % histogram(high_theta,30,'FaceAlpha',0.5)
    % hold on
    % histogram(low_theta,30,'FaceAlpha',0.5)
    
    fig = figure;
    histogram(theta_pow,40,'FaceAlpha',0.5,'Normalization','probability','EdgeColor','w')
    grid on
    hold on; xline(low_theta(end),'b');
    hold on; xline(high_theta(1),'r');
    legend('theta power','low theta','high theta')
    ylim([0 0.1])
    dir_ts_rec_fig = strcat(dir_Stim_Theta,sprintf('/high_low_theta',Sess));
    
    
    
    fname = strcat(dir_ts_rec_fig,'/modulator_theta_histogram.pdf');
    saveas(fig,fname);
    fname = strcat(dir_ts_rec_fig,'/modulator_theta_histogram.jpg');
    saveas(fig,fname);
    
    
    % plot LFP for a few trials
    
    for start = [0]
        fig_ts = figure;
        [ha] = tight_subplot(10,1,[.02 .02],[.05 .05],[.1 .1]);
        j = 0;
        for i = start:start+4
            
            axes(ha(j+1));
            plot(lfp_E_high(end-i,:),'color','b') % plot trials on the right tail of the distribution
            %set(gca,'xticklabel',[]);
            ylim([-330,330])
            grid on
            
            axes(ha(end-j));
            plot(lfp_E_low(i+1,:),'color','r') % plot trials on the left tail of the distribution
            grid on
            ylim([-330,330])
            
            j = j + 1;
            % set(ha(1:5),'XTickLabel',''); set(ha,'YTickLabel','')
            
            %     ylabel('Lfp','FontName','Arial','FontSize',12);
            set(gcf, 'Position',  [100, 600, 700, 1000]);
            
        end
        
        % %     save the plots
        dir_ts_rec_fig = strcat(dir_Stim_Theta,sprintf('/high_low_theta',Sess));
            
        reg = mod_rec_stim.RecordPairMRIlabels{el,1};
        %             reg = mod_rec_stim.mod_areas{cnt_m};
        fname = strcat(dir_ts_rec_fig,sprintf('/ctrl_OA_time_series_el_%d_reg_%s_sess_%d.pdf',el,reg,Sess));
        saveas(fig_ts,fname);
        fname = strcat(dir_ts_rec_fig,sprintf('/ctrl_OA_time_series_el_%d_reg_%s_sess_%d.jpg',el,reg,Sess))
        saveas(fig_ts,fname);
        
    end
    
end

M1  = find(strcmp(mod_rec_stim.RecordPairMRIlabels(:,1),'M1'));
CN  = find(strcmp(mod_rec_stim.RecordPairMRIlabels(:,1),'CN'));
OFC = find(strcmp(mod_rec_stim.RecordPairMRIlabels(:,1),'OFC'));
ACC = find(strcmp(mod_rec_stim.RecordPairMRIlabels(:,1),'ACC'));

OFC'
CN'
M1'

