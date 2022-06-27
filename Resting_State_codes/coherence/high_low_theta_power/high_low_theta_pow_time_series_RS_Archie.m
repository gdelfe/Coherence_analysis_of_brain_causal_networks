
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the time series (LFP) for low and high theta power of the
% modulator. High theta power time series should show oscillatory rhythms
% at the theta frequency. 
%
%    @ Gino Del Ferraro, November 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
dir_high_low_theta = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';

name_struct_input = '/session_data_lfp_movie.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_both = strcat(dir_main,sprintf('/both_monkeys/spectrograms/%s/',monkey));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

area_tot = {};
% for s = 1:length(sess_info{1})
s = 4; % 1:size(sess_info{1},1)
Sess = sess_info{1}(s); % Session number
cnt_m = 1;

Sess = sess_info{1}(s); % Session number
dir_Sess_mod_send_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_send/Data',Sess));
dir_Sess_mod_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Data',Sess));

dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));

load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed

% sess_data_lfp
% areas = sess_data_lfp.RecordPairMRIlabels(:,1);
% areas =  areas(~cellfun('isempty',areas));
% display(['Session ------- ',num2str(s)])
% area_sess = unique(areas)
% area_tot = [area_tot; area_sess];
    
% end

% unique(area_tot)


% for all the modulators in that session 
for cnt_m = 1 % 1:length(sess_data_lfp.mod_idx)
    
    cnt_m = 1;
    
%     
%     % General electrodes 
%     el = 5;
%     lfp_m = sq(sess_data_lfp.lfp_E(el,:,:));
    
    % Modulator 
    el = sess_data_lfp.mod_idx(cnt_m);
    lfp_m = sq(sess_data_lfp.lfp_E(el,:,:));
    
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
    
    sess_data_lfp
    
    sess_data_lfp
    areas = sess_data_lfp.RecordPairMRIlabels(:,1);
    areas =  areas(~cellfun('isempty',areas));
    unique(areas)
    
    Amyg  = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'Amyg'));
    GP  = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'GP'));
    NAc  = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'NAu'));
    Put = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'Put'));
    S1 = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'S1'));
    SMG = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'SMG'));
    SPL = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'SPL'));
    pCun = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'pCun'));
    
%     BFB = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'BFB'));
    
    
    lfp_m = sq(sess_data_lfp.lfp_E(el,:,:));
    
   
    id = 2;
    reg = {GP(1),NAc(1),Put(1),S1(1),pCun(1)};
    reg_name = {'GP','NAc','Put','S1','pCun'};


    reg = {Amyg(1),GP(1),NAc(1),Put(1),S1(1),SMG(1),SPL(1),pCun(1)};
    reg_name = {'Amyg','GP','NAc','Put','S1','SMG','SPL','pCun'};
%     reg = {BFB(1)}
%     reg_name = {'BFB'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALL THE PLOTS TOGETHER 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig_ts = figure;
    ha = tight_subplot(8,1,[.02 .02],[.2 .01],[.1 .1])
    for i = 1:8
        
        % chose the electrode channel and concatenate the splitted trials
        lfp_el = sq(sess_data_lfp.lfp_E(reg{i},:,:));
        lfp_T = lfp_el';
        X = lfp_T(:)';
        
        filter3 = thetafilter(X,600);
        
        axes(ha(i))
%         plot(lfp_ns(start:1000*stop),'color','b')
        plot(X-smooth(X,200),'color','b'); hold on 
        plot(filter3,'color','r'); hold on 
        ylabel(sprintf('%s',reg_name{i}),'FontName','Arial','FontSize',12);

        %set(gca,'xticklabel',[]);
        grid on

    end
    % set(ha(1:5),'XTickLabel',''); set(ha,'YTickLabel','')
    
    set(gcf, 'Position',  [100, 600, 2000, 1000]);
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plots separated 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:8
        
        % chose the electrode channel and concatenate the splitted trials
        lfp_el = sq(sess_data_lfp.lfp_E(reg{i},:,:));
        lfp_T = lfp_el';
        X = lfp_T(:)';
        
        filter3 = thetafilter(X,800);
        
        %         plot(lfp_ns(start:1000*stop),'color','b')
        fig_ts = figure
%         plot(X,'color','b'); hold on
        plot(X-smooth(X,200),'color','b'); hold on 
%         plot(smooth(X,5)-smooth(X,200),'color','k'); hold on 
        plot(filter3,'color','r'); hold on
        xlim([10000 13000])
%         set(gca,'xticklabel',[]);
        grid on
        ylabel(sprintf('%s',reg_name{i}),'FontName','Arial','FontSize',12);
        set(gcf, 'Position',  [100,100, 800, 150]);
%         set(gcf, 'Position',  [100, -400+i*300, 2000, 150]);
        
        fname = strcat(dir_RS_Theta,sprintf('/time_series/lfp_reg_%s_sess_%d.pdf',reg_name{i},Sess))
        saveas(fig_ts,fname);
        fname = strcat(dir_RS_Theta,sprintf('/time_series/lfp_reg_%s_sess_%d.fig',reg_name{i},Sess))
        saveas(fig_ts,fname);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Smooth signal, low pass filter, Plots separated 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for i = 1%:5
        
        % chose the electrode channel and concatenate the splitted trials
        lfp_el = sq(sess_data_lfp.lfp_E(reg{i},:,:));
        lfp_T = lfp_el';
        X = lfp_T(:)';
        X = lfpfilter_low(X,1000);
        
        
        filter3 = thetafilter(X,800);
        
        fig_ts = figure
        plot(X-smooth(X,200),'color','b'); hold on
        plot(filter3,'color','r'); hold on
        xlim([10000 13000])
        grid on
        ylabel(sprintf('%s',reg_name{i}),'FontName','Arial','FontSize',12);
        set(gcf, 'Position',  [100,100, 800, 150]);
        %         set(gcf, 'Position',  [100, -400+i*300, 2000, 150]);
        
%         fname = strcat(dir_RS_Theta,sprintf('/time_series/lfp_reg_%s_sess_%d.pdf',reg_name{i},Sess))
%         saveas(fig_ts,fname);
%         fname = strcat(dir_RS_Theta,sprintf('/time_series/lfp_reg_%s_sess_%d.fig',reg_name{i},Sess))
%         saveas(fig_ts,fname);
    end
    
    
    fk = [0 50];
    tapers = [0.7 3];
    dn = 0.005;
    fs = 1000;
    pad = 2;
    
    for i = 1%:8
        
        lfp_el = sq(sess_data_lfp.lfp_E(reg{i},:,:));
        % SPECTROGRAM
        [specRS, fRS , tiRS] = tfspec(lfp_el,tapers,fs,dn,fk,pad,0.05,1,0);
        
        fig = figure; tvimage(sq(log(specRS)));
        title(sprintf('%s',reg_name{i}))
        ticks = 0:10:100;
        ticklabels = 0:5:50;
        yticks(ticks)
        yticklabels(ticklabels)
        colorbar
        ylabel('frequency','FontName','Arial','FontSize',10);
        xlabel('time','FontName','Arial','FontSize',10);
        %         ax = gca;
        %         ax.CLim = [4,11];
        
        %         set(gcf, 'Position',  [100,100, 400, 300]);
        %         set(gcf, 'Position',  [100, -400+i*300, 2000, 150]);
        
        fname = strcat(dir_both,sprintf('spec_reg_%s_sess_%d.pdf',reg_name{i},Sess))
        saveas(fig,fname);
        fname = strcat(dir_both,sprintf('spec_reg_%s_sess_%d.fig',reg_name{i},Sess))
        saveas(fig,fname);
    end
    
    


    
end


% X = lfp_ns;
% filter3 = thetafilter(X,1000);
% figure;
% plot(X); hold on
% plot(filter)
% grid on
% figure
% plot(smooth(X,100))




