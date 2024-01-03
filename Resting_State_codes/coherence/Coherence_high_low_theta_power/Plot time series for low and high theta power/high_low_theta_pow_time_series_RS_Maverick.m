
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

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_both = strcat(dir_main,sprintf('/both_monkeys/spectrum/%s/',monkey));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% area_tot = {};
% for s = 1:length(sess_info{1})
s = 1 % 1:size(sess_info{1},1)
Sess = sess_info{1}(s); % Session number
cnt_m = 1;

Sess = sess_info{1}(s); % Session number
dir_Sess_mod_send_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_send/Data',Sess));
dir_Sess_mod_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Data',Sess));

dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));

load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed

sess_data_lfp
areas = sess_data_lfp.RecordPairMRIlabels(:,1);
areas =  areas(~cellfun('isempty',areas));
% display(['Session ------- ',num2str(s)])
area_sess = unique(areas);
% area_tot = [area_tot; area_sess];
%     
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
    areas = sess_data_lfp.RecordPairMRIlabels(:,1);
    areas =  areas(~cellfun('isempty',areas));
    unique(areas)

    M1  = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'M1'));
    CN  = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'CN'));
    OFC = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'OFC'));
    ACC = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'ACC'));
    dPFC = find(strcmp(sess_data_lfp.RecordPairMRIlabels(:,1),'dPFC'));

    
    lfp_m = sq(sess_data_lfp.lfp_E(el,:,:));
    
    
    % 1 figure
    figure;
    lfp_el = sq(sess_data_lfp.lfp_E(OFC(1),:,:));
    lfp_T = lfp_el';
    X = lfp_T(:)';
    filter = thetafilter(X,1000);
    filter2 = thetafilter(X,800);
    filter3 = thetafilter(X,1200);
    
    plot(X,'color','b'); hold on
%     plot(X-smooth(X,100)); hold on
%     plot(smooth(X,100)); hold on
    plot(filter,'color','r'); hold on;
    plot(filter2,'color','k'); hold on;
    plot(filter3,'color','g'); hold on 
%     plot(X-smooth(X,100)); hold on 
%     plot(X-smooth(X,200)); hold on
    legend
    xlim([34200 35800])
    
    %set(gca,'xticklabel',[]);
    grid on
    ylabel('Lfp','FontName','Arial','FontSize',12);
    set(gcf, 'Position',  [100, 600, 2000, 500]);
    
        
    id = 2;
    reg = [M1(1),CN(2),OFC(1),ACC(5),dPFC(2)];
    reg_name = {'M1','CN','OFC','ACC','dPFC'};
    start = 1;
    stop = 20;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALL PLOTS IN ONE FIGURE, BUT SEPARATED
    reg = [M1(1),dPFC(4),ACC(3), CN(5),OFC(3)];
    reg_name = {'M1','dPFC','ACC','CN','OFC'};


    fig_ts = figure;
    ha = tight_subplot(5,1,[.04 .04],[.2 .01],[.1 .1])
    for i = 1:5
        
        % chose the electrode channel and concatenate the splitted trials
        lfp_el = sq(sess_data_lfp.lfp_E(reg(i),:,:));
        lfp_T = lfp_el';
        X = lfp_T(:)';
        X = lfpfilter_low(X,1000);

        filter3 = thetafilter(X,800);
        
        axes(ha(i))
%         plot(lfp_ns(start:1000*stop),'color','b')
        plot(X-smooth(X,100),'color','b'); hold on
        plot(filter3,'color','r'); hold on 
        ylabel(sprintf('%s',reg_name{i}),'FontName','Arial','FontSize',12);

        %set(gca,'xticklabel',[]);
        grid on
        xlim([38500 41500])

    end
    % set(ha(1:5),'XTickLabel',''); set(ha,'YTickLabel','')
    
    set(gcf, 'Position',  [100, 600, 800, 600]);
      
    
    % MAVERICK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ONE FIGURE, ALL PLOTS TOGETHER SHIFT ON THE Y-AXIS
    
%     reg = [M1(4),48,CN(2),OFC(4),ACC(2),dPFC(3)];
    reg = [M1(1), CN(5),OFC(3),ACC(3),dPFC(4)];
    fig_ts = figure
    
%     lfp_S = sq(sess_data_lfp.lfp_S);
%     lfp_T = lfp_S';
%     X = lfp_T(:)';
%     
%     filter3 = thetafilter(X,800);
%     
%     plot(X-smooth(X,100),'color','b'); hold on
%     plot(filter3,'color','r'); hold on
%     grid on

    for i = 1:4
        
        % chose the electrode channel and concatenate the splitted trials
        lfp_el = sq(sess_data_lfp.lfp_E(reg(i),:,:));
        lfp_T = lfp_el';
        X = lfp_T(:)';
        X = lfpfilter_low(X,1000);
        
        
        filter3 = thetafilter(X,800);

        plot(X-smooth(X,100)+80*i,'color','b'); hold on 
        plot(filter3+80*i,'color','r'); hold on
        grid on 

    end
        xlim([38500 42000])
        set(gcf, 'Position',  [100,100, 1500, 800]);


    
    
    
    %         set(gca,'xticklabel',[]);
    grid on
    ylabel(sprintf('%s',reg_name{i}),'FontName','Arial','FontSize',12);
    set(gcf, 'Position',  [100, -400+i*300, 2000, 150]);
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MULTIPLE FIGURE, EACH PLOT IN A SEPARATE FIGURE 
    
    reg = [M1(1),CN(2),OFC(1),ACC(5),dPFC(5)];
    
    for i = 1:5
        
        % chose the electrode channel and concatenate the splitted trials
        lfp_el = sq(sess_data_lfp.lfp_E(reg(i),:,:));
        lfp_T = lfp_el';
        X = lfp_T(:)';
        X = lfpfilter_low(X,1000);
        
        
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
        
%         fname = strcat(dir_RS_Theta,sprintf('/time_series/lfp_reg_%s_sess_%d.pdf',reg_name{i},Sess))
%         saveas(fig_ts,fname);
%         fname = strcat(dir_RS_Theta,sprintf('/time_series/lfp_reg_%s_sess_%d.fig',reg_name{i},Sess))
%         saveas(fig_ts,fname);
    end
    
    
    
    reg = [M1(1),CN(2),OFC(1),ACC(5),dPFC(5)];
    reg_name = {'M1','CN','OFC','ACC','dPFC'};

    fk = [0 50];
    tapers = [0.7 3];
    dn = 0.005;
    fs = 1000;
    pad = 2;
    CLim = [4,11]
    
    for i = 1:5
        
        lfp_el = sq(sess_data_lfp.lfp_E(reg(i),:,:));
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
    
    
    % %%%%%%%%%%%%
    % SPECTRUM
    %%%%%%%%%%%%%%
    W = 8;
    for i = 1:5
        
        lfp_el = sq(sess_data_lfp.lfp_E(reg(i),:,:));
        % SPECTRUM
        [spec, f, err] = dmtspec(lfp_el, [1000/1e3,W],1e3,200,2,0.05,1);
        
        fig = figure;
        plot(f,log(spec));
        
        title(sprintf('%s',reg_name{i}))
        
        ylabel('log(PSD)','FontName','Arial','FontSize',10);
        xlabel('time','FontName','Arial','FontSize',10);
        xlim([0 95])
        grid on
        set(gca,'FontSize',13)
        %         ax = gca;
        %         ax.CLim = [4,11];
        
        %         set(gcf, 'Position',  [100,100, 400, 300]);
        %         set(gcf, 'Position',  [100, -400+i*300, 2000, 150]);
        
        fname = strcat(dir_both,sprintf('spectrum_reg_%s_sess_%d.pdf',reg_name{i},Sess))
        saveas(fig,fname);
        fname = strcat(dir_both,sprintf('spectrum_reg_%s_sess_%d.fig',reg_name{i},Sess))
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




