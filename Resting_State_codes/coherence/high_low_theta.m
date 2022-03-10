
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code ..
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

coh_all_c_ms_high = [];
coh_all_c_ms_low = [];
coh_all_c_mr_high = [];
coh_all_c_mr_low = [];

for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    clear mod_rec mod_send coh_all sess_data_lfp
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info 

    %     % ---  time parameter
    tot_time = 150001;
    
    % outliers time series in sender and receiver 
    outliers_S = sess_data_lfp.outliers_S;
    outliers_R = sess_data_lfp.outliers_R;
    

    % number of modulators in that session 
    n_mod = size(sess_data_lfp.mod_idx,2);
    cnt_m = 1; % counter for numb of modulators checked for the outlier trials remover
    
    % matrices to store high/low coherence wrt sender and receiver 
    coh_all_mod_ms_high = zeros(n_mod,409);
    coh_all_mod_ms_low = zeros(n_mod,409);
    coh_all_mod_mr_high = zeros(n_mod,409);
    coh_all_mod_mr_low = zeros(n_mod,409);
    
    
    for m = sess_data_lfp.mod_idx
        
        display(['-- Modulator ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);
        
        
        % When computing MS and MR coherence we have to remove the outliers
        % for the sender/receiver, and this are different time series for
        % the sender and the receiver 
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%
        % In relation to the sender
        
        lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
        lfp_S = sess_data_lfp.lfp_S;
        
        outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;
        
        outliers_ES = [outliers_S, outliers_E];
        outliers_ES = unique(outliers_ES);  % -- remove repeated entries in outliers
        
        lfp_E(outliers_ES,:) = [];
        lfp_S(outliers_ES,:) = [];
        
        % Compute the spectrum for each trial. Format: iTrial x times
        W = 3;
        [spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200);
        
        % Find low and high theta from the spectrum 
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
        
        mod_send.mod_idx = m;
        mod_send.all_theta_pow = theta_pow;
        mod_send.high_theta_pow = high_theta;
        mod_send.high_pow_idx = high_idx;
        mod_send.low_theta_pow = low_theta_pow;
        mod_send.low_pow_idx = low_idx;
        mod_send.all_trial = lfp_E;
        mod_send.lfp_E_high = lfp_E(high_idx,:);
        mod_send.lfp_E_low = lfp_E(low_idx,:);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Coherency Modulator-Sender        %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ---- parameters for the coherence-gram
        nt = tot_time;
        fs = 1000;
        fk = 200;
        pad = 2;
        N = 1;
        W = 5;
        
        % -- coherence calculation via coherency()
        [c_ms,f] = coherency(lfp_S,lfp_E,[N W],fs,fk,pad,0.05,1,1);
        [c_ms_high,f] = coherency(lfp_S(high_idx,:),lfp_E(high_idx,:),[N W],fs,fk,pad,0.05,1,1);
        [c_ms_low,f] = coherency(lfp_S(low_idx,:),lfp_E(low_idx,:),[N W],fs,fk,pad,0.05,1,1);
        
        
        mod_send.coh_ms = c_ms;
        mod_send.coh_ms_high = c_ms_high;
        mod_send.coh_ms_low = c_ms_low;
        mod_send.freq = f;
        
        
        
        fig_histo = figure;
        histogram(low_theta_pow,20,'FaceAlpha',.6); grid on
        hold on
        histogram(high_theta,20,'FaceAlpha',.6);
        legend('low theta','high theta')
        title('Theta power distribution Mod S','FontSize',12)
        ylabel('count')
        xlabel('Log Theta power (mean centered)')
        
        
        fig_coh = figure;
        plot(f,abs(c_ms));
        hold on
        plot(f,abs(c_ms_high));
        hold on
        plot(f,abs(c_ms_low));
        title('coherence MS - trials with high and low theta pow')
        legend('all trials','high pow trials','low pow trials')
        ylabel('cohernece')
        xlabel('frequency')
        grid on
        
        
        
        dir_Sess_mod_send_fig = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_send/Figures',Sess));
        if ~exist(dir_Sess_mod_send_fig, 'dir')
            mkdir(dir_Sess_mod_send_fig)
        end
        
        dir_Sess_mod_send_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_send/Data',Sess));
        if ~exist(dir_Sess_mod_send_data, 'dir')
            mkdir(dir_Sess_mod_send_data)
        end
        
        
        fname = strcat(dir_Sess_mod_send_fig,sprintf('/modulator_%d_theta_pow_histo_S.jpg',cnt_m));
        saveas(fig_histo,fname);
        fname = strcat(dir_Sess_mod_send_fig,sprintf('/modulator_%d_sender_coherence.jpg',cnt_m));
        saveas(fig_coh,fname);
        
        save(strcat(dir_Sess_mod_send_data,sprintf('/modulator_%d_send',cnt_m)),'mod_send');
        
        
        
        
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%
        % In relation to the receiver
        
        lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
        lfp_R = sess_data_lfp.lfp_R;
        
        outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;
        
        outliers_ER = [outliers_R, outliers_E];
        outliers_ER = unique(outliers_ER);  % -- remove repeated entries in outliers
        
        lfp_E(outliers_ER,:) = [];
        lfp_R(outliers_ER,:) = [];
        
        
        % Compute the spectrum for each trial. Format: iTrial x times
        W = 3;
        [spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200);
        theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
        theta_pow_mean = mean(theta_pow);
        theta_pow = theta_pow - theta_pow_mean;
        
        [sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
        cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution
        
        low_theta_pow = sort_theta(1:cut);
        low_idx = trial_idx(1:cut);
        
        high_theta = sort_theta(end-cut+1:end);
        high_idx = trial_idx(end-cut+1:end);
        
        mod_rec.mod_idx = m;
        mod_rec.all_theta_pow = theta_pow;
        mod_rec.high_theta_pow = high_theta;
        mod_rec.high_pow_idx = high_idx;
        mod_rec.low_theta_pow = low_theta_pow;
        mod_rec.low_pow_idx = low_idx;
        mod_rec.all_trial = lfp_E;
        mod_rec.lfp_E_high = lfp_E(high_idx,:);
        mod_rec.lfp_E_low = lfp_E(low_idx,:);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Coherency Modulator-Receiver        %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ---- parameters for the coherence-gram
        nt = tot_time;
        fs = 1000;
        fk = 200;
        pad = 2;
        N = 1;
        W = 5;
        
        % -- coherence calculation via coherency()
        [c_mr,f] = coherency(lfp_R,lfp_E,[N W],fs,fk,pad,0.05,1,1);
        [c_mr_high,f] = coherency(lfp_R(high_idx,:),lfp_E(high_idx,:),[N W],fs,fk,pad,0.05,1,1);
        [c_mr_low,f] = coherency(lfp_R(low_idx,:),lfp_E(low_idx,:),[N W],fs,fk,pad,0.05,1,1);
        
        
        mod_rec.coh_mr = c_mr;
        mod_rec.coh_mr_high = c_mr_high;
        mod_rec.coh_mr_low = c_mr_low;
        mod_rec.freq = f;
        
        
        
        fig_histo = figure;
        histogram(low_theta_pow,20,'FaceAlpha',.6); grid on
        hold on
        histogram(high_theta,20,'FaceAlpha',.6);
        legend('low theta','high theta')
        title('Theta power distribution Mod R','FontSize',12)
        ylabel('count')
        xlabel('Log Theta power (mean centered)')
        
        
        fig_coh = figure;
        plot(f,abs(c_ms));
        hold on
        plot(f,abs(c_ms_high));
        hold on
        plot(f,abs(c_ms_low));
        title('coherence MR - trials with high and low theta pow')
        legend('all trials','high pow trials','low pow trials')
        ylabel('cohernece')
        xlabel('frequency')
        grid on
        
        
        dir_Sess_mod_rec_fig = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Figures',Sess));
        if ~exist(dir_Sess_mod_rec_fig, 'dir')
            mkdir(dir_Sess_mod_rec_fig)
        end
        
        dir_Sess_mod_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Data',Sess));
        if ~exist(dir_Sess_mod_rec_data, 'dir')
            mkdir(dir_Sess_mod_rec_data)
        end
        
        
        fname = strcat(dir_Sess_mod_rec_fig,sprintf('/modulator_%d_theta_pow_histo_R.jpg',cnt_m));
        saveas(fig_histo,fname);
        fname = strcat(dir_Sess_mod_rec_fig,sprintf('/modulator_%d_receiver_coherence.jpg',cnt_m));
        saveas(fig_coh,fname);
        
        save(strcat(dir_Sess_mod_rec_data,sprintf('/modulator_%d_rec',cnt_m)),'mod_rec');
        
        
        coh_all_mod_ms_high(cnt_m,:) = mod_send.coh_ms_high;
        coh_all_mod_ms_low(cnt_m,:) = mod_send.coh_ms_low;
        
        coh_all_mod_mr_high(cnt_m,:) = mod_rec.coh_mr_high;
        coh_all_mod_mr_low(cnt_m,:) = mod_rec.coh_mr_low;
        
        cnt_m = cnt_m +1;
        
    end % end of for cycle for all the modulators in a given session 
    
    
    
    mean_coh_ms_high = mean(abs(coh_all_mod_ms_high),1);
    mean_coh_ms_low = mean(abs(coh_all_mod_ms_low),1);
    mean_coh_mr_high = mean(abs(coh_all_mod_mr_high),1);
    mean_coh_mr_low = mean(abs(coh_all_mod_mr_low),1);
    
    err_coh_ms_high = std(abs(coh_all_mod_ms_high),0,1)/sqrt(n_mod);
    err_coh_ms_low = std(abs(coh_all_mod_ms_low),0,1)/sqrt(n_mod);
    err_coh_mr_high = std(abs(coh_all_mod_mr_high),0,1)/sqrt(n_mod);
    err_coh_mr_low = std(abs(coh_all_mod_mr_low),0,1)/sqrt(n_mod);
    
    
    % Coherence of all the modulators for a given session 
    coh_all_mod.c_ms_high = coh_all_mod_ms_high;
    coh_all_mod.c_ms_low = coh_all_mod_ms_low;
    coh_all_mod.c_mr_high = coh_all_mod_mr_high;
    coh_all_mod.c_mr_low = coh_all_mod_mr_low;
    
    
    % Store coherence for all the modulators 
    coh_all_c_ms_high = [coh_all_c_ms_high; coh_all_mod_ms_high];
    coh_all_c_ms_low = [coh_all_c_ms_low; coh_all_mod_ms_low];
    coh_all_c_mr_high = [coh_all_c_mr_high; coh_all_mod_mr_high];
    coh_all_c_mr_low = [coh_all_c_mr_low; coh_all_mod_mr_low];
   
    
    save(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat',Sess)),'coh_all_mod');
    
    
    set(0,'DefaultFigureVisible','off')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));
    
    
    
    fig = figure;
    shadedErrorBar(f,mean_coh_ms_high,err_coh_ms_high,'lineprops',{'color',[28 199 139]/255 },'patchSaturation',0.5); hold on
    shadedErrorBar(f,mean_coh_ms_low,err_coh_ms_low,'lineprops',{'color',[50 250 93]/255 },'patchSaturation',0.5);
    
    grid on
    % title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
    set(gca,'FontSize',14)
    xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
    ylabel('Coherence','FontName','Arial','FontSize',15);
    title('Coherence MS high vs low theta power trial','FontSize',12)
    legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
    set(gcf, 'Position',  [100, 600, 898, 500])
    xlim([1 95])
    % ylim([0 0.25])
    grid on
    
    
    fname = strcat(dir_Sess_mod_rec_fig,sprintf('/MS_coherence_mod_mean.jpg',cnt_m));
    saveas(fig,fname);
    
    
    
    fig = figure;
    hold all
    shadedErrorBar(f,mean_coh_mr_high,err_coh_mr_high,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.5);
    shadedErrorBar(f,mean_coh_mr_low,err_coh_mr_low,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.5);
    
    grid on
    % title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
    set(gca,'FontSize',14)
    xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
    ylabel('Coherence','FontName','Arial','FontSize',15);
    title('Coherence MS high vs low theta power trial','FontSize',12)
    legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
    set(gcf, 'Position',  [100, 600, 898, 500])
    xlim([1 95])
    % ylim([0 0.25])
    grid on
    
    fname = strcat(dir_Sess_mod_rec_fig,sprintf('/MR_coherence_mod_mean.jpg',cnt_m));
    saveas(fig,fname);
    
   
    
end



save(strcat(dir_high_low_theta,'/coh_all_sess_ms_high.mat'),'coh_all_c_ms_high')
save(strcat(dir_high_low_theta,'/coh_all_sess_ms_low.mat'),'coh_all_c_ms_low')
save(strcat(dir_high_low_theta,'/coh_all_sess_mr_high.mat'),'coh_all_c_mr_high')
save(strcat(dir_high_low_theta,'/coh_all_sess_mr_low.mat'),'coh_all_c_mr_low')


% load(strcat(dir_high_low_theta,'/coh_all_sess_ms_high.mat'))
% load(strcat(dir_high_low_theta,'/coh_all_sess_ms_low.mat'))
% load(strcat(dir_high_low_theta,'/coh_all_sess_mr_high.mat'))
% load(strcat(dir_high_low_theta,'/coh_all_sess_mr_low.mat'))


mean_all_coh_ms_high = mean(abs(coh_all_c_ms_high),1);
mean_all_coh_ms_low = mean(abs(coh_all_c_ms_low),1);
mean_all_coh_mr_high = mean(abs(coh_all_c_mr_high),1);
mean_all_coh_mr_low = mean(abs(coh_all_c_mr_low),1);

err_all_coh_ms_high = std(abs(coh_all_c_ms_high),0,1)/sqrt(size(coh_all_c_ms_high,1));
err_all_coh_ms_low = std(abs(coh_all_c_ms_low),0,1)/sqrt(size(coh_all_c_ms_low,1));
err_all_coh_mr_high = std(abs(coh_all_c_mr_high),0,1)/sqrt(size(coh_all_c_mr_high,1));
err_all_coh_mr_low = std(abs(coh_all_c_mr_low),0,1)/sqrt(size(coh_all_c_mr_low,1));


% zscore_ms_high = (abs(coh_all_c_ms_high) - mean_all_coh_ms_high)./std(abs(coh_all_c_ms_high),0,1);
% zscore_ms_low = (abs(coh_all_c_ms_low) - mean_all_coh_ms_low)./std(abs(coh_all_c_ms_low),0,1);
% zscore_mr_high = (abs(coh_all_c_mr_high) - mean_all_coh_mr_high)./std(abs(coh_all_c_mr_high),0,1);
% zscore_mr_low = (abs(coh_all_c_mr_low) - mean_all_coh_mr_low)./std(abs(coh_all_c_mr_low),0,1);


% mean_zscore_ms_high = mean(zscore_ms_high,1);
% figure;
% plot(f,mean_all_coh_ms_high)
% hold on
% plot(f,abs(coh_all_c_ms_high))
% 


set(0,'DefaultFigureVisible','on')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));

fig = figure;
shadedErrorBar(f,mean_all_coh_mr_high,err_all_coh_mr_high,'lineprops',{'color',[28 199 139]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_mr_low,err_all_coh_mr_low,'lineprops',{'color',[50 250 93]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title('Coherence MR high vs low theta power trial','FontSize',12)
legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,sprintf('/MR_all_coherence_mean.jpg',cnt_m));
saveas(fig,fname);

