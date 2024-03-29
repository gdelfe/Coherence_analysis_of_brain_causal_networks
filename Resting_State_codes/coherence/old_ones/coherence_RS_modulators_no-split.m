
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the causal 
% modulators found by Shaoyu's and both the sender and the receiver
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators
%
% INPUT: file with session modulator info
%        .mat file with structure AM and MA information 
%
% OUTPUT: txt files with the values of the coherence MR, MS, SR and
% corresponding figures
%
%    @ Gino Del Ferraro, August 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_base = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
step = 110;

fid = fopen(strcat(dir_base,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

% -- load structure files
load(strcat(dir_base,'/session_AM.mat'))
load(strcat(dir_base,'/session_MA.mat'))


% -- print structures on stdout
%format short
for s=1:size(sess_info{1},1)
    session_AM(s)
    session_MA(s)
end

% -- matrices to store coherence for each mod and then compute
% average
ftot = 1638; % number of frequency points 
% 48 is the number of overall modulators across sections
coh_ms_AM = zeros(48,ftot); % to store abs mean coherence mod-sender. # of causal mod, # of freq bins
coh_mr_AM = zeros(48,ftot); % to store abs mean coherence mod-receiver.# of causal mod, # of freq bins

coh_ms_MA = zeros(48,ftot); % mean abs choerence # of causal mod, # of freq bins
coh_mr_MA = zeros(48,ftot); % mean abs coherence # of causal mod, # of freq bins

coh_sr_AM = zeros(20,ftot); % sender-receiver abs mean coherence # of causal mod, # of freq bin
coh_sr_MA = zeros(20,ftot); % sender-receiver mean abs coherence # of causal mod, # of freq bin

cnt = 1;
cnt_sr = 1;

for i=1:size(sess_info{1},1)  % For all the session with a modulator
    
    
    close all
    % addpath('/vol/sas8/Maverick_RecStim_vSUBNETS220/160125/004')
    addpath(sprintf('/vol/sas8/Maverick_RecStim_vSUBNETS220/%s/%s/',sess_info{2}{i},sess_info{3}{i})) % add path of the specific RS session
    
    % file = 'rec004.Frontal.lfp.dat'
    file = sprintf('rec%s.Frontal.lfp.dat',sess_info{3}{i})
    fid = fopen(file);
    format = 'float=>single';
    
    CH = 220; % tot number of channels
    FS = 1000; % sampling
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',   out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    data = fread(fid,[CH,inf],format); % load the RS data
    % h = fread(fid,[CH,diff(bn)*FS./1e3],format);
    % ---- bipolar referencing, pairs of electrodes
    dir_Sess = strcat(dir_base,sprintf('/Sess_%d',Sess));
    if ~exist(dir_Sess, 'dir')
        mkdir(dir_Sess)
    end
    
    % -- load list electrodes, sender, receiver
    electrode = importdata(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairs
    receiver = importdata(strcat(dir_Sess,sprintf('/receiver_Sess_%d.txt',Sess)));  % ---- receiver pair
    sender = importdata(strcat(dir_Sess,sprintf('/sender_Sess_%d.txt',Sess))); % ---- sender pair
    
    
    % ---  time parameter
    tot_time = 150000;
    % ---  freq parameter for the masking
    fmin = 10;
    fmax = 40;
    
    % ---- Lfp of the resting state for that specific pair of electrodes
    lfpRS = data(electrode(:,1),:) - data(electrode(:,2),:); % all potential modulators
    lfp_S = data(sender(1),:) - data(sender(2),:); % sender
    lfp_R = data(receiver(1),:) - data(receiver(2),:); % receiver
    
    % include signal up to time where signal is not corrupted
    lfpRS = lfpRS(:,1:tot_time);
    lfp_S = lfp_S(:,1:tot_time);
    lfp_R = lfp_R(:,1:tot_time);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- COHERENCE-GRAM Sender-Receiver -- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    display(['Computing coherence-gram...'])
    % ---- parameters for the coherence-gram
    tapers = [4 4];
    N = tapers(1);
    nt = tot_time;
    dn = 0.01;
    fs = 1000;
    fk = 200;
    pad = 2;
    nwin = single(floor((nt-N*fs)/(dn*fs)))
    
    % --- coherence
    
    mod_Ch = session_AM(i).mod_idx; % causal modulator channel
    
    
    % -- Coherogram sender-receiver calculation
    [c_sr,tf,f,spec_r,spec_s] = tfcoh_GINO(lfp_R,lfp_S,tapers,1e3,dn,fk,2,[],[],1); % coherence sender-receiver
   
    % --- FIGURE: SENDER-RECEIVER COHEROGRAM
    fig_sr = figure; tvimage(abs(c_sr(:,:))); colorbar; % coherence spectrum
    xticks = floor(linspace(1,length(tf),5));
    xticklabels = tf(xticks);
    xtickformat('%d')
    yticks = 1:100:length(f);
    yticklabels = floor(f(yticks));
    ytickformat('%.2f')
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
    title(sprintf('S-R Coherogram, sess = %d',Sess),'FontSize',12);
    xlabel('time (sec)');
    ylabel('freq (Hz)')
%     ylim([0,500])
    set(gcf, 'Position',  [100, 600, 1000, 600])
    
    fname = strcat(dir_Sess,sprintf('/SR_coherogram_fk_%d.jpg',fk));
    saveas(fig_sr,fname);
     
    % -- store coherence values sender-receiver 
    coh_sr_AM(cnt_sr,:) = abs(mean(c_sr(:,:),1)) ; % assign coherence S-R value 
    coh_sr_MA(cnt_sr,:) = mean(abs(c_sr(:,:)),1) ; % assign coherence S-R value 
    cnt_sr = cnt_sr + 1;
    
    for Ch = mod_Ch % for all the modulators in the session
        
        close all        
        
        % -- coherence-grams for modulator-sender, modulator-receiver 
        [c_ms,tf,f,spec_r,spec_s] = tfcoh_GINO(lfpRS(Ch,:),lfp_S,tapers,fs,dn,fk,2,[],[],1); % coherence modulator-sender
        [c_mr,tf,f,spec_r,spec_s] = tfcoh_GINO(lfpRS(Ch,:),lfp_R,tapers,fs,dn,fk,2,[],[],1); % coherence modulator-receiver
        
        
        %     dlmwrite(strcat(dir_Sess,'/sr_coherogram.txt'),c_sr,'delimiter',' ');
        % -- Figure: coherence spectrum
        
        % --- FIGURE: MODULATOR-SENDER COHEROGRAM
        fig_ms = figure; tvimage(abs(c_ms(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:100:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('M-S Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Sess,sprintf('/MS_coherogram_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig_ms,fname);
        
        % --- FIGURE: MODULATOR-RECEIVER COHEROGRAM
        fig_mr = figure; tvimage(abs(c_mr(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:100:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('M-R Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Sess,sprintf('/MR_coherogram_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig_mr,fname);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MEAN ABS COHERENCE            %%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % --- FIGURE --------- %%
        % -- Coherence vs time --- %
        fig = figure;
        plot(tf,mean(abs(c_sr(:,:)),2)); % sender-receiver
        hold on
        plot(tf,mean(abs(c_ms(:,:)),2)); % modulator-sender
        hold on
        plot(tf,mean(abs(c_mr(:,:)),2)); % modulator-receiver
        xlim([2000 50000])
        legend('S-R coherence','M-S coherence','M-R coherence')
        title(sprintf('Mean Abs coherence vs time, ch = %d, causal modulator',Ch),'FontSize',10);
        grid on
        hold off
        set(gcf, 'Position',  [100, 600, 1000, 500])
        
        fname = strcat(dir_Sess,sprintf('/coherence_vs_time_MA_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig,fname);
        
        
        % --- FIGURE --------- %%
        % -- Coherence vs frequency --- %
        fig = figure;
        plot(f,mean(abs(c_sr(:,:)),1))
        hold on
        plot(f,mean(abs(c_ms(:,:)),1))
        hold on
        plot(f,mean(abs(c_mr(:,:)),1))
        grid on
        title(sprintf('Mean Abs coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
        legend('S-R coherence','M-S coherence','M-R coherence')
%         xlim([0 60])
        set(gcf, 'Position',  [100, 600, 1000, 500])
    

        fname = strcat(dir_Sess,sprintf('/coherence_vs_freq_MA_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig,fname);
        
        coh_ms_MA(cnt,:) = mean(abs(c_ms(:,:)),1) ; % assign coherence M-S value for this modulator
        coh_mr_MA(cnt,:) = mean(abs(c_mr(:,:)),1) ; % assign coherence M-R value for this modulator
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ABS MEAN COHERENCE            %%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % --- FIGURE --------- %%
        % -- Coherence vs time --- %
        fig = figure;
        plot(tf,abs(mean(c_sr(:,:),2))); % sender-receiver
        hold on
        plot(tf,abs(mean(c_ms(:,:),2))); % modulator-sender
        hold on
        plot(tf,abs(mean(c_mr(:,:),2))); % modulator-receiver
        xlim([2000 50000])
        legend('S-R coherence','M-S coherence','M-R coherence')
        title(sprintf('Abs Mean coherence vs time, ch = %d, causal modulator',Ch),'FontSize',10);
        grid on
        hold off
        set(gcf, 'Position',  [100, 600, 1000, 500])
        
        fname = strcat(dir_Sess,sprintf('/coherence_vs_time_AM_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig,fname);
        
        
        % --- FIGURE --------- %%
        % -- Coherence vs frequency --- %
        fig = figure;
        plot(f,abs(mean(c_sr(:,:),1)))
        hold on
        plot(f,abs(mean(c_ms(:,:),1)))
        hold on
        plot(f,abs(mean(c_mr(:,:),1)))
        grid on
        title(sprintf('Abs Mean coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
        legend('S-R coherence','M-S coherence','M-R coherence')
%         xlim([0 60])
        set(gcf, 'Position',  [100, 600, 1000, 500])
        
        fname = strcat(dir_Sess,sprintf('/coherence_vs_freq_AM_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig,fname);
        
        coh_ms_AM(cnt,:) = abs(mean(c_ms(:,:),1)) ; % assign coherence M-S value for this modulator
        coh_mr_AM(cnt,:) = abs(mean(c_mr(:,:),1)) ; % assign coherence M-R value for this modulator
        
        cnt = cnt + 1;
        
    end
    
    
end

keyboard 

% Write out the data of the coherence vs frequency on txt files
dlmwrite(strcat(dir_base,sprintf('/coherence_ms_MA_fk_%d.txt',fk)),coh_ms_MA,'delimiter',' ');
dlmwrite(strcat(dir_base,sprintf('/coherence_mr_MA_fk_%d.txt',fk)),coh_mr_MA,'delimiter',' ');
dlmwrite(strcat(dir_base,sprintf('/coherence_sr_MA_fk_%d.txt',fk)),coh_sr_MA,'delimiter',' ');
dlmwrite(strcat(dir_base,sprintf('/coherence_ms_AM_fk_%d.txt',fk)),coh_ms_AM,'delimiter',' ');
dlmwrite(strcat(dir_base,sprintf('/coherence_mr_AM_fk_%d.txt',fk)),coh_mr_AM,'delimiter',' ');
dlmwrite(strcat(dir_base,sprintf('/coherence_sr_AM_fk_%d.txt',fk)),coh_sr_AM,'delimiter',' ');
dlmwrite(strcat(dir_base,sprintf('/frequency_range_fk_%d.txt',fk)),f,'delimiter',' ');

keyboard 

% --- load data 
fk = 200;
coh_ms_MA = importdata(strcat(dir_base,sprintf('/coherence_ms_MA_fk_%d.txt',fk)));
coh_mr_MA = importdata(strcat(dir_base,sprintf('/coherence_mr_MA_fk_%d.txt',fk)));
coh_sr_MA = importdata(strcat(dir_base,sprintf('/coherence_sr_MA_fk_%d.txt',fk)));
coh_ms_AM = importdata(strcat(dir_base,sprintf('/coherence_ms_AM_fk_%d.txt',fk)));
coh_mr_AM = importdata(strcat(dir_base,sprintf('/coherence_mr_AM_fk_%d.txt',fk)));
coh_sr_AM = importdata(strcat(dir_base,sprintf('/coherence_sr_AM_fk_%d.txt',fk)));
f = importdata(strcat(dir_base,sprintf('/frequency_range_fk_%d.txt',fk)));

% --- mean coherences
mean_cho_ms_MA = mean(coh_ms_MA(:,:));  % modulator - sender
mean_cho_mr_MA = mean(coh_mr_MA(:,:));  % modulator - receiver
mean_cho_sr_MA = mean(coh_sr_MA(:,:));  % sender - receiver 

mean_cho_ms_AM = mean(coh_ms_AM(:,:));
mean_cho_mr_AM = mean(coh_mr_AM(:,:));
mean_cho_sr_AM = mean(coh_sr_AM(:,:));

% --- std coherences
std_cho_ms_MA = std(coh_ms_MA(:,:));  % modulator - sender
std_cho_mr_MA = std(coh_mr_MA(:,:));  % modulator - receiver
std_cho_sr_MA = std(coh_sr_MA(:,:));  % modulator - receiver

std_cho_ms_AM = std(coh_ms_AM(:,:));
std_cho_mr_AM = std(coh_mr_AM(:,:));
std_cho_sr_AM = std(coh_sr_AM(:,:));

% --- Error bars
err_ms_MA = std_cho_ms_MA/sqrt(48);
err_mr_MA = std_cho_mr_MA/sqrt(48);
err_sr_MA = std_cho_sr_MA/sqrt(20);

err_ms_AM = std_cho_ms_AM/sqrt(48);
err_mr_AM = std_cho_mr_AM/sqrt(48);
err_sr_AM = std_cho_sr_AM/sqrt(20);


set(0,'DefaultFigureVisible','on')

% -- FIGURE: Plot average coherence across sessions for MR, SR, MS
fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all

% shadedErrorBar(f,mean_cho_ms_MA,err_ms_MA,'lineProps','b','patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_mr_MA,err_mr_MA,'lineprops',{'color',[255, 83, 26]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_sr_MA,err_sr_MA,'lineprops',{'color',[230, 184 , 0]/255},'patchSaturation',0.4); hold on


shadedErrorBar(f,mean_cho_ms_AM,err_ms_AM,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_mr_AM,err_mr_AM,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_sr_AM,err_sr_AM,'lineprops',{'color',[0 204 204]/255},'patchSaturation',0.4); hold on

grid on
title('Mean Abs coherence of MS, SR, SR','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('M-S mean abs','M-R mean abs','S-R mean abs','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% xlim([0 60])
set(gcf, 'Position',  [100, 600, 1000, 600])


fname = strcat(dir_base,sprintf('/coherence_MA_mean_causal_MR_SR_fk_%d.png',fk));
fname = strcat(dir_base,sprintf('/coherence_MA_mean_causal_MR_SR_fk_%d_zoom.png',fk));
saveas(fig,fname)


