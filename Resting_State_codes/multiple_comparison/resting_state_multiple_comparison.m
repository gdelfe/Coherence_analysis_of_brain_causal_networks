
clear all; close all;

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/vol/sas8/Maverick_RecStim_vSUBNETS220/160125/004')
addpath('/mnt/pesaranlab/People/Gino/DL-modulators/Gino_codes')

file = 'rec004.Frontal.lfp.dat'
fid = fopen(file);
format = 'float=>single';

CH = 220; % tot number of channels
FS = 1000; % sampling
Sess = 15; % Session number


data = fread(fid,[CH,inf],format); % load the RS data
% h = fread(fid,[CH,diff(bn)*FS./1e3],format);
% ---- bipolar referencing, pairs of electrodes
elect_dir = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d',Sess);
if ~exist(elect_dir, 'dir')
    mkdir(elect_dir)
end
electrode = importdata(strcat(elect_dir,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairs
receiver =  [196 199]; %[37 38];  %Data.StimResPairs;  % ---- receiver pair
sender = [43 44]; %Data.StimPairs.Syllable; % ---- sender pair


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
% --------- COHERENCE-GRAM --------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- parameters for the coherence-gram
tapers = [4 4];
N = tapers(1);
nt = tot_time;
dn = 0.01;
fs = 1000;
nwin = single(floor((nt-N*fs)/(dn*fs)))

% --- coherence
[c_sr,tf,f,spec_r,spec_s] = tfcoh_GINO(lfp_R,lfp_S,tapers,1e3,dn,60,2,[],[],1);

% Ch = 34;
% [c_ms,tf,f,spec_r,spec_s] = tfcoh_GINO(lfpRS(Ch,:),lfp_S,tapers,1e3,dn,60,2,[],[],1); % coherence modulator-sender 
% [c_mr,tf,f,spec_r,spec_s] = tfcoh_GINO(lfpRS(Ch,:),lfp_R,tapers,1e3,dn,60,2,[],[],1); % coherence modulator-receiver 

% -- Figure: coherence spectrum
fig = figure; tvimage(abs(c_sr(:,:))); colorbar; % coherence spectrum
% fig = figure; tvimage(abs(c_ms(:,:))); colorbar; % coherence spectrum
% fig = figure; tvimage(abs(c_mr(:,:))); colorbar; % coherence spectrum

xticks = floor(linspace(1,length(tf),5));
xticklabels = tf(xticks);
xtickformat('%d')
% yticks = floor(linspace(1,length(f),10));
yticks = 1:50:length(f);

yticklabels = floor(f(yticks));
ytickformat('%.2f')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
title('Coherence-gram S-R','FontSize',12);
xlabel('time (sec)');
ylabel('freq (Hz)')
% ylim([0,120])
set(gcf, 'Position',  [100, 600, 1000, 600])


% fname = strcat(dir,'/coherence-gram_SR.jpg');
% saveas(fig,fname);

% fname = strcat(elect_dir,sprintf('/coherence-gram_zoom_fq_%d_%d_2_thresh.jpg',fmin,fmax));
% saveas(fig,fname);
%
% fname = strcat(elect_dir,sprintf('/coherence-gram_zoom_fq_%d_%d_2_thresh_zoom.jpg',fmin,fmax));
% saveas(fig,fname);

% --- 1D array for p-values for FDR. Each value is a p-value cluster corrected for one channel
p_channel_AM = ones(1,size(lfpRS,1));
p_channel_MA = ones(1,size(lfpRS,1));

% keyboard

%% %%%%%%%%%%%%%%%%%%%%%
% -- Modulator Channel

Ch = 34; % Channel being analyzed

% for Ch = 3:52
    
    close all
    
    % directory path to save files
    dir = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d/p_th_0.05/Ch_%d',Sess,Ch);
    
    if ~exist(dir, 'dir')
        mkdir(dir)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Modulator's Spectrogram
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % parameters used by Shaoyu
    k = 4;
    fk = [0 60];
    tapers = [0.5 5];
    dn = 0.005;
    fs = FS;
    pad = 2;
    
    
    [specRS, fRS , tiRS] = tfspec_GINO(lfpRS(Ch,:),tapers,fs,dn,fk,pad,0.05,0,1);
    fig = figure; tvimage(sq(log(specRS(1:500,:)))); title(sprintf('RS - N = %.2f, W = %d, dn = %.4f, k = %d',tapers(1),tapers(2),dn,k)); colorbar;
    
    % fname = strcat(dir,sprintf('/spectrogram_Ch_%d.jpg',Ch));
    % saveas(fig,fname);
    
    
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ---------------------------------------  %%%%
    % -----          MODULATOR SCORE        -----   %
    %%% ---------------------------------------  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
    display([sprintf('Modulator score Ch %d',Ch)])
    % -- parameters
    fmin = 10;
    fmax = 40;
    Delta_ms = 0; % time (ms) at which to start computing spectrogram. Default is 0 which means 0 ms, beginning of the time series.
    
    N = tapers(1);
    nt = 1000; % time length (ms) for the baseline in the stimulation experiment
    nwin = single(floor((nt-N*fs)./(dn*fs))); % numb of moving window (# points) in the resultig
    
    
    bin = (nt-N*fs)/nwin; % how many milliseconds correspon to a bin in the spectrogram
    Delta_bin = Delta_ms/bin;
    Nstep = (N*fs)/bin; % how many step correspond to a temporal shift of N
    % step = Nstep + nwin; % shift of N + nwin. To go from end of the window to mid of next window
    step = 110; % with steps < 100 there are overlapping windows for the computation of the beta power
    ms = []; % modulator score
    for dt = (Delta_bin + nwin/2 +1):step:(size(specRS,1)-(nwin+Nstep)/2) % from mid window + 1, to end of time - N/2 - nwin/2, with step of nwin + N
        ms = [ms, mean2(log(specRS(dt:dt+nwin/2,fmin:fmax)))]; % mean of the log beta power in the range t = [500, 750] - keep in mind the spectrogram used a range [500,1000]
        %     ms = [ms, mean2(log(specRS(dt:dt+(nwin+Nstep)/2,fmin:fmax)))]; % mean of the log beta power in the range t = [500,1000] - keep in mind the spectrogram in this case used a range [500,1250]
        
    end
    
    dt = (Delta_bin + nwin/2 +1):step:(size(specRS,1)-(nwin+Nstep)/2); % BINS: from mid window + 1, to end of whole time - N/2 - nwin/2, with step of nwin + N
    ts = tiRS(dt+(nwin+Nstep)/2); % Starting time (ms) for the computation of the cross-correlation sender-receiver
    
    % just for the records, never used
    tspec_beg = tiRS(dt); % Initial window-time (ms) for the computation of the beta power
    tspec_end = tiRS(dt+nwin/2); % Final window-time (ms) for the computation of the beta power
    
    % % -- FIGURE -----------%
    % % Modulator score vs time
    % % --------------------- %
    % fig = figure;
    % plot(ts,ms,'-')
    % title('Modulator Score vs time ','FontSize',13);
    % xlabel('time (ms)','FontSize',13)
    % ylabel('# modulator score','FontSize',13)
    % % xticklabels = sprintf('%1.1e',0:2e+4:tot_time);
    % % % xticklabels = sprintf('%1.1e',xticklabels)
    % % xticks = xticklabels;
    % % set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    % grid on
    % legend('ms')
    % set(gcf, 'Position',  [100, 600, 800, 600])
    
    
    % fname = strcat(dir,sprintf('/modulator_score_vs_time_fq_%d_%d.jpg',fmin,fmax));
    % saveas(fig,fname);
    
    
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Histograms and CDF
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % % ------ Histogram
    % fig = figure;
    % h = histogram(ms,50,'Normalization','probability','FaceAlpha',.6)
    % title('Histogram - Modulator Score ','FontSize',13);
    % xlabel('modulator score','FontSize',13)
    % ylabel('# normalized counts','FontSize',13)
    % legend('ms')
    % grid on
    % hold off
    %
    % fname = strcat(dir,sprintf('/histogram_modulator_score_fq_%d_%d.jpg',fmin,fmax));
    % saveas(fig,fname);
    
    
    msSort = sort(ms); % sort ms to compute the cumulative
    
    % % -------- Cumulative distribution with bins
    % fig = figure
    % hcum = histogram(msSort,50,'Normalization','cdf')
    % title('Cumulative distribution ','FontSize',11);
    % xlabel('modulator score','FontSize',13)
    % ylabel('cumulative','FontSize',13)
    % grid on
    % %
    % % fname = strcat(dir,sprintf('/cdf_bins_modulator_score_fq_%d_%d.jpg',fmin,fmax));
    % % saveas(fig,fname);
    % % >>>>>> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% -------------------------------------------------------------- %%%%
    %   -----------           CROSS-CORRELATION  --------------------- %%%%
    %   ----  of Sender-Receiver for high and low modulator score --   %%%%
    %%% -------------------------------------------------------------- %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % keyboard
    
    %%%%%%%%%%%%%%%%%%%%%%
    % ---- One threshold
    % ------------------
    % th = 5.8; % threshold
    % msIndxHigh = find(ms >= th); % get the index for values of ms > threshold
    % msIndxLow = find(ms<th); % get the index for values of ms < threshold
    
    %%%%%%%%%%%%%%%%%%%%%%
    % ---- Two thresholds
    % -------------------
    
    msSort = sort(ms);
    
    X = floor(size(ms,2)*1/3); % 1/3 of the total modulator scores
    th1 = msSort(end-X+1); % higher threshold,
    th2 = msSort(X); % lower threshold,
    msIndxHigh = find(ms > th1); % get the index for values of ms > threshold
    msIndxLow = find(ms < th2);
    
    n = 10;
    delta = 100; % max lag shift in ms
    tinterval = floor(ts(n)):floor(ts(n))+delta; % from the moment the modulator has a certain value, i.e. ms(n), at ts(n), delta-ms onward
    XS = lfp_S(1,tinterval); % sender lfp delta-ms after modulator score is higher than threshold
    XR = lfp_R(1,tinterval); % receiver lfp delta-ms after modulator score is higher than threshold
    time = linspace(ts(n),ts(n)+delta,delta+1); % time for the x axis
    
    % % ------ cross-correlation (lagged)
    [C, lags] = xcorr(XS,XR,'coeff');
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Get cross-correlation for the whole time period
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % % ------- get lfp 100 ms after modulator score is high
    % delta = 100; % max lag shift in ms
    % C_high = [];
    % for n = msIndxHigh
    %     tinterval = floor(ts(n)):floor(ts(n))+delta;
    %     XS = lfp_S(1,tinterval); % sender lfp 100 ms after modulator score is higher than threshold
    %     XR = lfp_R(1,tinterval); % receiver lfp 100 ms after modulator score is higher than threshold
    %
    %     % cross-correlation (lagged)
    %     [C, ~] = xcorr(XS,XR,'coeff');
    %     C_high = [C_high; C];
    % end
    %
    % % -------- get lfp 100 ms after modulator score is low
    % C_low = [];
    % for n = msIndxLow
    %     tinterval = floor(ts(n)):floor(ts(n))+delta;
    %     XS = lfp_S(1,tinterval); % sender lfp 100 ms after modulator score is higher than threshold
    %     XR = lfp_R(1,tinterval); % receiver lfp 100 ms after modulator score is higher than threshold
    %
    %     % cross-correlation (lagged)
    %     [C, ~] = xcorr(XS,XR,'coeff');
    %     C_low = [C_low; C];
    % end
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Figures: cross-correlation for the whole period
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ------- cross-correlation high scores
    % fig = figure;
    % imagesc(C_high); colorbar; caxis([-1 1])
    % xticklabels = -100:20:100;
    % xticks = linspace(1, size(C_high, 2), numel(xticklabels));
    % set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    % title('cross-correlation for modulator high-scores','FontSize',12);
    % xlabel('lag (ms)','FontSize',13)
    % ylabel('cross-correlation','FontSize',13)
    % axis on
    %
    % % fname = strcat(dir,sprintf('/all_high_score_cross-correlations_fq_%d_%d.jpg',fmin,fmax));
    % % saveas(fig,fname);
    %
    % % ------- cross-correlation low scores
    % fig  = figure;
    % imagesc(C_low); colorbar; caxis([-1 1])
    % xticklabels = -100:20:100;
    % xticks = linspace(1, size(C_low, 2), numel(xticklabels));
    % set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    % title('cross-correlation for modulator low-scores','FontSize',12);
    % xlabel('lag (ms)','FontSize',13)
    % ylabel('cross-correlation','FontSize',13)
    % axis on
    
    % fname = strcat(dir,sprintf('/all_low_score_cross-correlations_fq_%d_%d.jpg',fmin,fmax));
    % saveas(fig,fname);
    
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Find max of the abs values for C_high and C_low and their index
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ccc
    % %%%%% --- Time intervals for correlations %%%%
    % tin = 71; % initial time (in steps) for the computation of the max in cross-correlation
    % tfin = 131; % final time (in steps) - as above
    % lags(tin:tfin)
    %
    % [Mh, indMh] = max(abs(C_high(:,tin:tfin)')); % high
    % maxSortHigh = fliplr(sort(Mh));
    % [Ml, indMl] = max(abs(C_low(:,tin:tfin)'));
    % maxSortLow = fliplr(sort(Ml));
    %
    % Xhigh = linspace(0,1,size(indMh,2));
    % Xlow = linspace(0,1,size(indMl,2));
    %
    %
    % % ----- Figure: max abs value of cross-correlation
    % fig = figure;
    % plot(Xhigh,maxSortHigh,'-o')
    % hold on
    % plot(Xlow,maxSortLow,'-o')
    % hold on
    % plot(Xhigh,repelem(mean(maxSortHigh),length(Xhigh)))
    % hold on
    % plot(Xlow,repelem(mean(maxSortLow),length(Xlow)))
    % legend('high-score','low-score','mean high','mean low')
    % title('max of the abs value for the cross-correlations (sorted)','FontSize',12);
    % ylabel('max(abs(cross-correlation))','FontSize',13)
    % xlabel('# of correlation points (normalized)','FontSize',13)
    % xlim([0,1.01])
    % grid on
    % hold off
    %
    %
    % % fname = strcat(dir,sprintf('/max_abs_values_cross-corr_sorted_fq_%d_%d_2_thresholds.jpg',fmin,fmax));
    % % saveas(fig,fname);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %  Find the extreme values for the cross-corr (either positive or negative)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % cLargeH = [];
    % for i = 1:size(C_high,1)
    %     cLargeH = [cLargeH, C_high(i,indMh(i))];
    % end
    %
    % cLargeL = [];
    % for i = 1:size(C_low,1)
    %     cLargeL = [cLargeL, C_low(i,indMl(i))];
    % end
    % % ---- sort the cross-correlation values
    % cLargeHsort = fliplr(sort(cLargeH));
    % cLargeLsort = fliplr(sort(cLargeL));
    %
    % % ----- Figure: extreme value of cross-correlation
    % fig = figure;
    % % plot(cLargeH)
    % % hold on
    % % plot(cLargeL)
    % % hold on
    % plot(Xhigh,cLargeHsort,'-o')
    % hold onc
    % plot(Xlow,cLargeLsort,'-o')
    % hold on
    % plot(Xhigh,repelem(mean(cLargeH),length(Xhigh)))
    % hold on
    % plot(Xlow,repelem(mean(cLargeL),length(Xlow)))
    % legend('high-score','low-score','mean high','mean low')
    % title('extreme values of cross-correlations (sorted)','FontSize',12);
    % ylabel('extreme(cross-correlation)','FontSize',13)
    % xlim([0,1.01])
    % grid on
    % hold off
    
    % % fname = strcat(dir,sprintf('/extreme_values_cross-corr_sorted_fq_%d_%d_th_%.1f.jpg',fmin,fmax,th));
    % % fname = strcat(dir,sprintf('/extreme_values_cross-corr_sorted_fq_%d_%d_2_thresholds.jpg',fmin,fmax));
    % % saveas(fig,fname);
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %  Find the max values for the cross-corr
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    % [Mh, indMh] = max(C_high(:,tin:tfin)'); % high
    % maxSortHigh = fliplr(sort(Mh));
    % [Ml, indMl] = max(C_low(:,tin:tfin)'); % low
    % maxSortLow = fliplr(sort(Ml));
    %
    % % ----- Figure: max value of cross-correlation
    % fig = figure;
    % plot(Xhigh,maxSortHigh,'-o')
    % hold on
    % plot(Xlow,maxSortLow,'-o')
    % hold on
    % plot(Xhigh,repelem(mean(maxSortHigh),length(Xhigh)))
    % hold on
    % plot(Xlow,repelem(mean(maxSortLow),length(Xlow)))
    % legend('high-score','low-score','mean high','mean low')
    % title('max value for the cross-correlations (sorted)','FontSize',12);
    % ylabel('max(abs(cross-correlation))','FontSize',13)
    % grid on
    % hold off
    %
    % % fname = strcat(dir,sprintf('/max_values_cross-corr_sorted_fq_%d_%d_th_%.1f.jpg',fmin,fmax,th));
    % % fname = strcat(dir,sprintf('/max_values_cross-corr_sorted_fq_%d_%d_2_thresholds.jpg',fmin,fmax));
    % % saveas(fig,fname);
    
    
    
    
    
    % Keyboard
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ---------------------------------------  %%%%
    %%% ----------- CHOERENCE ANALYSIS  -------- %%%%
    %%% ---------------------------------------  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    display([sprintf('Coherence Analysis, Ch %d',Ch)])

    
    fig = figure; plot(f,mean(abs(c_sr(:,:))));  grid on; title('Mean(abs(Coherence)) vs frequency S-R','FontSize',12); xlabel('freq (Hz)');
    % --- Fig: abs(mean()) plot mean coherence (averaged in time) vs frequency --
    hold on; plot(f,abs(mean(c_sr(:,:)))); % plot mean coherence in time vs frequency
    legend('mean(abs())','abs(mean))')
    
%     fig = figure; plot(f,mean(abs(c_ms(:,:))));  grid on; title('Mean(abs(Coherence)) vs frequency M-S','FontSize',12); xlabel('freq (Hz)');
%     % --- Fig: abs(mean()) plot mean coherence (averaged in time) vs frequency --
%     hold on; plot(f,abs(mean(c_ms(:,:)))); % plot mean coherence in time vs frequency
%     legend('mean(abs())','abs(mean))') 
%     
%     fig = figure; plot(f,mean(abs(c_mr(:,:))));  grid on; title('Mean(abs(Coherence)) vs frequency M-R','FontSize',12); xlabel('freq (Hz)');
%     % --- Fig: abs(mean()) plot mean coherence (averaged in time) vs frequency --
%     hold on; plot(f,abs(mean(c_mr(:,:)))); % plot mean coherence in time vs frequency
%     legend('mean(abs())','abs(mean))') 
    
    % % fname = strcat(dir,sprintf('coherence_vs_freq_MeanAbs_AbsMean_fq_%d_%d_2_thresh.jpg',fmin,fmax));
    % % saveas(fig,fname);
    %
    %
    % --- FIGURE --------- %%
    % -- coherence and modulator score --- %
    fig = figure;
    plot(tf,abs(c_sr(:,59))); % plot coherence at a given frequency (7 Hz)
    hold on
    plot(tf,mean(abs(c_sr(:,34:75)),2)); % plot coherence averaged in a range of frequencies vs time
    hold on
    plot(ts,ms - mean(ms));
    xlim([2000 50000])
    legend('coherence at 7Hz','coherence average in range ~[4,9] Hz','modulator score')
    title('Coherence vs time and modulator score','FontSize',12);
    grid on
    hold off
    set(gcf, 'Position',  [100, 600, 800, 600])
    
        
%     fname = strcat(dir,sprintf('/coherence_and_mod_score_vs_time_fq_%d_%d_2_thresh_step_%d.jpg',fmin,fmax,step));
%     saveas(fig,fname);
    
    
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ HIGH SCORE - LOW SCORE -- COHERENCE  -----
    %
    % Find the time tf of the sender-receiver at which the modulator is high or
    % low. The modulator time is ts, the sender-receiver time is indexed by tf
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % shift controls the time delay in milliseconds (forward), so shifth = 100
    % means that the coherence is computed 100 ms after the time at which the
    % modulator score is computed
    
    shift = 0;
    % -- for the high scores
    Hix = []; % time index for the coherence when ms is high
    for i = msIndxHigh
        [ d, ix ] = min( abs( tf-ts(i)-shift ) ); % find the closest tf time to ts (its index)
        Hix  = [Hix,ix];
    end
    
    % -- for the low scores
    Lix = []; % time index for the coherence when ms is low
    for i = msIndxLow
        [ d, ix ] = min( abs( tf-ts(i)-shift ) ); % find the closest tf time to ts (its index)
        Lix  = [Lix,ix];
    end
    
    
    
    % % -------------------------------- %
    % % --- FIGS: COHERENCE vs TIME ---- %
    %
    % -- Coherence (at a 7 Hz) vs time for high and score ms value
    
%     fig = figure;
%     plot(tf(Hix),abs(c_sr(Hix,59)) - mean(abs(c_sr(Hix,59))),'-'); % choerence high-ms at 7 Hz
%     hold on
%     plot(tf(Lix),abs(c_sr(Lix,59)) - mean(abs(c_sr(Lix,59))),'-');  % choerence low-ms at 7 Hz
%     hold on
%     plot(ts,ms - mean(ms),'-');
%     title(sprintf('Ch %d, abs of coherence (at 7hz) for high and low modulator score',Ch),'FontSize',12);
%     legend('coherence high-ms','coherence low-ms','modulator score')
%     xlabel('time (ms)','FontSize',13)
%     grid on
%     set(gcf, 'Position',  [100, 500, 800, 600])
%     hold off
    %
    %
    % % fname = strcat(dir,sprintf('/coherence_vs_time_at_7Hz_fq_%d_%d_2_thresh.jpg',fmin,fmax));
    % % saveas(fig,fname);
    
    
    
    % -------------------------------- %
    % --- FIGS: COHERENCE vs FREQ ---- %
    
    % ----- FIGURE: mean(abs(c_sr)) vs freq for high and low modulator scores
    fig = figure; plot(f,mean(abs(c_sr(Hix,:))))
    hold on
    plot(f,mean(abs(c_sr(Lix,:))))
    title(sprintf('mean(abs(...)) - delay: %d ms',shift));
    grid on
    legend('mean(abs(c_sr(High)))','mean(abs(c_sr(Low)))')
    xlim([0 30])
    
    fname = strcat(dir,sprintf('/mean_abs_cohe_vs_f_fq_%d_%d_2_thresh.jpg',fmin,fmax));
    saveas(fig,fname);
    
    % ----- FIGURE: abs(mean(c_sr)) for high and low modulator scores
    fig = figure; plot(f,abs(mean(c_sr(Hix,:))))
    hold on
    plot(f,abs(mean(c_sr(Lix,:))))
    grid on
    legend('abs(mean(c_sr(High)))','abs(mean(c_sr(Low)))')
    title(sprintf('abs(mean(...)) - delay: %d ms',shift));
    xlim([0 30])
    
    fname = strcat(dir,sprintf('/abs_mean_cohe_vs_f_fq_%d_%d_2_thresh.jpg',fmin,fmax));
    saveas(fig,fname);
    
    % ---- Difference c_sr(High) - c_sr(Low)
    result_AM = abs(mean(c_sr(Hix,:))) - abs(mean(c_sr(Lix,:))); % abs(mean(...))
    result_MA = mean(abs(c_sr(Hix,:))) - mean(abs(c_sr(Lix,:))); % mean(abs(...))
    
    % ----- FIGURE: Difference: coh(high) - coh(low)
    fig = figure;
    plot(f,result_AM)
    hold on
    plot(f,result_MA)
    title(sprintf('Coherence difference: High values - Low values, step = %.2f sec',step/200), 'FontSize',12)
    legend('diff[abs(mean(...))]','diff[mean(abs(..))]')
    grid on
    xlim([0 30])
    
    fname = strcat(dir,sprintf('/difference_cohe_vs_f_fq_%d_%d_step_%d.jpg',fmin,fmax,step));
    % fname = strcat(dir,sprintf('/difference_cohe.jpg'));
    saveas(fig,fname);
    
    % ------------------------  %%
    % -- PERMUTATION TEST ----  %%
    % ------------------------- %%
    display(['Running permutation test...'])
    iter = 1000; % number of iterations
    diff_AM = zeros(iter,size(c_sr,2));
    diff_MA = zeros(iter,size(c_sr,2));
    for j = 1:iter
        
        perm = randperm(length(tf));
        perm_Hix = perm(1:length(Hix));
        perm_Lix = perm(length(Hix)+1:length(tf));
        
        c_Hix_meanAbs = mean(abs(c_sr(perm_Hix,:)),1);
        c_Lix_meanAbs = mean(abs(c_sr(perm_Lix,:)),1);
        
        c_Hix_absMean = abs(mean(c_sr(perm_Hix,:),1));
        c_Lix_absMean = abs(mean(c_sr(perm_Lix,:),1));
        
        diff_absMean = c_Hix_absMean - c_Lix_absMean; % difference between high and low score values (abs(mean..))
        diff_meanAbs = c_Hix_meanAbs - c_Lix_meanAbs; % difference between high and low score values (mean(abs..))
        
        diff_AM(j,:) = diff_absMean;
        diff_MA(j,:) = diff_meanAbs;
        
    end
    display(['Permutation test DONE'])
    
   
    
    % FIGURE: histogram of the perm data
    % fig = figure; histogram(diff_AM(:,59),20,'Normalization','probability','FaceAlpha',.6); grid on
    % hold on; histogram(diff_MA(:,59),20,'Normalization','probability','FaceAlpha',.6); grid on
    % legend('diff AM','diff MA')
    % title('Distribution of perm data','FontSize',12)
    %
    % fname = strcat(dir,'/distributions_of_perm_data.jpg');
    % saveas(fig,fname);
    %
    % fname = strcat(dir,sprintf('/distribution_of_permuted_data_fq_%d_%d.jpg',fmin,fmax));
    % saveas(fig,fname);
    
    % --- FIGURE: mean of the differences across samples vs frequency
    fig = figure;
    plot(f,mean(diff_AM));
    hold on
    plot(f,mean(diff_MA));
    legend('diff[abs(mean(..))]', 'diff[mean(abs(...)]')
    grid on
    title(sprintf('mean of the diff coh(high) - coh(low), perm data, step = %.2f sec',step/200), 'FontSize',10);
    
    fname = strcat(dir,sprintf('/mean_of_coherence_difference_fq_%d_%d_step_%d.jpg',fmin,fmax,step));
    saveas(fig,fname);
    
    % --- FIGURE: std of the differences across samples vs frequency
    fig = figure;
    plot(f,std(diff_AM));
    hold on
    plot(f,std(diff_MA));
    title(sprintf('Std of the difference coh(high) - coh(low), perm data, step = %.2f sec',step/200), 'FontSize',10);
    titleHandle = get( gca ,'Title' );
    pos  = get( titleHandle , 'position' );
    set( titleHandle , 'position' , pos + [0. 0.0005 0] );
    legend('abs mean', 'mean abs')
    grid on
    
    fname = strcat(dir,sprintf('/std_of_coherence_difference_fq_%d_%d_step_%d.jpg',fmin,fmax,step));
    saveas(fig,fname);
    
    
    
    % ------------------------  %%
    % -----   Z-SCORES     ----  %%
    % ------------------------- %%
    
    % -- zscores for each frequency
    zscore_AM = (result_AM - mean(diff_AM))./std(diff_AM);
    zscore_MA = (result_MA - mean(diff_MA))./std(diff_MA);
    
    dlmwrite(strcat(dir,sprintf('/zscore_vs_freq_AM_step_%d.txt',step)),zscore_AM');
    dlmwrite(strcat(dir,sprintf('/zscore_vs_freq_MA_step_%d.txt',step)),zscore_MA');
    
    % -- FIGURE: zscores vs frequency
    fig = figure;
    plot(f,zscore_AM)
    hold on
    plot(f,zscore_MA)
    hold on
    title(sprintf('zscore vs frequency, step = %d = %.2f sec',step,step/200),'FontSize',12)
    legend('zscore AM','zscore MA')
    ylabel('z-score')
    xlabel('frequency')
    xlim([0,30])
    grid on
    
    fname = strcat(dir,sprintf('/zscore_vs_freq_fq_%d_%d_step_%d.jpg',fmin,fmax,step));
    saveas(fig,fname);
    
    
    
    % -- FIGURE: histogram of zscores across frequencies
    % fig = figure; histogram(zscore_AM,30,'Normalization','probability','FaceAlpha',.6); grid on
    % hold on; histogram(zscore_MA,30,'Normalization','probability','FaceAlpha',.6); grid on
    % legend('zscore AM','zscore MA')
    % title('z distribution across freq','FontSize',12)
    
    % fname = strcat(dir,sprintf('/zscore_histog_fq_%d_%d.jpg',fmin,fmax));
    % saveas(fig,fname);
    
    
    
    % ------------------------  %%
    % -----   P-VALUES     ---  %%
    % ------------------------- %%
    
    % -- Assuming the zscore distribution is gaussian, two tails test
    pval_AM = 2*normcdf(-abs(zscore_AM))
    pval_MA = 2*normcdf(-abs(zscore_MA))
    
    
    %-- FIGURE: p-value vs frequency, assuming gaussian
    fig = figure;
    semilogy(f,pval_AM)
    hold on
    semilogy(f,pval_MA)
    grid on
    title(sprintf('p-value vs frequency, step = %.2f sec, assuming gaussian',step/200),'FontSize',10)
    legend({'p-val AM','p-val MA'},'Location','southeast')
    xlim([0,30])
    
    fname = strcat(dir,sprintf('/pval_vs_freq_fq_%d_%d_step_%d.jpg',fmin,fmax,step));
    saveas(fig,fname);
    
    % keyboard
    
    % -- FIGURE: histogram of zscores across frequencies
    fig = figure;
    histogram(pval_AM,20,'Normalization','probability','FaceAlpha',.6); grid on
    hold on; histogram(pval_MA,20,'Normalization','probability','FaceAlpha',.6); grid on
    legend('p-val AM','p-val MA')
    title('p-val distribution','FontSize',12)
    
    fname = strcat(dir,sprintf('/pval_histog_fq_%d_%d_step_%d.jpg',fmin,fmax,step));
    saveas(fig,fname);
    
    
    % % --- Empirical P-VALUE
    % pcount_AM = zeros(1,size(diff_AM,2));
    % pcount_MA = zeros(1,size(diff_MA,2));
    % for freq = 1:size(diff_AM,2)
    %
    %     % -- two tailed, symmetric
    %     pcount_AM(freq) = nnz(diff_AM(:,freq) > abs(result_AM(freq)) | diff_AM(:,freq) < -abs(result_AM(freq)) )/iter;
    %     pcount_MA(freq) = nnz(diff_MA(:,freq) > abs(result_MA(freq)) | diff_MA(:,freq) < -abs(result_MA(freq)) )/iter;
    %
    % end
    %
    % %-- FIGURE: p-value vs frequency, empirical
    % fig = figure;
    % semilogy(f,pcount_AM)
    % hold on
    % semilogy(f,pcount_MA)
    % grid on
    % title(sprintf('p-val vs frequency, empirical, step = %.2f',step/200),'FontSize',12)
    % legend({'p-val AM','p-val MA'},'Location','southeast')
    % xlim([0,30])
    % % ylim([0.0005,1])
    %
    %
    % fname = strcat(dir,sprintf('/pval_vs_freq_empirical_fq_%d_%d_step_%d.jpg',fmin,fmax,step));
    % saveas(fig,fname);
    
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CLUSTER CORRECTION  %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    p_th = 1e-03; % p-value threshold
    z_th =norminv((0.5*p_th)) % get zscore threshold from p-val treshold
    p = @(zs) 2*normcdf(-abs(zs)) % p-value from z-score, function
    
    f_lim = 246; % f bin limit: -> consider only frequencies up to 30 Hz
    
    % ---- OBSERVED VALUE ---- %
    
    [zmax_AM_Obs,id_zmax_AM_Obs,blocks_F_AM] = cluster_correction(result_AM,diff_AM,f,f_lim,z_th); % Abs Mean
    [zmax_MA_Obs,id_zmax_MA_Obs,blocks_F_MA] = cluster_correction(result_MA,diff_MA,f,f_lim,z_th); % Mean Abs
    
    % --- Print data of the clusters: value of max(sum(z-score)), index of the
    % block with max z-score, frequency range of the max z-score cluster
    
    cluster_data_AM = [zmax_AM_Obs,id_zmax_AM_Obs,blocks_F_AM{id_zmax_AM_Obs}];
    fid = fopen(strcat(dir,sprintf('/cluster_data_AM_step_%d.txt',step)),'w');
    fprintf(fid, '%.3f\n',cluster_data_AM) ;
    fclose(fid) ;
    
    cluster_data_MA = [zmax_MA_Obs,id_zmax_MA_Obs,blocks_F_MA{id_zmax_MA_Obs}];
    fid = fopen(strcat(dir,sprintf('/cluster_data_MA_step_%d.txt',step)),'w');
    fprintf(fid, '%.3f\n',cluster_data_MA) ;
    fclose(fid) ;
    
    
    
    % fig = figure;
    % plot(f,zscore_AM)
    % hold on
    % plot(f,zscore_MA)
    % hold on
    % plot(blocks_F_AM{id_zmax_Obs_AM},repelem(abs(z_th),size(blocks_F_AM{id_zmax_Obs_AM},2)))
    % hold on
    % plot(blocks_F_AM{id_zmax_Obs_AM},repelem(-abs(z_th),size(blocks_F_AM{id_zmax_Obs_AM},2)))
    % hold on
    % plot(blocks_F_MA{id_zmax_Obs_MA},repelem(abs(z_th),size(blocks_F_MA{id_zmax_Obs_MA},2)))
    % hold on
    % plot(blocks_F_MA{id_zmax_Obs_MA},repelem(-abs(z_th),size(blocks_F_MA{id_zmax_Obs_MA},2)))
    % grid on
    % title(sprintf('zscore vs frequency, step = %d = %.2f sec',step,step/200),'FontSize',12)
    % legend('zscore AM','zscore MA')
    % ylabel('z-score')
    % xlabel('frequency')
    % xlim([0,30])
    
    display([sprintf('p-values cluster, Ch %d',Ch)])

    % ------------------------ %
    % --- PERMUTED VALUES ---- %
    % ------------------------ %
    
    zmaxAM_Perm = zeros(1,iter);
    zmaxMA_Perm = zeros(1,iter);
    for i=1:iter
        [zmaxAM_Perm(i),id_zmax] = cluster_correction(diff_AM(i,:),diff_AM,f,f_lim,z_th);
        [zmaxMA_Perm(i),id_zmax] = cluster_correction(diff_MA(i,:),diff_MA,f,f_lim,z_th);
    end
    
    dlmwrite(strcat(dir,sprintf('/zmax_AM_perm_step_%d.txt',step)),zmaxAM_Perm');
    dlmwrite(strcat(dir,sprintf('/zmax_MA_perm_step_%d.txt',step)),zmaxMA_Perm');
    
%     zmaxAM_Perm = importdata(strcat(dir,sprintf('/Ch_%dzmax_AM_perm.txt',Ch)));
%     zmaxMA_Perm = importdata(strcat(dir,sprintf('/Ch_%dzmax_MA_perm.txt',Ch)));
    
    % Histogram of the z-score max
    fig = figure;
    histogram(abs(zmaxAM_Perm),30,'Normalization','probability','FaceAlpha',.6); grid on
    hold on; histogram(abs(zmaxMA_Perm),30,'Normalization','probability','FaceAlpha',.6); grid on
    legend('z-score cluster size AM','z-score cluster size MA')
    title('z-score cluster size','FontSize',12)
    
    fname = strcat(dir,sprintf('/zscore_max_cluster_distribution_step_%d.png',step));
    saveas(fig,fname);
    
    pClust_AM = nnz(zmaxAM_Perm > abs(zmax_AM_Obs) | zmaxAM_Perm < -abs(zmaxAM_Perm))/iter;
    pClust_MA = nnz(zmaxMA_Perm > abs(zmax_MA_Obs) | zmaxMA_Perm < -abs(zmaxMA_Perm))/iter;
    
    
    p_channel_AM(Ch) = pClust_AM;
    p_channel_MA(Ch) = pClust_MA;


% dlmwrite(strcat(elect_dir,sprintf('/p_channel_AM_step_%d.txt',step)),p_channel_AM');
% dlmwrite(strcat(elect_dir,sprintf('/p_channel_MA_step_%d.txt',step)),p_channel_MA');

% p_channel_AM = importdata(strcat(elect_dir,sprintf('/p_channel_AM_step_%d.txt',step)));
% p_channel_MA = importdata(strcat(elect_dir,sprintf('/p_channel_MA_step_%d.txt',step)));







