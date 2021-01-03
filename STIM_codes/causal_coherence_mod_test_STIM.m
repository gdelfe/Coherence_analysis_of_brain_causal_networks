
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the STIM coherence between the causal modulators found by
% Shaoyu's and the receiver
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators and plots this coherence for the
% hits and the misses separately
%
% INPUT: file with session modulator info
%        .mat file with structure AM and MA information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all

set(0,'DefaultLineLineWidth',2)

subjects = {'maverick','archie'};

%for
iSubject = 1% : length(subjects) % Loop on the animals
%     clearvars -except subjects iSubject
if strcmp(subjects{iSubject},'archie')
    archie_vSUBNETS220_rig3
else
    maverick_vSUBNETS220_rig3
end
PreStimSess = PreStimResponseAll_Database_NetworkEdge;

for iSess = 1 : numel(PreStimSess)
    %         pAccLLRperm = PreStimSess{iSess}{end-2};
    %         useSessIndx(iSess) = pAccLLRperm <= 0.05 & ~isnan(pAccLLRperm);
    
    pAccLLRperm = PreStimSess{iSess}{end-2};
    pFDR_logic = PreStimSess{iSess}{end-1};
    useSessIndx(iSess) = pFDR_logic & ~isnan(pAccLLRperm);
end

UsedSess = find(useSessIndx);

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Stim = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';
step = 110;

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% -- load structure files
newAM = load(strcat(dir_RS,'/session_AM.mat'))
session_AM = newAM.session_AM;

newMA = load(strcat(dir_RS,'/session_MA.mat'))
session_MA = newMA.session_MA;


% % -- print structures on stdout
% %format short
% for s=1:size(sess_info{1},1)
%     session_AM(s)
%     session_MA(s)
% end

% -- matrices to store coherence for each mod and then compute
% average
t_tot = 60;
ftot = 122;
coh_mr = zeros(48,t_tot,ftot); % coherence mod-receiver.# of causal mod, # time, # frequency
coh_mr_H = zeros(48,t_tot,ftot); % Hits
coh_mr_M = zeros(48,t_tot,ftot); % Misses

S_m = zeros(48,t_tot,ftot); % spectrum modulator, # of causal mod, # time, # frequency
S_r = zeros(48,t_tot,ftot); % spectrum receiver 
S_m_H = zeros(48,t_tot,ftot); % Hits
S_r_H = zeros(48,t_tot,ftot); % Misses
S_m_M = zeros(48,t_tot,ftot); % Hits
S_r_M = zeros(48,t_tot,ftot); % Misses


cnt = 1;

for i=1:size(sess_info{1},1)  % For all the session with a modulator
    
    Sess = sess_info{1}(i); % Session number
    
    disp(['Session ' num2str(Sess) ' out of ' num2str(length(PreStimSess)) ' ...'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD DATA ...                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    RespPair = sessElectrode(PreStimSess{Sess}); % responding channel
    stimName = PreStimSess{Sess}{9};
    stimTask = PreStimSess{Sess}{7};
    day = sessDay(PreStimSess{Sess});
    
    
    % % loading Pre data
    dataDir_Pre = sprintf('%s/AccLLR/%sStimAllSess/StimResponseSessions/',DATADIR,stimName);
    switch stimTask
        case 'StimSinglePulse'
            fileName_Pre = sprintf('%sSess%03d_%s_AccLLR_Elec%03d-Elec%03d_%s_1stPulse.mat',dataDir_Pre,Sess,day,RespPair(1),RespPair(2),stimName);
            
        case 'StimBlockShort'
            fileName_Pre = sprintf('%sSess%03d_%s_AccLLR_Elec%03d-Elec%03d_%s_grouped.mat',dataDir_Pre,Sess,day,RespPair(1),RespPair(2),stimName);
    end
    
    tic
    disp('Loading Pre data ...')
    load(fileName_Pre)
    disp('Done with Pre data loading')
    toc
    
    % % parameters %%%%%%
    fs = Data.Fs.lfp;% lfp sampling rate
    Fs = Data.Fs.raw;% raw sampling rate
    
    
    AnalParams = Data.Params.Anal;
    AnalParams.Tapers = [0.5,2];
    AnalParams.TestSpecDiff.fk = [10 40];
    
    
    fkNames = {'\beta'};
    Data.Params.Anal = AnalParams;
    Data.Spec.ROC.fk = AnalParams.TestSpecDiff.fk;
    StimTrials = Data.StimTrials(Data.goodTrials_index);
    sys = StimTrials(1).MT;
    bn_Pre = [-1005 -5]; % ms
    
    
    
    % % extract AccLLR results
    Results = Data.AccLLR.Results;
    EventST = Results.EventST;
    
    
    if ~isequal(sum(~isnan(EventST)),0) % if detected
        nFreqBands = size(AnalParams.TestSpecDiff.fk,1);
        
        if strcmp(Data.Spec.recordConfig,'Bipolar')
            lfp_Detected = Data.spec.lfp.bipolar.Detected;
            lfp_notDetected = Data.spec.lfp.bipolar.notDetected;
            
            Lfp_Pre = [];
            for iFreqBand = 1 : nFreqBands
                pval = [];
                pval_roc = [];
                AnalParams.Spec.Test.fk = AnalParams.TestSpecDiff.fk(iFreqBand,:);
                Fk = AnalParams.Spec.Test.fk;
                
                
                %                 dlmwrite(sprintf('recorded_pairs_modulators_Sess_%d.txt',iSess),Data.RecordPair,'delimiter',',');
                
                %  ---- Data.Spec.ROC.sigChIndx{1}
                % ---- find(Data.Spec.ROC.sigChIndx{1})
                if ~isempty(Data.Spec.ROC.sigChs{iFreqBand})
                    % load LFPs
                    if isempty(Lfp_Pre)
                        try
                            disp('Loading lfp data ...')
                            [Lfp_Pre] = trialStimPulseLfp(StimTrials, sys, [], [], 'PulseStarts', bn_Pre); % returns monopolar recording
                        catch
                            disp('Loading raw data ...')
                            [Raw_Pre] = trialStimPulseRaw(StimTrials, sys, [], [], 'PulseStarts', bn_Pre); % returns monopolar recording
                            
                            for iCh = 1 : size(Raw_Pre,2)
                                lfp_Pre(:,iCh,:) = mtfilter(sq(Raw_Pre(:,iCh,:)),[0.0025,400],Fs,0);
                            end
                            Lfp_Pre = lfp_Pre(:,:,1:Fs/fs:end);
                            clear Raw_Pre lfp_Pre
                        end
                    end
                    
                    fprintf('\n');
                    disp(['Start ROC on PSDs @ ' num2str(AnalParams.TestSpecDiff.fk(iFreqBand,1)) '-' num2str(AnalParams.TestSpecDiff.fk(iFreqBand,2)) ' Hz'])
                    fprintf('\n\n');
                    reverseStr = '';
                    
                end
            end
        end
    end
    
    % %%%%%%%%%%%%%%%%% --->>> LFP data loaded at this point 
    
    keyboard
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',   out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    
    % ---- bipolar referencing, pairs of electrodes
    dir_Sess = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data/Sess_%d',Sess);
    dir_RS_Sess = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state/Sess_%d',Sess);
    
    if ~exist(dir_Sess, 'dir')
        mkdir(dir_Sess)
    end
    
    load(strcat(dir_RS_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    
    % -- load list electrodes, sender, receiver
    electrode = sess_data.RecordPair; % ---- all electrode pairs
    receiver = sess_data.receiver_pair;  % ---- receiver pair
    sender = sess_data.sender_pair; % ---- sender pair

    
    % ---  freq parameter for the masking
    fmin = 10;
    fmax = 40;
    
    %% %%%%%% Assign LFP %%%%%%
    lfp_E = sq(Lfp_Pre(:,electrode(:,1),:) - Lfp_Pre(:,electrode(:,2),:)); % modulator lfp
    lfp_R = sq(Lfp_Pre(:,receiver(1),:) - Lfp_Pre(:,receiver(2),:)); % receiver lfp
    lfp_S = sq(Lfp_Pre(:,sender(1),:) - Lfp_Pre(:,sender(2),:)); % sender lfp
     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --------- COHERENCE-GRAM --------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ---- parameters for the coherence-gram
    tapers = [0.7 8];
    N = tapers(1);
    nt = 1000;
    dn = 0.005;
    fs = 1000;
    fk = 60;
    nwin = single(floor((nt-N*fs)/(dn*fs)))
    k = floor(2*N*tapers(2)-1)
    
    % --- coherence
    
    mod_Ch = session_AM(i).mod_idx; % causal modulators channel
    
    indx_list = [];
    
    for Ch = mod_Ch % for all the modulators in the session
        
        close all
        
        indx_list = [indx_list, cnt]; % store the cnt number --- needed for multiple plotting 
        
        dir_Ch = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d/p_th_0.05/Ch_%d',Sess,Ch);
        
        % -- labels Hits and Misses
        hitIndx = Data.spec.lfp.DetectedIndx{Ch}; % labels for the hits (which trial was a hit)
        missIndx = Data.spec.lfp.notDetectedIndx{Ch}; % labels for the misses (which trial was a miss)
        
        
        
        display(['Computing coherence-gram...'])

        % Compute coherences and spectrums, all, hits, and misses 
        [c_mr,tf,f,spec_m,spec_r] = tfcoh_GINO(sq(lfp_E(:,Ch,:)),lfp_R(:,:),tapers,1e3,dn,fk,2,0.05,11); % coherence modulator-receiver
        [c_mr_H,tf,f,spec_m_H,spec_r_H] = tfcoh_GINO(sq(lfp_E(hitIndx,Ch,:)),lfp_R(hitIndx,:),tapers,1e3,dn,fk,2,0.05,11); % Hits
        [c_mr_M,tf,f,spec_m_M,spec_r_M] = tfcoh_GINO(sq(lfp_E(missIndx,Ch,:)),lfp_R(missIndx,:),tapers,1e3,dn,fk,2,0.05,11); % Misses

        
        %     dlmwrite(strcat(dir_Sess,'/sr_coherogram.txt'),c_sr,'delimiter',' ');
        % -- Figure: coherence spectrum
        
        
        % --- FIGURE: MODULATOR-RECEIVER COHEROGRAM
        fig_mr = figure; tvimage(abs(c_mr(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:10:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('M-R Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Sess,sprintf('/MR_coherogram_ch_%d_N_%.2f_W_%d.jpg',Ch,tapers(1),tapers(2)));
        saveas(fig_mr,fname);
        

        % --- FIGURE: HITS MODULATOR-RECEIVER COHEROGRAM
        fig_mr = figure; tvimage(abs(c_mr_H(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:10:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('HITS M-R Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Sess,sprintf('/MR_coherogram_HITS_ch_%d_N_%.2f_W_%d.jpg',Ch,tapers(1),tapers(2)));
        saveas(fig_mr,fname);
        
        
        % --- FIGURE: MISSES MODULATOR-RECEIVER COHEROGRAM
        fig_mr = figure; tvimage(abs(c_mr_M(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:10:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('HITS M-R Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Sess,sprintf('/MR_coherogram_MISSES_ch_%d_N_%.2f_W_%d.jpg',Ch,tapers(1),tapers(2)));
        saveas(fig_mr,fname);
        
        
        
        % -- write coherogram
        dlmwrite(strcat(dir_Sess,sprintf('/coherogram_MR_ch_%d_N_%.2f_W_%d.txt',Ch,tapers(1),tapers(2))),c_mr,'delimiter',' ');
        dlmwrite(strcat(dir_Sess,sprintf('/coherogram_MR_Hits_ch_%d_N_%.2f_W_%d.txt',Ch,tapers(1),tapers(2))),c_mr_H,'delimiter',' ');
        dlmwrite(strcat(dir_Sess,sprintf('/coherogram_MR_Misses_ch_%d_N_%.2f_W_%d.txt',Ch,tapers(1),tapers(2))),c_mr_M,'delimiter',' ');
        
        % -- write spectrogram 
        dlmwrite(strcat(dir_Sess,sprintf('/spectrogram_M_ch_%d_N_%.2f_W_%d.txt',Ch,tapers(1),tapers(2))),spec_m,'delimiter',' ');
        dlmwrite(strcat(dir_Sess,sprintf('/spectrogram_M_Hits_ch_%d_N_%.2f_W_%d.txt',Ch,tapers(1),tapers(2))),spec_m_H,'delimiter',' ');
        dlmwrite(strcat(dir_Sess,sprintf('/spectrogram_M_Misses_ch_%d_N_%.2f_W_%d.txt',Ch)),spec_m_M,'delimiter',' ');
        
        dlmwrite(strcat(dir_Sess,sprintf('/spectrogram_R_ch_%d_N_%.2f_W_%d.txt',Ch,tapers(1),tapers(2))),spec_r,'delimiter',' ');
        dlmwrite(strcat(dir_Sess,sprintf('/spectrogram_R_Hits_ch_%d_N_%.2f_W_%d.txt',Ch,tapers(1),tapers(2))),spec_r_H,'delimiter',' ');
        dlmwrite(strcat(dir_Sess,sprintf('/spectrogram_R_Misses_ch_%d_N_%.2f_W_%d.txt',Ch,tapers(1),tapers(2))),spec_r_M,'delimiter',' ');
        
        % -- store coherence 
        coh_mr(cnt,:,:) = c_mr;
        coh_mr_H(cnt,:,:) = c_mr_H;
        coh_mr_M(cnt,:,:) = c_mr_M;
        
        % -- store spectrum modularor and receiver 
        S_m(cnt,:,:) = spec_m;
        S_m_H(cnt,:,:) = spec_m_H;
        S_m_M(cnt,:,:) = spec_m_M;
        
        S_r(cnt,:,:) = spec_r;
        S_r_H(cnt,:,:) = spec_r_H;
        S_r_M(cnt,:,:) = spec_r_M;
        
        cnt = cnt + 1;
        
    end % --- end of all modulator channels
    
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MODULATOR - RECEIVER COHERENCE FIGURES   %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % --- FIGURE --------- %%
    % -- Coherence vs time --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1
    for cnt_m = indx_list 
        plot(tf,mean(abs(sq(coh_mr(cnt_m,:,10:40))),2)); hold on % mean abs
        plot(tf,abs(mean(sq(coh_mr(cnt_m,:,10:40)),2))); % abs mean
        leg{cnt_lgd} = sprintf('Mean Abs MR coh, ch %d',mod_Ch(cnt)) %%% change this 
        leg{cnt_lgd + 1} = sprintf('Abs Mean MR coh, ch %d',mod_Ch(cnt))
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    legend(leg,'FontSize',10)
    title('M-R Coherence vs time (f avg 10-40 Hz)','FontSize',10);
    grid on
    xlim([250 750])
    hold off
    set(gcf, 'Position',  [100, 600, 800, 400])

    
    fname = strcat(dir_Sess,sprintf('/coherence_vs_time_Sess_%d.jpg',Sess));
    saveas(fig,fname);
    
    
    % -- Coherence vs time HITS/MISSES Mean Abs --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(tf,mean(abs(sq(coh_mr_H(cnt_m,:,10:40))),2)); hold on % mean abs -- HITS
        plot(tf,mean(abs(sq(coh_mr_M(cnt_m,:,10:40))),2));         % mean abs -- MISSES
        leg{cnt_lgd} = sprintf('HITS - Mean Abs MR coh, ch %d',mod_Ch(cnt))
        leg{cnt_lgd + 1} = sprintf('MISSES - Mean Abs MR coh, ch %d',mod_Ch(cnt))
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    legend(leg,'FontSize',10)
    title('HITS/MISSES M-R Mean Abs Coherence vs time (f avg 10-40 Hz)','FontSize',10);
    grid on
    xlim([250 750])
    hold off
    set(gcf, 'Position',  [100, 600, 800, 400])
    
    
    fname = strcat(dir_Sess,sprintf('/coherence_vs_time_Mean_Abs_Hits_Misses_Sess_%d.jpg',Sess));
    saveas(fig,fname);
    
    
    
    % -- Coherence vs time HITS/MISSES Abs Mean --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(tf,abs(mean(sq(coh_mr_H(cnt_m,:,10:40)),2))); hold on % abs mean -- HITS
        plot(tf,abs(mean(sq(coh_mr_M(cnt_m,:,10:40)),2))); % abs mean -- MISSES
        leg{cnt_lgd} = sprintf('HITS - Abs Mean MR coh, ch %d',mod_Ch(cnt))
        leg{cnt_lgd + 1} = sprintf('MISSES - Abs Mean MR coh, ch %d',mod_Ch(cnt))
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    legend(leg,'FontSize',10)
    title('HITS/MISSES M-R Abs Mean Coherence vs time (f avg 10-40 Hz)','FontSize',10);
    grid on
    xlim([250 750])
    hold off
    set(gcf, 'Position',  [100, 600, 800, 400])
    
    
    fname = strcat(dir_Sess,sprintf('/coherence_vs_time_Abs_Mean_Hits_Misses_Sess_%d.jpg',Sess));
    saveas(fig,fname);
    
    
    keyboard 
    % --- FIGURE --------- %%
    % -- Coherence vs frequency --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(f,mean(abs(sq(coh_mr(cnt_m,:,:))),1)); hold on % mean abs
        plot(f,abs(mean(sq(coh_mr(cnt_m,:,:)),1)));  % abs mean
        leg{cnt_lgd} = sprintf('Mean Abs MR coh, ch %d',mod_Ch(cnt))
        leg{cnt_lgd + 1} = sprintf('Abs Mean MR coh, ch %d',mod_Ch(cnt))
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    legend(leg,'FontSize',10)
    title('M-R Coherence vs freq','FontSize',10);
    grid on
    hold off
    set(gcf, 'Position',  [100, 600, 800, 400])

    
    fname = strcat(dir_Sess,sprintf('/coherence_vs_freq_Sess_%d.jpg',Sess));
    saveas(fig,fname);
    
    
    % --- FIGURE --------- %%
    % -- Coherence vs frequency HITS/MISSES Mean Abs --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(f,mean(abs(sq(coh_mr_H(cnt_m,:,:))),1)); hold on % mean abs
        plot(f,mean(abs(sq(coh_mr_M(cnt_m,:,:))),1));  % abs mean
        leg{cnt_lgd} = sprintf('HITS - Mean Abs MR coh, ch %d',mod_Ch(cnt))
        leg{cnt_lgd + 1} = sprintf('MISSES - Mean Abs MR coh, ch %d',mod_Ch(cnt))
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    legend(leg,'FontSize',10)
    title('HITS/MISSES Mean Abs M-R Coherence vs freq','FontSize',10);
    grid on
    hold off
    set(gcf, 'Position',  [100, 600, 1000, 400])

    
    fname = strcat(dir_Sess,sprintf('/coherence_vs_freq_Mean_Abs_Hits_Misses_Sess_%d.jpg',Sess));
    saveas(fig,fname);
    
    
    % --- FIGURE --------- %%
    % -- Coherence vs frequency HITS/MISSES Abs Mean --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(f,abs(mean(sq(coh_mr_H(cnt_m,:,:)),1))); hold on %  abs mean
        plot(f,abs(mean(sq(coh_mr_M(cnt_m,:,:)),1)));  % abs mean 
        leg{cnt_lgd} = sprintf('HITS - Abs Abs MR coh, ch %d',mod_Ch(cnt))
        leg{cnt_lgd + 1} = sprintf('MISSES - Abs Abs MR coh, ch %d',mod_Ch(cnt))
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    legend(leg,'FontSize',10)
    title('HITS/MISSES Abs Mean M-R Coherence vs freq','FontSize',10);
    grid on
    hold off
    set(gcf, 'Position',  [100, 600, 1000, 400])

    
    fname = strcat(dir_Sess,sprintf('/coherence_vs_freq_Abs_Mean_Hits_Misses_Sess_%d.jpg',Sess));
    saveas(fig,fname);
    
    
end

keyboard 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS - MEAN COHERENCES vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- mean coherences across modulator (vs frequency)
mean_cho_MA = mean(sq(mean(abs(coh_mr(:,:,:)),2)));  % mean abs
mean_cho_H_MA = mean(sq(mean(abs(coh_mr_H(:,:,:)),2)));  % HITS mean abs
mean_cho_M_MA = mean(sq(mean(abs(coh_mr_M(:,:,:)),2)));  % MISSES mean abs


std_cho_MA = std(sq(mean(abs(coh_mr(:,:,:)),2)));  % mean abs
std_cho_H_MA = std(sq(mean(abs(coh_mr_H(:,:,:)),2)));  % HITS mean abs
std_cho_M_MA = std(sq(mean(abs(coh_mr_M(:,:,:)),2)));  % MISSES mean abs

% --- Error bars
err_MA = std_cho_MA/48;
err_H_MA = std_cho_H_MA/48;
err_M_MA = std_cho_M_MA/48;

set(0,'DefaultLineLineWidth',2)


fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all
errorbar(f,mean_cho_H_MA,err_H_MA); hold on
errorbar(f,mean_cho_M_MA,err_H_MA); hold on
% errorbar(f,mean_cho_ms_AM,err_ms_AM,'Color',[0.4940, 0.1840, 0.5560]); hold on
% errorbar(f,mean_cho_mr_AM,err_mr_AM,'Color',[102/255, 153/255 0]); hold on
% errorbar(f,mean_cho_sr_AM,err_sr_AM,'color',[26/255 198/255 1]);
grid on
title('STIM: Hits/Misses Mean Abs MR coh of all the caus mod','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('HITS M-R mean abs','MISSES M-R mean abs','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 700, 500])

fname = strcat(dir_Stim,'/mean_coherence_MR_MA_causal_modulators.png');
saveas(fig,fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABS MEAN - MEAN COHERENCES vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- mean coherences across modulator (vs frequency)
mean_cho_AM = mean(sq(abs(mean(coh_mr(:,:,:),2))));  % mean abs
mean_cho_H_AM = mean(sq(abs(mean(coh_mr_H(:,:,:),2))));  % HITS mean abs
mean_cho_M_AM = mean(sq(abs(mean(coh_mr_M(:,:,:),2))));  % MISSES mean abs


std_cho_AM = std(sq(abs(mean(coh_mr(:,:,:),2))));   % mean abs
std_cho_H_AM = std(sq(abs(mean(coh_mr_H(:,:,:),2))));  % HITS mean abs
std_cho_M_AM = std(sq(abs(mean(coh_mr_M(:,:,:),2))));  % MISSES mean abs

% --- Error bars
err_AM = std_cho_AM/48;
err_H_AM = std_cho_H_AM/48;
err_M_AM = std_cho_M_AM/48;

set(0,'DefaultLineLineWidth',2)


fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
hold all
errorbar(f,mean_cho_H_MA,err_H_MA); hold on
errorbar(f,mean_cho_M_MA,err_H_MA); hold on
% errorbar(f,mean_cho_ms_AM,err_ms_AM,'Color',[0.4940, 0.1840, 0.5560]); hold on
% errorbar(f,mean_cho_mr_AM,err_mr_AM,'Color',[102/255, 153/255 0]); hold on
% errorbar(f,mean_cho_sr_AM,err_sr_AM,'color',[26/255 198/255 1]);
grid on
title('STIM: Hits/Misses Abs Mean MR coh of all the caus mod','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('HITS M-R abs mean','MISSES M-R abs mean','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 700, 500])

fname = strcat(dir_Stim,'/mean_coherence_MR_AM_causal_modulators.png');
saveas(fig,fname)


% --- mean coherences across modulator (vs frequency)
mean_Sm = mean(sq(mean(S_m(:,:,:),2)));  % mean abs
mean_Sm_H = mean(sq(mean(S_m_H(:,:,:),2)));  % HITS mean abs
mean_Sm_M = mean(sq(mean(S_m_M(:,:,:),2)));  % MISSES mean abs


std_Sm = std(sq(mean(S_m(:,:,:),2)));  % mean abs
std_Sm_H = std(sq(mean(S_m_H(:,:,:),2)));  % HITS mean abs
std_Sm_M = std(sq(mean(S_m_M(:,:,:),2)));  % MISSES mean abs

% --- Error bars
err_Sm = std_Sm/48;
err_Sm_H = std_Sm_H/48;
err_Sm_M = std_Sm_M/48;



fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all
errorbar(f,mean_Sm_H,err_Sm_H); hold on
errorbar(f,mean_Sm_M,err_Sm_M); hold on
% errorbar(f,mean_cho_ms_AM,err_ms_AM,'Color',[0.4940, 0.1840, 0.5560]); hold on
% errorbar(f,mean_cho_mr_AM,err_mr_AM,'Color',[102/255, 153/255 0]); hold on
% errorbar(f,mean_cho_sr_AM,err_sr_AM,'color',[26/255 198/255 1]);
grid on
title('STIM: Hits/Misses Spectrum caus mod','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('HITS Mod Spectrum','MISSES Mod Spectrum','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 700, 500])

fname = strcat(dir_Stim,'/mean_spectrum_causal_modulators.png');
saveas(fig,fname)
