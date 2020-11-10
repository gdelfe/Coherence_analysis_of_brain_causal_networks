
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
close all;

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureVisible','off')

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

addpath('/mnt/pesaranlab/People/Gino/DL-modulators/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State';
dir_Stim = '/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Stim_data';
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


% --------- COHERENCE-GRAM parameters --%
tapers = [0.7 8];
N = tapers(1);
nt = 1000;
dn = 0.005;
fs = 1000;
fk = 200;
nwin = single(floor((nt-N*fs)/(dn*fs)))
k = floor(2*N*tapers(2)-1)
    
cnt_mod = 1; % -- count for the total number of modulators 

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
    
    
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',   out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    
    % ---- bipolar referencing, pairs of electrodes
    dir_Stim_Sess = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Stim_data/Sess_%d',Sess);
    dir_RS_Sess = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d',Sess);
    
    if ~exist(dir_Stim_Sess, 'dir')
        mkdir(dir_Stim_Sess)
    end
    
    % -- load list electrodes, sender, receiver
    electrode = importdata(strcat(dir_RS_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairs
    receiver = importdata(strcat(dir_RS_Sess,sprintf('/receiver_Sess_%d.txt',Sess)));  % ---- receiver pair
    sender = importdata(strcat(dir_RS_Sess,sprintf('/sender_Sess_%d.txt',Sess))); % ---- sender pair
    
    
    
    %% %%%%%% Assign LFP %%%%%%
    lfp_E = sq(Lfp_Pre(:,electrode(:,1),:) - Lfp_Pre(:,electrode(:,2),:)); % modulator lfp
    lfp_R = sq(Lfp_Pre(:,receiver(1),:) - Lfp_Pre(:,receiver(2),:)); % receiver lfp
    lfp_S = sq(Lfp_Pre(:,sender(1),:) - Lfp_Pre(:,sender(2),:)); % sender lfp
    
    % --- coherence
    
    mod_Ch = session_AM(i).mod_idx; % causal modulators channel
    
    indx_list = [];
    
    for Ch = mod_Ch % for all the modulators in the session
        
        close all
        
        indx_list = [indx_list, cnt_mod]; % store the cnt number --- needed for multiple plotting 
                
        % -- labels Hits and Misses
        hitIndx = Data.spec.lfp.DetectedIndx{Ch}; % labels for the hits (which trial was a hit)
        missIndx = Data.spec.lfp.notDetectedIndx{Ch}; % labels for the misses (which trial was a miss)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --------- COHERENCE-GRAM ----------- --%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        display(['Computing coherence-gram...'])
        

        % Compute coherences and spectrums, all, hits, and misses 
        [c_mr,tf,f,spec_m,spec_r] = tfcoh_GINO(sq(lfp_E(:,Ch,:)),lfp_R(:,:),tapers,1e3,dn,fk,2,0.05,11); % coherence modulator-receiver
        [c_mr_H,tf,f,spec_m_H,spec_r_H] = tfcoh_GINO(sq(lfp_E(hitIndx,Ch,:)),lfp_R(hitIndx,:),tapers,1e3,dn,fk,2,0.05,11); % Hits
        [c_mr_M,tf,f,spec_m_M,spec_r_M] = tfcoh_GINO(sq(lfp_E(missIndx,Ch,:)),lfp_R(missIndx,:),tapers,1e3,dn,fk,2,0.05,11); % Misses 
        
        % -- store coherence 
        coh_mr(cnt_mod).all = c_mr;
        coh_mr(cnt_mod).hits = c_mr_H;
        coh_mr(cnt_mod).misses = c_mr_M;
        coh_mr(cnt_mod).HitsLabels = hitIndx;
        coh_mr(cnt_mod).MissesLabels = missIndx;
        coh_mr(cnt_mod).tapers = tapers;
        coh_mr(cnt_mod).dn = dn;
        
        % -- store spectrum modulator 
        S_m(cnt_mod).all = spec_m;
        S_m(cnt_mod).hits = spec_m_H;
        S_m(cnt_mod).misses = spec_m_M;
        S_m(cnt_mod).HitsLabels = hitIndx;
        S_m(cnt_mod).MissesLabels = missIndx;
        S_m(cnt_mod).tapers = tapers;
        S_m(cnt_mod).dn = dn;
        
        % -- store spectrum receiver  
        S_r(cnt_mod).all = spec_r;
        S_r(cnt_mod).hits = spec_r_H;
        S_r(cnt_mod).misses = spec_r_M;
        S_r(cnt_mod).HitsLabels = hitIndx;
        S_r(cnt_mod).MissesLabels = missIndx;
        S_r(cnt_mod).tapers = tapers;
        S_r(cnt_mod).dn = dn;
        
         
        %     dlmwrite(strcat(dir_Sess,'/sr_coherogram.txt'),c_sr,'delimiter',' ');
        % -- Figure: coherence spectrum
        
        
        % --- FIGURE: MODULATOR-RECEIVER COHEROGRAM
        fig_mr = figure; tvimage(abs(c_mr(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:50:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('M-R Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Stim_Sess,sprintf('/MR_coherogram_ch_%d_N_%.2f_W_%d_fk_%d.jpg',Ch,tapers(1),tapers(2),fk));
        saveas(fig_mr,fname);
        

        % --- FIGURE: HITS MODULATOR-RECEIVER COHEROGRAM
        fig_mr = figure; tvimage(abs(c_mr_H(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:50:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('HITS M-R Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Stim_Sess,sprintf('/MR_coherogram_HITS_ch_%d_N_%.2f_W_%d_fk_%d.jpg',Ch,tapers(1),tapers(2),fk));
        saveas(fig_mr,fname);
        
        
        % --- FIGURE: MISSES MODULATOR-RECEIVER COHEROGRAM
        fig_mr = figure; tvimage(abs(c_mr_M(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:50:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('MISSES M-R Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Stim_Sess,sprintf('/MR_coherogram_MISSES_ch_%d_N_%.2f_W_%d_fk_%d.jpg',Ch,tapers(1),tapers(2),fk));
        saveas(fig_mr,fname);
        
        cnt_mod = cnt_mod + 1;
        
    end % --- end of all modulator channels
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MODULATOR - RECEIVER COHERENCE FIGURES   %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fmin = 10;
    fmax = 40;
    % --- FIGURE --------- %%
    % -- Coherence vs time --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1
    for cnt_m = indx_list 
        plot(tf,mean(abs(sq(coh_mr(cnt_m).all(:,fmin:fmax))),2)); hold on % mean abs
        plot(tf,abs(mean(sq(coh_mr(cnt_m).all(:,fmin:fmax)),2))); % abs mean
        leg{cnt_lgd} = sprintf('Mean Abs MR coh, ch %d',mod_Ch(cnt)); %%% change this 
        leg{cnt_lgd + 1} = sprintf('Abs Mean MR coh, ch %d',mod_Ch(cnt));
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    xlabel('time (sec)');
    ylabel('coherence')
    legend(leg,'FontSize',10,'Location','southoutside')
    title('M-R Coherence vs time (f avg 10-40 Hz)','FontSize',10);
    grid on
    xlim([320 680])
    hold off
    set(gcf, 'Position',  [100, 600, 800, 400])

    
    fname = strcat(dir_Stim_Sess,sprintf('/coherence_vs_time_Sess_%d_N_%.2f_W_%d_fk_%d.jpg',Sess,tapers(1),tapers(2),fk));
    saveas(fig,fname);
    
    
    % -- Coherence vs time HITS/MISSES Mean Abs --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(tf,mean(abs(sq(coh_mr(cnt_m).hits(:,fmin:fmax))),2)); hold on % mean abs -- HITS
        plot(tf,mean(abs(sq(coh_mr(cnt_m).misses(:,fmin:fmax))),2));         % mean abs -- MISSES
        leg{cnt_lgd} = sprintf('HITS - Mean Abs MR coh, ch %d',mod_Ch(cnt));
        leg{cnt_lgd + 1} = sprintf('MISSES - Mean Abs MR coh, ch %d',mod_Ch(cnt));
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    xlabel('time (sec)');
    ylabel('coherence')
    legend(leg,'FontSize',10,'Location','southoutside')
    title('HITS/MISSES M-R Mean Abs Coherence vs time (f avg 10-40 Hz)','FontSize',10);
    grid on
    xlim([250 750])
    hold off
    set(gcf, 'Position',  [100, 600, 800, 400])
    
    
    fname = strcat(dir_Stim_Sess,sprintf('/coherence_vs_time_MA_Hit_Miss_Sess_%d_N_%.2f_W_%d_fk_%d.jpg',Sess,tapers(1),tapers(2),fk));
    saveas(fig,fname);
    
    
    
    % -- Coherence vs time HITS/MISSES Abs Mean --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(tf,abs(mean(sq(coh_mr(cnt_m).hits(:,fmin:fmax)),2))); hold on % abs mean -- HITS
        plot(tf,abs(mean(sq(coh_mr(cnt_m).misses(:,fmin:fmax)),2))); % abs mean -- MISSES
        leg{cnt_lgd} = sprintf('HITS - Abs Mean MR coh, ch %d',mod_Ch(cnt));
        leg{cnt_lgd + 1} = sprintf('MISSES - Abs Mean MR coh, ch %d',mod_Ch(cnt));
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    xlabel('time (sec)');
    ylabel('coherence');
    legend(leg,'FontSize',10,'Location','southoutside')
    title('HITS/MISSES M-R Abs Mean Coherence vs time (f avg 10-40 Hz)','FontSize',10);
    grid on
    xlim([250 750])
    hold off
    set(gcf, 'Position',  [100, 600, 800, 400])
    
    
    fname = strcat(dir_Stim_Sess,sprintf('/coherence_vs_time_AM_Hit_Miss_Sess_%d_N_%.2f_W_%d_fk_%d.jpg',Sess,tapers(1),tapers(2),fk));
    saveas(fig,fname);
    
     
    % --- FIGURE --------- %%
    % -- Coherence vs frequency --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(f,mean(abs(sq(coh_mr(cnt_m).all)),1)); hold on % mean abs
        plot(f,abs(mean(sq(coh_mr(cnt_m).all),1)));  % abs mean
        leg{cnt_lgd} = sprintf('Mean Abs MR coh, ch %d',mod_Ch(cnt));
        leg{cnt_lgd + 1} = sprintf('Abs Mean MR coh, ch %d',mod_Ch(cnt));
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    xlabel('freq (Hz)');
    ylabel('coherence');
    legend(leg,'FontSize',10,'Location','southoutside')
    title('M-R Coherence vs freq','FontSize',10);
    grid on
    hold off
    set(gcf, 'Position',  [100, 600, 800, 400])

    
    fname = strcat(dir_Stim_Sess,sprintf('/coherence_vs_freq_Sess_%d_N_%.2f_W_%d_fk_%d.jpg',Sess,tapers(1),tapers(2),fk));
    saveas(fig,fname);
    
    
    % --- FIGURE --------- %%
    % -- Coherence vs frequency HITS/MISSES Mean Abs --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(f,mean(abs(sq(coh_mr(cnt_m).hits)),1)); hold on % mean abs
        plot(f,mean(abs(sq(coh_mr(cnt_m).misses)),1));  % abs mean
        leg{cnt_lgd} = sprintf('HITS - Mean Abs MR coh, ch %d',mod_Ch(cnt));
        leg{cnt_lgd + 1} = sprintf('MISSES - Mean Abs MR coh, ch %d',mod_Ch(cnt));
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    xlabel('freq (Hz)');
    ylabel('coherence');
    legend(leg,'FontSize',10,'Location','southoutside')
    legend(leg,'FontSize',10)
    title('HITS/MISSES Mean Abs M-R Coherence vs freq','FontSize',10);
    grid on
    hold off
    set(gcf, 'Position',  [100, 600, 1000, 400])

    
    fname = strcat(dir_Stim_Sess,sprintf('/coherence_vs_freq_MA_Hit_Miss_Sess_%d_N_%.2f_W_%d_fk_%d.jpg',Sess,tapers(1),tapers(2),fk));
    saveas(fig,fname);
    
    
    % --- FIGURE --------- %%
    % -- Coherence vs frequency HITS/MISSES Abs Mean --- %
    leg = cell(2*size(mod_Ch,2),1); % dynamic legend
    fig = figure;
    cnt_lgd = 1;
    cnt = 1;
    for cnt_m = indx_list 
        plot(f,abs(mean(sq(coh_mr(cnt_m).hits),1))); hold on %  abs mean
        plot(f,abs(mean(sq(coh_mr(cnt_m).misses),1)));  % abs mean 
        leg{cnt_lgd} = sprintf('HITS - Abs Abs MR coh, ch %d',mod_Ch(cnt));
        leg{cnt_lgd + 1} = sprintf('MISSES - Abs Abs MR coh, ch %d',mod_Ch(cnt));
        cnt_lgd = cnt_lgd + 2;
        cnt = cnt + 1;
    end
    xlabel('freq (Hz)');
    ylabel('coherence');
    legend(leg,'FontSize',10,'Location','southoutside')
    title('HITS/MISSES Abs Mean M-R Coherence vs freq','FontSize',10);
    grid on
    hold off
    set(gcf, 'Position',  [100, 600, 1000, 400])

    
    fname = strcat(dir_Stim_Sess,sprintf('/coherence_vs_freq_Abs_Mean_Hit_Miss_Sess_%d_N_%.2f_W_%d_fk_%d.jpg',Sess,tapers(1),tapers(2),fk));
    saveas(fig,fname);
    
    
end


% -- Write structures files
save(strcat(dir_Stim,sprintf('/coherence_N_%.2f_W_%d_fk_%d.mat',tapers(1),tapers(2),fk)),'coh_mr');
save(strcat(dir_Stim,sprintf('/spec_m_N_%.2f_W_%d_fk_%d.mat',tapers(1),tapers(2),fk)),'S_m');
save(strcat(dir_Stim,sprintf('/spec_r_N_%.2f_W_%d_fk_%d.mat',tapers(1),tapers(2),fk)),'S_r');

keyboard 

nmod = 48; % -- number of modulators 
t_tot = size(coh_mr(1).all,1); 
ftot = size(coh_mr(1).all,2);

% -- create matrices to perform the averages 
coh_MR = zeros(nmod,t_tot,ftot);
coh_MR_H = zeros(nmod,t_tot,ftot);
coh_MR_M = zeros(nmod,t_tot,ftot);
spec_M = zeros(nmod,t_tot,ftot);
spec_M_H = zeros(nmod,t_tot,ftot);
spec_M_M = zeros(nmod,t_tot,ftot);
spec_R = zeros(nmod,t_tot,ftot);
spec_R_H = zeros(nmod,t_tot,ftot);
spec_R_M = zeros(nmod,t_tot,ftot);

% -- assign structure values to matrices 
for mod = 1:nmod

    coh_MR(mod,:,:) = coh_mr(mod).all;
    coh_MR_H(mod,:,:) = coh_mr(mod).hits;
    coh_MR_M(mod,:,:) = coh_mr(mod).misses;
    
    spec_M(mod,:,:) = S_m(mod).all;
    spec_M_H(mod,:,:) = S_m(mod).hits;
    spec_M_M(mod,:,:) = S_m(mod).misses;
    
    spec_R(mod,:,:) = S_r(mod).all;
    spec_R_H(mod,:,:) = S_r(mod).hits;
    spec_R_M(mod,:,:) = S_r(mod).misses;

end


% -- load structures files 
% load_coh = load(strcat(dir_Stim,sprintf('/coherence_N_%.2f_W_%d_fk_%d.mat',tapers(1),tapers(2),fk)));
% coh_mr = load_coh.coh_mr;
% 
% load_spec_m = load(strcat(dir_Stim,sprintf('/spec_m_N_%.2f_W_%d_fk_%d.mat',tapers(1),tapers(2),fk)));
% spec_m = load_spec_m.spec_m;
% 
% load_spec_r = load(strcat(dir_Stim,sprintf('/spec_r_N_%.2f_W_%d_fk_%d.mat',tapers(1),tapers(2),fk)));
% spec_r = load_spec_r.spec_r;


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS - MEAN COHERENCES vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tin = 1;
tfin = 60;
% --- mean coherences across modulator (vs frequency)
mean_cho_MA = mean(sq(mean(abs(coh_MR),2)));  % mean abs
mean_cho_H_MA = mean(sq(mean(abs(coh_MR_H(:,tin:tfin,:)),2)));  % HITS mean abs
mean_cho_M_MA = mean(sq(mean(abs(coh_MR_M(:,tin:tfin,:)),2)));  % MISSES mean abs


std_cho_MA = std(sq(mean(abs(coh_MR),2)));  % mean abs
std_cho_H_MA = std(sq(mean(abs(coh_MR_H(:,tin:tfin,:)),2)));  % HITS mean abs
std_cho_M_MA = std(sq(mean(abs(coh_MR_M(:,tin:tfin,:)),2)));  % MISSES mean abs

% --- Error bars
err_MA = std_cho_MA/sqrt(48);
err_H_MA = std_cho_H_MA/sqrt(48);
err_M_MA = std_cho_M_MA/sqrt(48);

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureVisible','on')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- FIGURE -- Mean Abs/Abs Mean - %

fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='linear'
hold all
shadedErrorBar(f,mean_cho_MA,err_MA,'lineprops',{'color',[170, 57, 57]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_AM,err_AM,'lineprops',{'color',[85, 170, 85]/255},'patchSaturation',0.4); hold on
grid on
title('STIM: Mean Abs/Abs Mean MR coh of all the caus mod','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('M-R mean abs','M-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 800, 500])

fname = strcat(dir_Stim,sprintf('/mean_coherence_MR_MA_AM_causal_modulators_N_%.1f_W_%d_dn_%.3f.png',tapers(1),tapers(2),dn));
saveas(fig,fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- FIGURE - Hits/Misses   %

fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='linear'
hold all
shadedErrorBar(f,mean_cho_H_MA,err_H_MA,'lineprops',{'color',[0, 102, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_M_MA,err_M_MA,'lineprops',{'color',[255, 102, 0]/255},'patchSaturation',0.4); hold on
grid on
title('STIM: Hits/Misses Mean Abs MR coh of all the caus mod','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('HITS -- M-R mean abs','MISSES -- M-R mean abs','FontSize',10)
set(gcf, 'Position',  [100, 600, 800, 500])

fname = strcat(dir_Stim,sprintf('/mean_coherence_MR_MA_causal_modulators_N_%.1f_W_%d_dn_%.3f.png',tapers(1),tapers(2),dn));
saveas(fig,fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABS MEAN - MEAN COHERENCES vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- mean coherences across modulator (vs frequency)
mean_cho_AM = mean(sq(abs(mean(coh_MR,2))));  % abs mean
mean_cho_H_AM = mean(sq(abs(mean(coh_MR_H,2))));  % HITS abs mean 
mean_cho_M_AM = mean(sq(abs(mean(coh_MR_M,2)))); % MISSES abs mean 


std_cho_AM = std(sq(abs(mean(coh_MR,2))));  % abs mean
std_cho_H_AM = std(sq(abs(mean(coh_MR_H,2))));  % HITS abs mean 
std_cho_M_AM = std(sq(abs(mean(coh_MR_M,2)))); % MISSES abs mean 

% --- Error bars
err_AM = std_cho_AM/sqrt(48);
err_H_AM = std_cho_H_AM/sqrt(48);
err_M_AM = std_cho_M_AM/sqrt(48);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- FIGURE -- Hits/Misses  %

fig = figure;
hold all

shadedErrorBar(f,mean_cho_H_MA,err_H_AM,'lineprops',{'color',[0, 102, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_M_MA,err_H_AM,'lineprops',{'color',[255, 102, 0]/255},'patchSaturation',0.4); hold on
grid on
title('STIM: Hits/Misses Abs Mean MR coh of all the caus mod','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('HITS -- M-R abs mean','MISSES -- M-R abs mean','FontSize',10)
set(gcf, 'Position',  [100, 600, 800, 500])

fname = strcat(dir_Stim,sprintf('/mean_coherence_MR_AM_causal_modulators_N_%.1f_W_%d_dn_%.3f.png',tapers(1),tapers(2),dn));
saveas(fig,fname)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----- SPECTRUM - HITS AND MISSES ------  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  %%%%%%%%%%%% MODULATORS

% --- mean spectrum modulator
mean_spec_m = mean(sq(mean(spec_M,2))); 
mean_spec_m_H = mean(sq(mean(spec_M_H,2)));  % HITS
mean_spec_m_M = mean(sq(mean(spec_M_M,2))); % MISSES 

% --- std 
std_Sm = std(sq(mean(spec_M,2)));
std_Sm_H = std(sq(mean(spec_M_H,2)));  % HITS
std_Sm_M = std(sq(mean(spec_M_M,2))); % MISSES 

% --- Error bars
err_m_AM = std_Sm/sqrt(48);
err_m_H_AM = std_Sm_H/sqrt(48);
err_m_M_AM = std_Sm_M/sqrt(48);


fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'
hold all
shadedErrorBar(f,mean_spec_m_H,err_m_H_AM,'lineprops',{'color',[0, 102, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_spec_m_M,err_m_M_AM,'lineprops',{'color',[255, 102, 0]/255},'patchSaturation',0.4); hold on
grid on
title('STIM: Hits/Misses Spectrum caus mod','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
legend('HITS --- Mod Spectrum','MISSES -- Mod Spectrum','FontSize',10)
set(gcf, 'Position',  [100, 600, 900, 500])

fname = strcat(dir_Stim,sprintf('/mean_spectrum_causal_modulators_N_%.1f_W_%d_dn_%.3f.png',tapers(1),tapers(2),dn));
saveas(fig,fname)



% %% RECEIVER %%%%%%%%%%%%%%%%%%%%%%%%%

% --- mean spectrum receiver
mean_spec_r = mean(sq(mean(spec_R,2))); 
mean_spec_r_H = mean(sq(mean(spec_R_H,2)));  % HITS
mean_spec_r_M = mean(sq(mean(spec_R_M,2))); % MISSES 

% --- std 
std_Sr = std(sq(mean(spec_R,2)));
std_Sr_H = std(sq(mean(spec_R_H,2)));  % HITS
std_Sr_M = std(sq(mean(spec_R_M,2))); % MISSES 

% --- Error bars
err_r_AM = std_Sr/sqrt(48);
err_r_H_AM = std_Sr_H/sqrt(48);
err_r_M_AM = std_Sr_M/sqrt(48);


fig = figure;
hAx=axes;
hAx.XScale='linear'
hAx.YScale='log'

shadedErrorBar(f,mean_spec_r_H,err_r_H_AM,'lineprops',{'color',[0, 102, 255]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_spec_r_M,err_r_M_AM,'lineprops',{'color',[255, 102, 0]/255},'patchSaturation',0.4); hold on
grid on
title('STIM: Hits/Misses Spectrum receiver','FontSize',11);
xlabel('freq (Hz)');
ylabel('spectrum');
legend('HITS --- Mod Spectrum','MISSES -- Mod Spectrum','FontSize',10)
set(gcf, 'Position',  [100, 600, 900, 500])

fname = strcat(dir_Stim,sprintf('/mean_spectrum_receivers_N_%.1f_W_%d_dn_%.3f.png',tapers(1),tapers(2),dn));
saveas(fig,fname)

