% This routine calculates the RMS of the Pre and Probe pulse signals
%
% Written by Shaoyu Qiao, March 26, 2019

clear all;
close all

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



for iSess = 6 %UsedSess
    %         clearvars -except iSess PreStimSess DATADIR FIGUREDIR MONKEYDIR iSubject subjects UsedSess
    
    disp(['Session ' num2str(iSess) ' out of ' num2str(length(PreStimSess)) ' ...'])
    
    RespPair = sessElectrode(PreStimSess{iSess}); % responding channel
    stimName = PreStimSess{iSess}{9};
    stimTask = PreStimSess{iSess}{7};
    day = sessDay(PreStimSess{iSess});
    
    %% loading Pre data
    dataDir_Pre = sprintf('%s/AccLLR/%sStimAllSess/StimResponseSessions/',DATADIR,stimName);
    switch stimTask
        case 'StimSinglePulse'
            fileName_Pre = sprintf('%sSess%03d_%s_AccLLR_Elec%03d-Elec%03d_%s_1stPulse.mat',dataDir_Pre,iSess,day,RespPair(1),RespPair(2),stimName);
            
        case 'StimBlockShort'
            fileName_Pre = sprintf('%sSess%03d_%s_AccLLR_Elec%03d-Elec%03d_%s_grouped.mat',dataDir_Pre,iSess,day,RespPair(1),RespPair(2),stimName);
    end
    
    tic
    disp('Loading Pre data ...')
    load(fileName_Pre)
    disp('Done with Pre data loading')
    toc
    
    %% identifying beta modulators %%%%%%
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
    
    
    
    %% extract AccLLR results
    Results = Data.AccLLR.Results;
    EventST = Results.EventST;
    %
    %         TargCh = Results.Ch;
    %         LPRawEvent = Results.LPRawEvent;
    %         LPRawNull = Results.LPRawNull;
    %         StimArtifactBlankWin = Results.StimArtifactBlankWin; % ms
    %         StimArtifactStartInd = Results.StimArtifactStartInd;
    %         AccLLRwin = Results.AccLLRwin;
    %         nTr = length(Results.EventST);
    %         AccLLRRawNull = Results.NoHist.Null.AccLLR;
    %         mAccLLRRawNull = mean(AccLLRRawNull);
    %
    %         [dum,ind] = sort(EventST);
    %         nTr_CorrectDetect = sum(~isnan(dum));
    
    %         mLPRawEvent = mean(LPRawEvent,1);
    %         if mLPRawEvent(round(mean(EventST(~isnan(EventST))))) > 0
    %             LPRawEvent = -LPRawEvent;
    %         end
    
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
                    
                    % get modulator channel
                    %                         chs_FreqBand = Data.Spec.ROC.sigChs{iFreqBand};
                    %                         nCh_Use = numel(chs_FreqBand);
                    keyboard
                    %                         nCh_Use = size(Data.RecordPair,1)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       se = size(Data.RecordPair,1)
                    for iCh = 34 %:2
                        
                        display(['Electrode # ......',num2str(iCh)]);
                        
                        %% %%%%%% Spectrogram, permutation test   GINO %%%%%%
                        hitIndx = Data.spec.lfp.DetectedIndx{iCh}; % labels for the hits (which trial was a hit)
                        missIndx = Data.spec.lfp.notDetectedIndx{iCh}; % labels for the misses (which trial was a miss)
                        ModulatorPair = Data.RecordPair(iCh,:); % get the modulator pair
                        ReceiverPair = Data.StimResPairs;
                        
                        ModulatorLFP_Pre = sq(Lfp_Pre(:,ModulatorPair(1),:) - Lfp_Pre(:,ModulatorPair(2),:)); % modulator lfp
                        ReceiverLFP_Pre = sq(Lfp_Pre(:,ReceiverPair(1),:) - Lfp_Pre(:,ReceiverPair(2),:)); % receiver lfp
                        SenderLFP_Pre = sq(Lfp_Pre(:,Data.StimPairs.Syllable(1),:) - Lfp_Pre(:,Data.StimPairs.Syllable(2),:)); % sender lfp
                        
                        % directory path to save files
                        dir_stim = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Stim_data/Sess_%d/Ch_%d',iSess,iCh);
                        
                        if ~exist(dir_stim, 'dir')
                            mkdir(dir_stim)
                        end
                        
                        time = 1000
                        trial = 1;
                        figure;
                        plot(ModulatorLFP_Pre(trial,1:time));
                        % ylim([-200,200])
                        hold on
                        plot(SenderLFP_Pre(trial,1:time));
                        hold on
                        plot(ReceiverLFP_Pre(trial,1:time));
                        legend('Modulator','Sender','Receiver','FontSize',10)
                        title('Lfp Stimulation - baseline','FontSize',13);
                        xlabel('time (ms)','FontSize',13)
                        ylabel('','FontSize',13)
                        grid on
                        hold off
                     
                        
                        %
                        %                         % Complex analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %                         [speck_hit, fhit] = tfsp_proj_GINO(ModulatorLFP_Pre(hitIndx,:), tapers,fs, dn, fk, pad);
                        %                         [speck_miss, fhit] = tfsp_proj_GINO(ModulatorLFP_Pre(missIndx,:), tapers,fs, dn, fk, pad);
                        %
                        %                         % HITS
                        %                         spec_RH = mean(real(speck_hit),3);
                        %                         spec_IH = mean(imag(speck_hit),3);
                        %                         size(spec_R)
                        %
                        %                         figure; tvimage(sq(spec_RH(indx,:,1,:))); title(sprintf('Hit Real - N = %.2f, W = %d, dn = %.4f, k = %d',tapers(1),tapers(2),dn,k)); colorbar;
                        %                         size(spec_RH)
                        %
                        %                         figure; tvimage(sq(spec_IH(indx,:,1,:))); title(sprintf('Hit Im - N = %.2f, W = %d, dn = %.4f, k = %d',tapers(1),tapers(2),dn,k)); colorbar;
                        %                         size(specK_hit)
                        %
                        %                         % MISS
                        %                         spec_RM = mean(real(speck_miss),3);
                        %                         spec_IM = mean(imag(speck_miss),3);
                        %                         size(spec_RM)
                        %
                        %                         figure; tvimage(sq(spec_RM(indx,:,1,:))); title(sprintf('Miss Real - N = %.2f, W = %d, dn = %.4f, k = %d',tapers(1),tapers(2),dn,k)); colorbar;
                        %
                        %                         figure; tvimage(sq(spec_IM(indx,:,1,:))); title(sprintf('Miss Im - N = %.2f, W = %d, dn = %.4f, k = %d',tapers(1),tapers(2),dn,k)); colorbar;
                        %
                        
                        
                        % 1. %%% Standard spectrogram, i.e. homogeneous weighted sum with 1/K weights
                        %                         homogeneous_proj(ModulatorLFP_Pre,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh)
                        
                        % 1. %%% Complex analysis, homogeneous projections
                        homogeneous_proj_complex(ModulatorLFP_Pre,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh)
                        
                        % 2. %%% Projections on unitary vectors, i.e. taper 1, taper 2, ..., taper K
                        %                         orthogonal_proj(ModulatorLFP_Pre,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh)
                        %                         keyboard
                        
                        % 3. %%% Projections on unitary vectors, only on k-1 taper and not k, sum of 2 tapers only
                        %                         orthogonal_proj_combo(ModulatorLFP_Pre,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh)
                        R = 100
                        % 4. %%% Random projections, i.e. inhomogeneous weighted sum
                        %                         random_proj(ModulatorLFP_Pre,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh,R)
                        random_proj_complex(ModulatorLFP_Pre,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh,R)
                        
                        %                         random_proj_subset_tapers(ModulatorLFP_Pre,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh,R)
                        keyboard
                        %                         for freq = 25:30
                        %                             fk = [freq,freq+60];
                        %                             homogeneous_proj_freq(ModulatorLFP_Pre,hitIndx,missIndx,bn_Pre,tapers,fs,dn,fk,pad,nPerm,alpha,iSubject,iSess,iCh)
                        %                         end
                        
                        
                        
                        
                        
                        
                        %                             [SpecTest,ftest] = tfspec(ModulatorLFP_Pre,tapers,fs,dn,fk,pad);
                        %                             figure(1); tvimage(log(sq(SpecTest(1,:,:)))); colorbar;
                        %
                        % Projection of Xk on the K tapers
                        %                             [SpecK, f] = tfsp_proj_GINO(ModulatorLFP_Pre, tapers,fs, dn, fk, pad);
                        %                             % Projection of the spectrogram on the K tapers, i.e. |X_k|^2,
                        %                             Sx_Length = SpecK.*conj(SpecK); % compute the length of the complex number
                        
                        %                             figure(3); tvimage(log(specA(1,:,:))); colorbar; % for trial 1, taper 1
                        %                             figure(3); tvimage(log(mean(specA(:,:,:),1))); colorbar; % for trial 1, taper 1
                        %
                        %                             figure(4); tvimage(log(sq(Sx_Length(1,:,1,:)))); colorbar; % for trial 1, taper 1
                        %                             Sx_Length_mean = sq(mean(Sx_Length,3)); % average across tapers
                        %                             figure(5); tvimage(log(Sx_Length_mean(1,:,:))); colorbar; % spectrogram for trail 1
                        
                        
                        
                    end
                end
            end
        end
    end
    keyboard
    clear Data
    
end
%end