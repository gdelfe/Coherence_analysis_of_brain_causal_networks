 % This routine calculates the RMS of the Pre and Probe pulse signals
%
% Written by Shaoyu Qiao and Gino Del Ferraro, March, 2021
% 
% This code finds the theta modulator in the STIM experiment across all
% sessions using Shaoyu's method (frequency band in a given range to predict hits and misses)
% 

clear all;
close all

addpath('/mnt/pesaranlab/Matlab/monkeys');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));


% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')

subjects = {'maverick','archie'};

for iSubject = 1% : length(subjects)
    clearvars -except subjects iSubject dir_Stim
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
    
    for iSess = 5 %UsedSess
        clearvars -except iSess PreStimSess DATADIR FIGUREDIR MONKEYDIR iSubject subjects UsedSess dir_Stim
        
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
        
%         tic
%         disp('Loading Pre data ...')
%         load(fileName_Pre)
%         disp('Done with Pre data loading')
%         toc
      
        dir_Sess = strcat(dir_Stim,sprintf('/Sess_%d',iSess));
        load(strcat(dir_Sess,'/Data_with_theta_band.mat'));


        
%         dataB = Data;
%         clear Data
%         dataT.Spec = rmfield(dataT.Spec,'ROC'); % --remove ROC field 
%         
%         
        %% identifying theta modulators %%%%%%
        fs = Data.Fs.lfp;% lfp sampling rate
        Fs = Data.Fs.raw;% raw sampling rate
        
        AnalParams = Data.Params.Anal;
        AnalParams.Tapers = [0.5,2];
        AnalParams.TestSpecDiff.fk = [4 8]; 
        
        fkNames = {'\theta'};
        Data.Params.Anal = AnalParams;
        Data.Spec.ROC.fk = AnalParams.TestSpecDiff.fk;
        StimTrials = Data.StimTrials(Data.goodTrials_index);
        sys = StimTrials(1).MT;
        bn_Pre = [-1000 -5]; % ms
        
        %% extract AccLLR results
        Results = Data.AccLLR.Results;
        EventST = Results.EventST;
        
        TargCh = Results.Ch;
        LPRawEvent = Results.LPRawEvent;
        LPRawNull = Results.LPRawNull;
        StimArtifactBlankWin = Results.StimArtifactBlankWin; % ms
        StimArtifactStartInd = Results.StimArtifactStartInd;
        AccLLRwin = Results.AccLLRwin;
        nTr = length(Results.EventST);
        AccLLRRawNull = Results.NoHist.Null.AccLLR;
        mAccLLRRawNull = mean(AccLLRRawNull);
        
        [dum,ind] = sort(EventST);
        nTr_CorrectDetect = sum(~isnan(dum));
        
        mLPRawEvent = mean(LPRawEvent,1);
        if mLPRawEvent(round(mean(EventST(~isnan(EventST))))) > 0
            LPRawEvent = -LPRawEvent;
        end
        
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
                    
%                     if ~isempty(Data.Spec.ROC.sigChs{iFreqBand})
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
                        
                        %                         get modulator channel
% %                         chs_FreqBand = Data.Spec.ROC.sigChs{iFreqBand};
% %                         nCh_Use = numel(chs_FreqBand);
                        
                        nCh_Use = size(Data.RecordPair,1)
                        for iCh = 52 %1 : nCh_Use
                            
% %                             iCh = find(chs_FreqBand(indx)==Data.RecordPair(:,1)); %-- ask Shaoyu if we should have this line here
                            
%                             fprintf('\n\n');
%                             msg = sprintf('%d/%d Channels',iCh,nCh_Use);
%                             fprintf([reverseStr,msg]);
%                             reverseStr = repmat(sprintf('\b'),1,length(msg));
%                             fprintf('\n\n');

                                                          
%                             figure('Position',[100 100 2000 800])
%                             subplot(2,7,1)
%                             plotElectrodeLocationMRIoverlay(day,TargCh(1));
%                             title(['Receiver e' num2str(TargCh(1))])
%                             
% % %                             subplot(2,7,8)
% % %                             plotElectrodeLocationMRIoverlay(day,chs_FreqBand(indx));
% % %                             title(['Modulator e' num2str(chs_FreqBand(indx))])
% %                         

                            % RECEIVER HITS vs MISSES ------------------
                            h1 = figure;
                            pos1 = get(h1,'Position');
                            ERPplotRange = [-round(max(abs(LPRawEvent(:))),-1) round(max(abs(LPRawEvent(:))),-1)];
                            %ERPplotRange = [-30 30];
                            imagesc((1:size(LPRawEvent,2))/fs*1e3+StimArtifactBlankWin,1:nTr,LPRawEvent(ind,:),ERPplotRange);
                            box off;
                            xlabel('Time from pulse onset (ms)','Fontname','Arial');
                            ylabel('Sortetd events','Fontname','Arial');
                            xlim([0 StimArtifactBlankWin+AccLLRwin])
                            c = colorbar;
                            c.Label.String = 'Voltage (\muV)';
                            
                            % colormap -----------------------
                            a = [0, 17, 250]/255;
                            b = [6, 199, 131]/255;
%                             b = [255, 141, 7]/255;
                            start = 5;
                            maptop = [linspace(1,a(1),64)', linspace(1,a(2),64)', linspace(1,a(3),64)'];
                            mapbot = [linspace(1,b(1),64)', linspace(1,b(2),64)', linspace(1,b(3),64)'];
                            map = [flipud(maptop(start:end,:));mapbot(start:end,:)];
                            colormap(gca,polarmap)
                            colormap(map)
                            % -------------------------------------                            
                            for i = 1 : nTr_CorrectDetect
                                text(dum(i)/fs*1e3+StimArtifactBlankWin,i-0.5,'x','color','k');
                            end
%                             title('Receiver(AccLLR)')
                            
                            hold on
                            h = line([0 size(LPRawEvent,2)/fs*1e3+StimArtifactBlankWin],[nTr_CorrectDetect+0.5 nTr_CorrectDetect+0.5]);
                            set(h,'Linewidth',1,'Linestyle','--','Color','k');
                            
                            text(3,mean(1:nTr_CorrectDetect)+2,'Hit','Rotation',90);
                            text(3,mean(nTr_CorrectDetect+1:nTr)+3,'Miss','Rotation',90);
                            

                            hitIndx = Data.spec.lfp.DetectedIndx{iCh};
                            missIndx = Data.spec.lfp.notDetectedIndx{iCh};
                            
                     
                            
                            %% modulator decoder hit vs miss
                            % FIGURE: plot ROC curve of modulator activity in hit vs miss events
                            
                            stdThresh = Data.spec.lfp.stdThresh;
                            X1 = sq(lfp_Detected(:,iCh,:));
                            [X1,goodInd1] = removeNoisyLfpTrials(X1,stdThresh);
                            
                            X2 = sq(lfp_notDetected(:,iCh,:));
                            [X2,goodInd2] = removeNoisyLfpTrials(X2,stdThresh);
                            fig = figure;
                            [auc,se,S1,S2,roc_Thresh,maxYoudenIndex] = calcRocSpecDiff_HistAUC(X1,X2,AnalParams);                        
                            
                            Data.Spec.ROC.auc{iFreqBand}(iCh) = auc;
                            Data.Spec.ROC.se{iFreqBand}(iCh) = se;    
                            
                            
                            S_all = [];
                            
                            S_all(hitIndx) = S1;
                            S_all(missIndx) = S2;
                            
                            % get different decoding rates across different ROC thresholds
                            [modDecodeHitRate,modDecodeMissRate,rocThresh,~] = runModulatorDecoder(S1,S2,S_all,EventST,roc_Thresh);
                                                        
                            %find optimal decoding rate
                            ind = modDecodeHitRate > modDecodeMissRate;
                            if isequal(sum(ind),0) % modulator hit decoding rates across roc thresholds are all greater than miss decoding rate
                                iThresh = find(roc_Thresh==maxYoudenIndex);
                            elseif strcmp(subjects{iSubject},'archie') && iSess == 13
                                [~,iThresh] = min(abs(roc_Thresh-4.277));
                            else
                                if ~isequal(sum(modDecodeMissRate>0.5),0)
                                    indd = find(ind & modDecodeMissRate>0.5);
                                    modDecodeHitPlusMissRate = modDecodeHitRate(ind & modDecodeMissRate>0.5) + modDecodeMissRate(ind & modDecodeMissRate>0.5);
                                else
                                    indd = find(ind);
                                    modDecodeHitPlusMissRate = modDecodeHitRate(ind) + modDecodeMissRate(ind);
                                end
                                [~,modDecodeHitPlusMissRateInd] = max(modDecodeHitPlusMissRate);
                                [~,iThresh] = min(abs(roc_Thresh-rocThresh(indd(modDecodeHitPlusMissRateInd))));
                            end

                            %iThresh = find(roc_Thresh==maxYoudenIndex);
                            [optModDecodeHitRate,optModDecodeMissRate,optRocThresh,DecoderTrials] = runModulatorDecoder(S1,S2,S_all,EventST,roc_Thresh,iThresh);
                            
                            % FIGURE: Histogram hits vs misses 
                            figure;
                            hg1=histogram(S1,'BinMethod','fd','Normalization','count','FaceColor',[0, 17, 250]/255,'facealpha',0.7);hold on
                            hg2=histogram(S2,'BinMethod','fd','Normalization','count','FaceColor',[6, 199, 131]/255,'facealpha',0.7); %[255, 141, 7]/255
                            plot([optRocThresh,optRocThresh],[0,max([hg1.Values hg2.Values])+1],'b--')
                            legend('Hit','Miss')
                            xlabel(['Mean log ' fkNames{iFreqBand} ' power']);
                            ylabel('Count')
                            title('Modulator')
                            xlim([floor(min([S1;S2])),ceil(max([S1;S2]))]);
                            
                            
                            
                            % FIGURE: Receiver evoked response 
                            figure;
                            avgERP_RecHit_tr = mean(LPRawEvent(hitIndx,:));
                            semERP_RecHit_tr = std(LPRawEvent(hitIndx,:))/sqrt(numel(hitIndx));
                            plotSTA((1:size(LPRawEvent,2))/fs*1e3+StimArtifactBlankWin, avgERP_RecHit_tr,semERP_RecHit_tr,[0, 17, 250]/255);
                            
                            hold on
                            avgERP_RecMiss_tr = mean(LPRawEvent(missIndx,:));
                            semERP_RecMiss_tr = std(LPRawEvent(missIndx,:))/sqrt(numel(missIndx));
%                             plotSTA((1:size(LPRawEvent,2))/fs*1e3+StimArtifactBlankWin,avgERP_RecMiss_tr,semERP_RecMiss_tr,[255, 141, 7]/255);
                            plotSTA((1:size(LPRawEvent,2))/fs*1e3+StimArtifactBlankWin,avgERP_RecMiss_tr,semERP_RecMiss_tr,[6, 199, 131]/255);
                            xlabel('Time after stim onset (ms)');
                            ylabel('Amplitude (\muV)')
                            title('ERP')
                            %legend('Hit',' ', 'Miss',' ','Location','SouthEast')
                            xlim([0 StimArtifactBlankWin+AccLLRwin])
                            %ylim([-20 10])
                            
                            
                            
                            
                        end
                        
                        BrainArea = 'all';
                        saveFigFlag = 0;
                        plotFigFlag = 0;
                        plotElecGrid = 1;
                        plotMRIoverlay = 1;
                        PlotStimElec = 0;
                        %cmap = cbrewer('seq','Blues',6);
                        cmap = parula;
                        
                        display(['Computing p-value for all the electrodes... '])
                        Data = plotModulationNetworkAucMRIoverlay(Data,BrainArea,saveFigFlag,plotFigFlag,plotElecGrid,plotMRIoverlay,PlotStimElec,cmap);
                        
                        if ~isempty(Data.Spec.ROC.sigChs{iFreqBand})
                            
                            dir_Sess = strcat(dir_Stim,sprintf('/Sess_%d',iSess));
                            if ~exist(dir_Sess, 'dir')
                                mkdir(dir_Sess)
                            end
                            
                            save(strcat(dir_Sess,'/Data_with_theta_band.mat'),'Data','-v7.3');
                        end
%                     end
                    
                end
            end
        end
 
        clear Data
    end
end


