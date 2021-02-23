function [Lfp_Pre, Raw_Pre, Data] = load_LFP_STIM(PreStimSess,Sess,bn_Pre,DATADIR)

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

end