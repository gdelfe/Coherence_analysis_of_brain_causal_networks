

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code loads RS and STIM data (lfp and all the relevant infos about
% electrodes, p-values, etc...)
%
%
% @ Gino Del Ferraro, December 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all

set(0,'DefaultLineLineWidth',2)

subjects = {'maverick','archie'};

%for
iSubject = 1 % : length(subjects) % Loop on the animals
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
dir_Stim = '//mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';


fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);



for i=1:size(sess_info{1},1)  % For all the session with a modulator
    
    clear dataG Data
    Sess = sess_info{1}(i); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));

    disp(['STIMULATION: Session ' num2str(Sess) ' out of ' num2str(length(PreStimSess)) ' ...'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD DATA STIM DATA ...                %
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
    
    % %%%%%%%%%%%%%%%%% --->>> STIMULATION LFP data loaded at this point 
    
    
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',   out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     LOAD STIMULATION DATA INTO STRUCTURE 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    dataG.sess_idx = [i, Sess];
    dataG.MRIlabels = Data.MRIlabels;
    dataG.lfp_STIM_Pre = Lfp_Pre;
    dataG.fs = Data.Fs.lfp;
    dataG.day_STIM = Data.day;
    dataG.rec_STIM = Data.rec;
    dataG.RecordPair = Data.RecordPair;
    dataG.RecordPairMRIlabels = Data.RecordPairMRIlabels;
    dataG.StimPairs = Data.StimPairs;
    dataG.Spec = Data.Spec;

    dataG.receiver =  Data.StimResPairs;  % ---- receiver pair
    dataG.sender = Data.StimPairs.Syllable; % ---- sender pair
    dataG.receiver_idx = find(Data.RecordPair == Data.StimResPairs(1)); % --- receiver label 
    dataG.modulators_idx = find(Data.Spec.ROC.sigChIndx{1}); % modulators index 

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD RESTING STATE DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % addpath('/vol/sas8/Maverick_RecStim_vSUBNETS220/160125/004')
    addpath(sprintf('/vol/sas8/Maverick_RecStim_vSUBNETS220/%s/%s/',sess_info{2}{i},sess_info{3}{i})) % add path of the specific RS session
    
    % file = 'rec004.Frontal.lfp.dat'
    file = sprintf('rec%s.Frontal.lfp.dat',sess_info{3}{i})
    fid = fopen(file);
    format = 'float=>single';
    
    CH = 220; % tot number of channels
    FS = 1000; % sampling
    
    Sess = sess_info{1}(i); % Session number
    display(['-- RESTING STATE: Session ',num2str(i),' -- label: ',num2str(Sess),',   out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    dataRS = fread(fid,[CH,inf],format); % load the RS data
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     LOAD RESTING STATE DATA INTO STRUCTURE 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    dataG.day_RS = sess_info{2}{i}; % -- date of RS acquisition 
    dataG.rec_RS = sess_info{3}{i};
    dataG.lfpRS = dataRS; % -- all the 220 recorded electrodes
    
    
    save(strcat(dir_Sess,'/data_RS_and_STIM.mat'),'dataG');

    
end 
    