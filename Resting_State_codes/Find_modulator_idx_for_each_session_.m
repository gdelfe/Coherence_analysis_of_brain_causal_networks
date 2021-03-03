

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
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state'
dlmwrite(strcat(dir_RS,'/Sessions_list_theta.txt'),UsedSess); % write the label list of all the used Sessions 
 
Session = Rest_Database; % Sessions for the RS 

Sess_modulator = {}; % list of Sessions with at least one modulator
Summary = {};
Sess_cnt = 1;

keyboard 

for iSess = UsedSess % For all the sessions 
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
    
    elect_dir = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d',iSess);
    if ~exist(elect_dir, 'dir')
        mkdir(elect_dir)
    end
    
    receiver =  Data.StimResPairs;  % ---- receiver pair
    sender = Data.StimPairs.Syllable; % ---- sender pair
    
    receiver_label = find(Data.RecordPair == receiver(1));
    
    dlmwrite(strcat(elect_dir,sprintf('/receiver_Sess_%d.txt',iSess)),receiver); % write receiver pair
    dlmwrite(strcat(elect_dir,sprintf('/receiver_label_Sess_%d.txt',iSess)),receiver_label); % write receiver pair
    dlmwrite(strcat(elect_dir,sprintf('/sender_Sess_%d.txt',iSess)),sender); % write sender pair
    dlmwrite(strcat(elect_dir,sprintf('/recorded_pairs_modulators_Sess_%d.txt',iSess)),Data.RecordPair);   % ---- all recorded pair
    
    modulator_idx = find(Data.Spec.ROC.sigChIndx{1}); % modulator index 

   
    
    if ~isempty(modulator_idx) % if there is at least one modulator 
        
        for i=1:length(Session)
            if ~isempty(find(strcmp(Session{i}{1},Data.day))) % find the RS session associated to that date (Data.day)
                RS_session = Session{i}{2}{1}; % resting state session 
            end
        end
        
        % --- Session info 
        Sess_modulator{Sess_cnt,1} = iSess; % number of the Session
        Sess_modulator{Sess_cnt,2} = Data.day; % day of the session 
        Sess_modulator{Sess_cnt,3} = RS_session; % corresponding RS session to this stim session
        Sess_cnt = Sess_cnt + 1;
        
        dlmwrite(strcat(elect_dir,'/Modulators_idx.txt'),modulator_idx'); % write the idx of the modulator(s)

    else
        dlmwrite(strcat(elect_dir,'/Modulators_idx.txt'),modulator_idx'); % write empty file
        dlmwrite(strcat(elect_dir,'/NO_Modulator.txt'),modulator_idx'); % write a file to warn there is no modulator
    end
    
    
    % --- Session info
    Summary{iSess,1} = iSess; % number of the Session
    Summary{iSess,2} = Data.day; % day of the session
    Summary{iSess,3} = size(Data.RecordPair,1); % number of tot channels for that Session 
    
    
end
  
dir = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State');
if ~exist(dir, 'dir')
    mkdir(dir)
end

% ---- write the list of Sessions with at least one modulator
writecell(Sess_modulator,strcat(dir,'/Sessions_with_modulator_info.txt'),'delimiter','\t'); % format: label_session, date session, label RS corresponding session 
% ---- Session summary 
writecell(Sess_modulator,strcat(dir,'/Summary_sessions.txt'),'delimiter','\t'); % format: label_session, date session, number of channels 

% [unic,idx_unic] = unique(Sess_modulator(:,2)); % get RS session date without repetitions
% % -- write the RS sessions without repetitions
% writecell(Sess_modulator(idx_unic,:),strcat(dir,'/Sessions_with_mod_no_repetition_info.txt'),'delimiter','\t'); % format: label_session, date session, label RS corresponding session 


    