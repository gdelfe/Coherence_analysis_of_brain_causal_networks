

% This routine calculates the RMS of the Pre and Probe pulse signals
%
% Written by Shaoyu Qiao, March 26, 2019

clear all;
close all

subjects = {'maverick','archie'};

%for
iSubject = 2% : length(subjects) % Loop on the animals
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
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
dir_RS_beta = strcat(dir_main,'/Archie/Resting_state/beta_band')
dlmwrite(strcat(dir_RS_beta,'/Sessions_list.txt'),UsedSess); % write the label list of all the used Sessions

Session = Rest_Database; % Sessions for the RS
Movie_Database

Sess_modulator = {}; % list of Sessions with at least one modulator
Summary = {};
cnt_sess = 1;



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
    
    
    modulator_idx = find(Data.Spec.ROC.sigChIndx{1}); % modulator index
    
    if ~isempty(modulator_idx) % if there is at least one modulator
        
        
        dir_Sess = strcat(dir_RS_beta, sprintf('/Sess_%d',iSess));
        if ~exist(dir_Sess, 'dir')
            mkdir(dir_Sess)
        end
        
        sess_info(cnt_sess).session = [cnt_sess, iSess]; % -- session index
        sess_info(cnt_sess).day = Data.day;
        sess_info(cnt_sess).modulator_idx = modulator_idx;
        
        receiver =  Data.StimResPairs;  % ---- receiver pair
        sender = Data.StimPairs.Syllable; % ---- sender pair
        
        receiver_label = find(Data.RecordPair == receiver(1));
        
        sess_info(cnt_sess).receiver = receiver;
        sess_info(cnt_sess).rec_label = receiver_label;
        sess_info(cnt_sess).sender  = sender;
        sess_info(cnt_sess).pairs  = Data.RecordPair; % -- all recorded pairs
        
        
        dlmwrite(strcat(dir_Sess,sprintf('/receiver_Sess_%d.txt',iSess)),receiver); % write receiver pair
        dlmwrite(strcat(dir_Sess,sprintf('/receiver_label_Sess_%d.txt',iSess)),receiver_label); % write receiver pair
        dlmwrite(strcat(dir_Sess,sprintf('/sender_Sess_%d.txt',iSess)),sender); % write sender pair
        dlmwrite(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',iSess)),Data.RecordPair);   % ---- all recorded pair
        
        
        % --- find the corresponding RS session and date
        RS_check = 0;
        for i=1:length(Session)
            if ~isempty(find(strcmp(Session{i}{1},Data.day))) % find the RS session associated to that date (Data.day)
                RS_session = Session{i}{2}{1}; % resting state session
                RS_check = 1;
                sess_info(cnt_sess).RS = 'YES'
            end
        end
        
        % --- Session info
        Sess_modulator{cnt_sess,1} = iSess; % number of the Session
        Sess_modulator{cnt_sess,2} = Data.day; % day of the session
        if RS_check == 1
            Sess_modulator{cnt_sess,3} = RS_session; % corresponding RS session to this stim session
        else
            Sess_modulator{cnt_sess,3} = '000';
            sess_info(cnt_sess).RS = 'NOT Available'
        end
        
        
        dlmwrite(strcat(dir_Sess,'/Modulators_idx.txt'),modulator_idx'); % write the idx of the modulator(s)
        
        cnt_sess = cnt_sess + 1;
        
        
    else
        display(['--- NO modulator in Session ',num2str(iSess)]);
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     LOAD STIMULATION DATA INTO STRUCTURE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dataG.sess_idx = [cnt_sess, i Sess];
        dataG.day_STIM = Data.day;
        dataG.rec_STIM = Data.rec;
        
        dataG.fs = Data.Fs.lfp;
        dataG.lfp_STIM_Pre = Lfp_Pre;
        
        dataG.MRIlabels = Data.MRIlabels;
        dataG.RecordPair = Data.RecordPair;
        dataG.RecordPairMRIlabels = Data.RecordPairMRIlabels;
        dataG.StimPairs = Data.StimPairs;
        dataG.Spec = Data.Spec;
        
        dataG.receiver =  Data.StimResPairs;  % ---- receiver pair
        dataG.sender = Data.StimPairs.Syllable; % ---- sender pair
        dataG.receiver_idx = find(Data.RecordPair == Data.StimResPairs(1)); % --- receiver label
        dataG.modulators_idx = find(Data.Spec.ROC.sigChIndx{1}); % modulators index
        
        save(strcat(dir_Sess,'/data_STIM.mat'),'dataG');
        
        
        % --- Session info
        Summary{iSess,1} = iSess; % number of the Session
        Summary{iSess,2} = Data.day; % day of the session
        Summary{iSess,3} = size(Data.RecordPair,1); % number of tot channels for that Session
        
    end
    
 
    
end


% ---- write the list of Sessions with at least one modulator
writecell(Sess_modulator,strcat(dir_RS_beta,'/Sessions_with_modulator_info.txt'),'delimiter','\t'); % format: label_session, date session, label RS corresponding session
% ---- Session summary
writecell(Sess_modulator,strcat(dir_RS_beta,'/Summary_sessions.txt'),'delimiter','\t'); % format: label_session, date session, number of channels

% [unic,idx_unic] = unique(Sess_modulator(:,2)); % get RS session date without repetitions
% % -- write the RS sessions without repetitions
% writecell(Sess_modulator(idx_unic,:),strcat(dir,'/Sessions_with_mod_no_repetition_info.txt'),'delimiter','\t'); % format: label_session, date session, label RS corresponding session


for i = 1:11
    sess_info(i)
end



