%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% @ Gino Del Ferraro, March 2021, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all

set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - PATHS     --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'beta_band';
monkey = 'Archie'

subjects = {'maverick','archie'};
iSubject = 2;

dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));


%     clearvars -except subjects iSubject
if strcmp(subjects{iSubject},'archie')
    archie_vSUBNETS220_rig3
    Session = Movie_Database; % Sessions for the RS -- for Archie we use the Movie Resting State
else
    maverick_vSUBNETS220_rig3
    Session = Rest_Database; % Sessions for the RS
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


cnt_sess = 1;
for Sess = UsedSess  % For all the sessions

    
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
    
    
   
    
    modulator_idx = find(Data.Spec.ROC.sigChIndx{1}); % modulator index
    
    if ~isempty(modulator_idx) % if there is at least one modulator
        
        display(['-- Session ', num2str(cnt_sess), '  -- label: ',num2str(Sess),',  with modulator(s)'])

        
        dir_Sess = strcat(dir_RS, sprintf('/Sess_%d',Sess));
        if ~exist(dir_Sess, 'dir')
            mkdir(dir_Sess)
        end
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     LOAD STIMULATION DATA INTO STRUCTURE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -- days and dates
        sess_data.sess_idx = [cnt_sess, Sess];
        sess_data.day = Data.day;
        sess_data.rec_STIM = Data.rec;        
        
        % -- MRI infos
        sess_data.RecordPair = Data.RecordPair;
        sess_data.MRIlabels = Data.MRIlabels;
        sess_data.RecordPairMRIlabels = Data.RecordPairMRIlabels;
        sess_data.Spec = Data.Spec; % -- p-values and stats
        
        
        % -- indexes and pairs
        sess_data.sender_pair = Data.StimPairs.Syllable; % ---- sender pair
        sess_data.sender_area = PreStimSess{Sess}{11}{1};
        sess_data.receiver_pair = Data.StimResPairs;  % ---- receiver pair
        sess_data.receiver_idx = find(Data.RecordPair == Data.StimResPairs(1)); % --- receiver label
        sess_data.receiver_area = Data.RecordPairMRIlabels(sess_data.receiver_idx,1)';
        
        sess_data.mod_idx = find(Data.Spec.ROC.sigChIndx{1}); % modulators index
        sess_data.mod_areas = Data.RecordPairMRIlabels(sess_data.mod_idx,1)'
        
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LOAD RESTING STATE DATA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % --- find the corresponding RS session and date
        RS_check = 0;
        for i=1:length(Session)
            if ~isempty(find(strcmp(Session{i}{1},Data.day))) % find the RS session associated to that date (Data.day)
                RS_session = Session{i}{2}{1}; % resting state recording
                RS_check = 1;
                sess_data.RS = 'YES'
                sess_data.rec_RS = RS_session;

            end
        end
        
        % --- Resting state Session info
        Sess_modulator{cnt_sess,1} = Sess; % number of the Session
        Sess_modulator{cnt_sess,2} = Data.day; % day of the session
        if RS_check == 1
            Sess_modulator{cnt_sess,3} = RS_session; % corresponding RS session to this stim session
        else
            Sess_modulator{cnt_sess,3} = '000';
            sess_info.RS = 'NOT Available'
        end
        
        save(strcat(dir_Sess,'/session_data_info.mat'),'sess_data');

       
    end
    cnt_sess = cnt_sess + 1;
    
end


% ---- write the list of Sessions with at least one modulator
writecell(Sess_modulator,strcat(dir_RS,'/Sessions_with_modulator_info.txt'),'delimiter','\t'); % format: label_session, date session, label RS corresponding session






