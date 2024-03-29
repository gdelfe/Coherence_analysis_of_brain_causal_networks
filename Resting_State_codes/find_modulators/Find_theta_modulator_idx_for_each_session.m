

% This routine calculates the RMS of the Pre and Probe pulse signals
%
% Written by Shaoyu Qiao, March 26, 2019

clear all;
close all

subjects = {'maverick','archie'};

monkey = 'Maverick';
iSubject = 1;

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
dir_RS_Theta = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/%s/Resting_state/theta_band',monkey);
dir_Stim_Theta = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/%s/Stim_data/Theta_band',monkey);
sess_list = importdata(strcat(dir_Stim_Theta,'/Sessions_with_modulators_list.txt')); % write the label list of all the used Sessions 
 
Session = Rest_Database; % Sessions for the RS 

Sess_modulator = {}; % list of Sessions with at least one modulator
Summary = {};
cnt_sess = 1;
 

for iSess = 62 %sess_list % For all the sessions 
    %         clearvars -except iSess PreStimSess DATADIR FIGUREDIR MONKEYDIR iSubject subjects UsedSess
    disp(['Session ' num2str(iSess) ' out of ' num2str(length(PreStimSess)) ' ...'])
    
    RespPair = sessElectrode(PreStimSess{iSess}); % responding channel
    stimName = PreStimSess{iSess}{9};
    stimTask = PreStimSess{iSess}{7};
    day = sessDay(PreStimSess{iSess});
    
    % ---- loading Stim Data for theta band 
    dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d/',iSess));
    tic
    disp('Loading Theta band data ...')
    load(strcat(dir_Stim_Sess,'/Data_with_theta_band.mat'));
    disp('Done with Theta band data loading')
    toc

        
    modulator_idx = find(Data.Spec.ROC.sigChIndx{1}); % modulator index 

    if ~isempty(modulator_idx) % if there is at least one modulator 
                
        % -- create directory for the resting state theta band analysis
        dir_Sess = strcat(dir_RS_Theta,sprintf('/Sess_%d',iSess));
        if ~exist(dir_Sess, 'dir')
            mkdir(dir_Sess)
        end
        
        sess_info(cnt_sess).session = [cnt_sess, iSess]; % -- session index  
        sess_info(cnt_sess).modulator_idx = modulator_idx;

        receiver =  Data.StimResPairs;  % ---- receiver pair
        sender = Data.StimPairs.Syllable; % ---- sender pair
        
        receiver_label = find(Data.RecordPair == receiver(1)); % -- receiver label
        
        sess_info(cnt_sess).receiver = receiver;
        sess_info(cnt_sess).rec_label = receiver_label;
        sess_info(cnt_sess).sender  = sender;
        sess_info(cnt_sess).pairs  = Data.RecordPair; % -- all recorded pairs 
        
        
%         dlmwrite(strcat(dir_Sess,sprintf('/receiver_Sess_%d.txt',iSess)),receiver); % write receiver pair
%         dlmwrite(strcat(dir_Sess,sprintf('/receiver_label_Sess_%d.txt',iSess)),receiver_label); % write receiver pair
%         dlmwrite(strcat(dir_Sess,sprintf('/sender_Sess_%d.txt',iSess)),sender); % write sender pair
%         dlmwrite(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',iSess)),Data.RecordPair);   % ---- all recorded pair
%         
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
            Sess_modulator{cnt_sess,3} = 000;
            sess_info(cnt_sess).RS = 'NOT Available'
        end 
            
            
            
%         dlmwrite(strcat(dir_Sess,'/Modulators_idx.txt'),modulator_idx'); % write the idx of the modulator(s)
        
        cnt_sess = cnt_sess + 1;
        
    else 
        display(['--- NO theta modulator in Session ',num2str(iSess)]);
    end
    
    
    % --- Session info
    Summary{iSess,1} = iSess; % number of the Session
    Summary{iSess,2} = Data.day; % day of the session
    Summary{iSess,3} = size(Data.RecordPair,1); % number of tot channels for that Session 
    
    
end
  

keyboard 
% ---- write the list of Sessions with at least one modulator
writecell(Sess_modulator,strcat(dir_RS_Theta,'/Sessions_with_modulator_info.txt'),'delimiter','\t'); % format: label_session, date session, label RS corresponding session 
save(strcat(dir_RS_Theta,'/session_theta_modulator_info.mat'),'sess_info');


load(strcat(dir_RS_Theta,'/session_theta_modulator_info.mat'));

cnt_m = 0;
for i = 1:size(sess_info,2)
    sess_info(i)
    cnt_m = cnt_m + size(sess_info(i).modulator_idx,2)

    if sess_info(i).RS(1) == 'Y'
%        sess_info(i)
        cnt_m = cnt_m + size(sess_info(i).modulator_idx,2)
    end
end 


% [unic,idx_unic] = unique(Sess_modulator(:,2)); % get RS session date without repetitions
% % -- write the RS sessions without repetitions
% writecell(Sess_modulator(idx_unic,:),strcat(dir,'/Sessions_with_mod_no_repetition_info.txt'),'delimiter','\t'); % format: label_session, date session, label RS corresponding session 


    