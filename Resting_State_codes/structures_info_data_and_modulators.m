

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

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

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

% -- load structure files
load(strcat(dir_RS,'/session_AM.mat'))
load(strcat(dir_RS,'/session_MA.mat'))



% -- print structures on stdout
%format short
for s=1:size(sess_info{1},1)
    
    Sess = sess_info{1}(s); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));
    load(strcat(dir_Sess,'/data_RS_and_STIM.mat')); % --- dataG: all data info and LFP
    session_AM(s)
    session_MA(s);
    
    % -- days and dates
    session_data(s).sess_idx = session_AM(s).session_idx;
    session_data(s).day = dataG.day_STIM;
    session_data(s).rec_STIM = dataG.rec_STIM;
    session_data(s).rec_RS = dataG.rec_RS;
    
    % -- MRI infos 
    session_data(s).RecordPair = dataG.RecordPair;
    session_data(s).MRIlabels = dataG.MRIlabels;
    session_data(s).RecordPairMRIlabels = dataG.RecordPairMRIlabels;
    session_data(s).Spec = dataG.Spec; % -- p-values and stats 
    
    % -- indexes and pairs
    session_data(s).sender_pair = session_AM(s).sender;
    session_data(s).sender_area = PreStimSess{Sess}{11}{1};
    session_data(s).receiver_pair = session_AM(s).receiver;
    session_data(s).receiver_idx = dataG.receiver_idx;
    session_data(s).receiver_area = dataG.RecordPairMRIlabels(session_data(s).receiver_idx,1)';
    session_data(s).receiver_area2 = PreStimSess{Sess}{11}{2};

    session_data(s).mod_idx = session_MA(s).mod_idx;
    session_data(s).mod_areas = dataG.RecordPairMRIlabels(session_data(s).mod_idx,1)'
        
    sess_data = session_data(s);
    save(strcat(dir_Sess,'/session_data_info.mat'),'sess_data');
   
    
end


