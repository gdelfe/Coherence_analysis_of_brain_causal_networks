%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates the structures session_data_info.mat containing all the
% info about the data except for the LFPs for the case of the
% modulators' data (not the controls)
% INPUT: structure_AM or structure_MA
% OUTPUT: session_data_info.mat in each Session folder containing a modulator 
% 
% @ Gino Del Ferraro, December 2020, Pesaran Lab, NYU



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
    
    % -- days and dates
    sess_data.sess_idx = session_AM(s).session_idx;
    sess_data.day = dataG.day_STIM;
    sess_data.rec_STIM = dataG.rec_STIM;
    sess_data.rec_RS = dataG.rec_RS;
    
    % -- MRI infos 
    sess_data.RecordPair = dataG.RecordPair;
    sess_data.MRIlabels = dataG.MRIlabels;
    sess_data.RecordPairMRIlabels = dataG.RecordPairMRIlabels;
    sess_data.Spec = dataG.Spec; % -- p-values and stats 
    
    % -- indexes and pairs
    sess_data.sender_pair = session_AM(s).sender;
    sess_data.sender_area = PreStimSess{Sess}{11}{1};
    sess_data.receiver_pair = session_AM(s).receiver;
    sess_data.receiver_idx = dataG.receiver_idx;
    sess_data.receiver_area = dataG.RecordPairMRIlabels(sess_data.receiver_idx,1)';

    sess_data.mod_idx = session_MA(s).mod_idx;
    sess_data.mod_areas = dataG.RecordPairMRIlabels(sess_data.mod_idx,1)'
        
    save(strcat(dir_Sess,'/session_data_info.mat'),'sess_data');
   
    % -- print out 
    session_AM(s)
    sess_data
    
end


