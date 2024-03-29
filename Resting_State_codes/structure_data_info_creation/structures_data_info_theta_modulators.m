%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates the structures session_data_info.mat containing all the
% info about the data except for the LFPs for the case of the
% modulators' data (not the controls)
%
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
iSubject = 2;
monkey = 'Archie';

if strcmp(subjects{iSubject},'archie')
    archie_vSUBNETS220_rig3
else
    maverick_vSUBNETS220_rig3
end
PreStimSess = PreStimResponseAll_Database_NetworkEdge;

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_RS_Theta = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/%s/Resting_state/theta_band',monkey);
dir_Stim_Theta = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/%s/Stim_data/theta_band',monkey);

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)


% -- print structures on stdout
%format short
for s=1%:size(sess_info{1},1)
    
    Sess = sess_info{1}(s); % Session number
    dir_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d/',Sess));
    load(strcat(dir_Sess,'/Data_with_theta_band.mat')); % --- dataG: all data info and LFP
    
    % -- days and dates
    sess_data.sess_idx = [s Sess];
    sess_data.day = Data.day;
    sess_data.rec_STIM = Data.rec; % -- recording stim 
    sess_data.rec_RS = sess_info{3}{s};
    
    % -- MRI infos 
    sess_data.RecordPair = Data.RecordPair;
    sess_data.MRIlabels = Data.MRIlabels;
    sess_data.RecordPairMRIlabels = Data.RecordPairMRIlabels;
    sess_data.Spec = Data.Spec; % -- p-values and stats 
    
    % -- indexes and pairs
    sess_data.sender_pair = Data.StimPairs.Syllable;
    sess_data.sender_area = PreStimSess{Sess}{11}{1};
    sess_data.receiver_pair = Data.StimResPairs;
    sess_data.receiver_idx = find(Data.RecordPair == Data.StimResPairs(1)); % --- receiver label;
    sess_data.receiver_area = Data.RecordPairMRIlabels(sess_data.receiver_idx,1)';

    sess_data.mod_idx = find(Data.Spec.ROC.sigChIndx{1});
    sess_data.mod_areas = Data.RecordPairMRIlabels(sess_data.mod_idx,1)'
    
    dir_Sess_RS = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    save(strcat(dir_Sess_RS,'/session_data_info.mat'),'sess_data');
   
    % -- print out 
    sess_data
    
end


