

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

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
    
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));
    load(strcat(dir_Sess,'/data_RS_and_STIM.mat')); % --- dataG: all data info and LFP
    session_AM(s)
    session_MA(s)
    
    session_data(s).day = dataG.day_STIM;
    session_data(s).rec_STIM = dataG.rec_STIM;
    session_data(s).rec_RS = dataG.rec_RS;
    session_data(s).RecordPair = dataG.RecordPair;
    session_data(s).MRIlabels = dataG.MRIlabels;
    session_data(s).RecordPairMRIlabels = dataG.RecordPairMRIlabels;
    session_data(s).Spec = dataG.Spec;
    session_data(s).sess_idx = session_AM(s).session_idx;
    session_data(s).sender_pair = session_AM(s).sender;
    session_data(s).receiver_pair = session_AM(s).receiver;
    session_data(s).receiver_idx = dataG.receiver_idx;
    session_data(s).mod_idx = session_MA(s).mod_idx;
    
end


