%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates the structures session_control_info.mat containing all the
% info about the data except for the LFPs for the case of the
% the controls electrodes - controls are chosen in the same area as the
% modulator's (all the available electrodes in that area)
%
% INPUT:  session_data_info.mat
% OUTPUT: session_control_info.mat in each Session folder containing a modulator
%
% @ Gino Del Ferraro, December 2020, Pesaran Lab, NYU


clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

% % -- exclude bad sessions
% excluded_sess = [8,22,30,31];
% excluded_idx = [2,5,8,9];
% sess_list = 1:size(sess_info{1},1);
% sess_list(excluded_idx) = [];

% -- exclude bad sessions -- modified (V2)
excluded_sess = [8,30,31];
excluded_idx = [2,8,9];
sess_list = 1:size(sess_info{1},1);
sess_list(excluded_idx) = [];

for s = sess_list
    
    Sess = sess_info{1}(s); % Session number
    dir_Sess = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    
    mod_Ch = sess_data.mod_idx;
    RecordPairMRIlabels = sess_data.RecordPairMRIlabels; % -- MRI labels of the recorder pars
    MRIlabels = sess_data.MRIlabels; % -- all the available MRI labels
    receiver_idx = sess_data.receiver_idx; % -- receiver idx
    
    % --- exclude bad channels
    if Sess == 19
        mod_Ch(mod_Ch == 29) = [];   % Exclude bad channel
    elseif Sess == 29
        mod_Ch(mod_Ch == 68) = [];   % Exclude bad channel
    elseif Sess == 41
        mod_Ch(mod_Ch == 8) = [];   % Exclude bad channel
    end
    
    [mod_Ch_rand,area_Ch_rand] = choose_ALL_control_other_Regions(RecordPairMRIlabels,MRIlabels,receiver_idx,mod_Ch);
    
    sess_All_controls_other_areas = sess_data;
    sess_All_controls_other_areas.ctrl_idx = mod_Ch_rand;
    sess_All_controls_other_areas.ctrl_area = area_Ch_rand(:)';
    sess_All_controls_other_areas
    
    dir_Ctrl = strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_other_areas',Sess));
    if ~exist(dir_Ctrl, 'dir')
        mkdir(dir_Ctrl)
    end
    save(strcat(dir_Ctrl,'/session_controls_other_areas_info_rec001_002_removed_artifacts_V2.mat'),'sess_All_controls_other_areas');
    
    
end





