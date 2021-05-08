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
set(0,'DefaultLineLineWidth',2)

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


% % -- define list of sessions
% if strcmp(monkey,'Maverick')
%     list_sess = 1:19;
%     list_sess(17) = [];
% else
%     list_sess = 1:length(sess_info{1});
% end

% % -- exclude bad sessions 
% excluded_sess = [8,22,30,31];
% excluded_idx = [2,5,8,9];
% sess_list = 1:size(sess_info{1},1);
% sess_list(excluded_idx) = [];

% -- exclude bad sessions --- modified (v2)
excluded_sess = [8,30,31];
excluded_idx = [2,8,9];
sess_list = 1:size(sess_info{1},1);
sess_list(excluded_idx) = [];

% -- print structures on stdout
%format short

for s= sess_list %length(sess_info{1}) %list_sess
    
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
    
    [mod_Ch_rand,area_Ch_rand] = choose_ALL_control_same_Region(RecordPairMRIlabels,MRIlabels,receiver_idx,mod_Ch);

    sess_All_controls_same_area = sess_data;
    sess_All_controls_same_area.ctrl_idx = mod_Ch_rand;
    sess_All_controls_same_area.ctrl_area = area_Ch_rand(:)';
    sess_All_controls_same_area 
    
    dir_Ctrl = strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_same_area',Sess));
    if ~exist(dir_Ctrl, 'dir')
        mkdir(dir_Ctrl)
    end
    save(strcat(dir_Ctrl,'/session_controls_same_area_info_rec001_002_removed_artifacts_V2.mat'),'sess_All_controls_same_area');
   
    
end





