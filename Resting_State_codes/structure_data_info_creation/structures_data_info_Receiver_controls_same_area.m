%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates the structures session_Receiver_controls_same_area_info 
% containing all the info about the data except for the LFPs for the case of the
% the Receiver controls electrodes chosen in the same brain area as the
% receiver 
%
% INPUT:  session_data_info.mat
% OUTPUT: session_Receiver_controls_same_area_info in each Session 
% 
% @ Gino Del Ferraro, February 2021, Pesaran Lab, NYU


clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

monkey = 'Archie';
dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_RS = fullfile(dir_main, monkey, 'Resting_state', 'theta_band');
fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

% -- load structure files

% -- print structures on stdout
%format short
for s=1:size(sess_info{1},1)
    
    Sess = sess_info{1}(s); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info an
    RecordPairMRIlabels = sess_data.RecordPairMRIlabels; % -- MRI labels of the recorder pars 
    MRIlabels = sess_data.MRIlabels; % -- all the available MRI labels 
    receiver_idx = sess_data.receiver_idx; % -- receiver idx 
    mod_Ch = sess_data.mod_idx; % -- modulator(s)' index
    
    [rec_Ch,area_Ch_rand] = choose_Receiver_control_same_Region(RecordPairMRIlabels,MRIlabels,receiver_idx,mod_Ch);

    sess_Rec_ctrl_same_area = sess_data;
    sess_Rec_ctrl_same_area.ctrl_idx = rec_Ch;
    sess_Rec_ctrl_same_area.ctrl_area = area_Ch_rand(:)';
    sess_Rec_ctrl_same_area 
    
    dir_Receiver = strcat(dir_Sess,'/Receiver_controls_same_area');
    if ~exist(dir_Receiver, 'dir')
        mkdir(dir_Receiver)
    end
    
    save(strcat(dir_Receiver,'/session_Receiver_controls_same_area_info.mat'),'sess_Rec_ctrl_same_area');
   
    
end





