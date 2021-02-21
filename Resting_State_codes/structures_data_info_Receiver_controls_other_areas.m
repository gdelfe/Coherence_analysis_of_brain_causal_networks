%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates the structures session_control_info.mat containing all the
% info about the data except for the LFPs for the case of the
% the controls electrodes - controls are for the RECEIVER and are chosen in all the brain areas
% except the receiver's area itself, the sender's, and the modulators' area
%
% INPUT:  session_data_info.mat
% OUTPUT: session_All_control_other_areas_info.mat in each Session folder containing a modulator 
% 
% @ Gino Del Ferraro, December 2020, Pesaran Lab, NYU


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

% -- print structures on stdout
%format short
for s=1:size(sess_info{1},1)
    
    Sess = sess_info{1}(s); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));
    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    
    mod_Ch = sess_data.mod_idx;
    RecordPairMRIlabels = sess_data.RecordPairMRIlabels; % -- MRI labels of the recorder pairs 
    MRIlabels = sess_data.MRIlabels; % -- all the available MRI labels 
    receiver_idx = sess_data.receiver_idx; % -- receiver idx 
    send_area = sess_data.sender_area; 

    % -- same function both for the sender and the receiver's controls all other areas
    [mod_Ch_rand,area_Ch_rand] = choose_Send_Rec_control_other_Regions(RecordPairMRIlabels,MRIlabels,send_area,receiver_idx,mod_Ch,Sess);

    sess_Rec_ctrl_other_areas = sess_data;
    sess_Rec_ctrl_other_areas.ctrl_idx = mod_Ch_rand;
    sess_Rec_ctrl_other_areas.ctrl_area = area_Ch_rand(:)';
    sess_Rec_ctrl_other_areas 
    
    dir_Receiver = strcat(dir_Sess,'/Receiver_controls_other_areas');
    if ~exist(dir_Receiver, 'dir')
        mkdir(dir_Receiver)
    end
    
    
    save(strcat(dir_Receiver,'/session_Receiver_controls_other_areas_info.mat'),'sess_Rec_ctrl_other_areas');
   
    
end





