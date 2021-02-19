
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
    
    RecordPairMRIlabels = sess_data.RecordPairMRIlabels; % -- MRI labels of the recorder pars 
    MRIlabels = sess_data.MRIlabels; % -- all the available MRI labels 
    send_area = sess_data.sender_area; % -- receiver idx 

    send_ctrl = MRIlabels.(send_area).ElecIndx;  % -- get the indexes of the electrodes in the same brain region as the sender
    
    sess_Send_ctrl_same_area = sess_data;
    sess_Send_ctrl_same_area.send_ctrl_idx = send_ctrl;
    sess_Send_ctrl_same_area.send_ctrl_area = send_area;
    sess_Send_ctrl_same_area 
    
    save(strcat(dir_Sess,'/session_Sender_controls_same_area_info.mat'),'sess_Send_ctrl_same_area');
   
    
end





