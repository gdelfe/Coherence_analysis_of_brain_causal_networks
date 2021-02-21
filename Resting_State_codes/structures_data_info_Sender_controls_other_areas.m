%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates the structures session_control_info.mat containing all the
% info about the data except for the LFPs for the case of the
% the controls electrodes - controls are for the RECEIVER and are chosen in all the brain areas
% except the receiver's area itself and the modulators'
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
    RecordPairMRIlabels = sess_data.RecordPairMRIlabels; % -- MRI labels of the recorder pars 
    MRIlabels = sess_data.MRIlabels; % -- all the available MRI labels 
    send_area = sess_data.sender_area; % -- sender's brain area
     
    receiver_idx = sess_data.receiver_idx; % -- receiver idx 

    % -- NOTE: For the receivers controls, you can use the same function as
    % the modulators' controls, since you want to esclude also the
    % modulators' brain regions from the "controls other areas"
    [mod_Ch_rand,area_Ch_rand] = choose_Send_Rec_control_other_Regions(RecordPairMRIlabels,MRIlabels,send_area,receiver_idx,mod_Ch,Sess);

    sess_Send_ctrl_other_areas = sess_data;
    sess_Send_ctrl_other_areas.ctrl_idx = mod_Ch_rand;
    sess_Send_ctrl_other_areas.ctrl_area = area_Ch_rand(:)';
    sess_Send_ctrl_other_areas 
    
    dir_Sender = strcat(dir_Sess,'/Sender_controls_other_areas');
    if ~exist(dir_Sender, 'dir')
        mkdir(dir_Sender)
    end
    
    
    save(strcat(dir_Sender,'/session_Sender_controls_other_areas_info.mat'),'sess_Send_ctrl_other_areas');
   
    
end





