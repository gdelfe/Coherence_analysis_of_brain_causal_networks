%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates the structures session_control_info.mat containing all the
% info about the data except for the LFPs for the case of the
% the controls electrodes - controls are chosen in all the brain areas
% except the modulator(s)' and the receiver's
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
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'beta_band';
monkey = 'Maverick';
dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));

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
    receiver_idx = sess_data.receiver_idx; % -- receiver idx 

    [mod_Ch_rand,area_Ch_rand] = choose_ALL_control_other_Regions(RecordPairMRIlabels,MRIlabels,receiver_idx,mod_Ch);

    sess_All_controls_other_areas = sess_data;
    sess_All_controls_other_areas.ctrl_idx = mod_Ch_rand;
    sess_All_controls_other_areas.ctrl_area = area_Ch_rand(:)';
    sess_All_controls_other_areas 
    
    save(strcat(dir_Sess,'/session_all_controls_other_areas_info.mat'),'sess_All_controls_other_areas');
   
    
end





