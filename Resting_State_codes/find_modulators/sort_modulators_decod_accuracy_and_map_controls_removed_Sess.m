%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code sorts the modulators according to their decoding accuracy
% For all the controls, same area and different area, it lists them
% together with the session label and the i indx
% 
% @ Gino Del Ferraro, NYU, March 2021



clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
monkey = 'Maverick';

dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_out = strcat(dir_main,sprintf('%s/Resting_state/%s/Modulators_controls',monkey,freq_band));

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); %: USE MOVIE also for Maverick (bad sess removed) load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

filename = '_AUC'; % -- write out file for the modulators 
filename_ctrl = ''; % -- write out file for the controls 

% % -- exclude bad sessions 
% excluded_sess = [8,22,30,31];
% excluded_idx = [2,5,8,9];
sess_list = 1:size(sess_info{1},1);
% sess_list(excluded_idx) = [];

name_structure = '/modulators_decod_accuracy.mat';
name_struct_input_SA = '/session_controls_same_area_info.mat'; % -- structure controls Same Area
name_struct_input_OA = '/session_controls_other_areas_info.mat'; % -- structure controls Other Areas


modulators = [];

cnt_el = 1;
cnt = 0;
for i= sess_list  % For each session with at least one modulator
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Sess,name_structure))
    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    
    rec_idx = sess_data.receiver_idx; 
    mod_list = mod_accuracy.mod_idx;
    accuracy = mod_accuracy.auc';
%     accuracy = mod_accuracy.Decod_Accuracy';
    
    idx = find(mod_list == rec_idx); % get the index for the modulator-receiver 

    if ~isempty(idx) % remove element if it is receiver
        cnt = cnt + 1;
        mod_list(idx) = [];
        accuracy(idx) = [];
    end 
    
    sessions = repmat(Sess,1,size(mod_list,2));
    mod_sess = [double(sessions)', double(mod_list)', accuracy] % session, modulator idx, decod accuracy
    modulators = [modulators; mod_sess];
     
end


% modulators(37,:) = [];
% modulators(30,:) = [];
% modulators(9,:) = [];

% modulators(37,:) = [];
% modulators(30,:) = [];
% modulators(9,:) = [];


modulators = [modulators, double(1:size(modulators,1))'] % session, modulator idx, decod accuracy, order in coherence avg structure
% modulators = sortrows(modulators,3,'descend');

dlmwrite(strcat(dir_out,sprintf('/modulators_unsorted%s.txt',filename)),modulators,'delimiter','\t'); % session, modulator idx, decod accuracy, order index i
 
 
%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS SAME AREA
%%%%%%%%%%%%%%%%%%%%%%%%

ctrl_list = [];

for i = sess_list  % For each session with at least one modulator

    Sess = sess_info{1}(i); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Controls_same_area',Sess));
    load(strcat(dir_Sess,name_struct_input_SA)); % load data structure info 
    sess_controls = sess_All_controls_same_area;
%     sess_controls = sess_control_lfp;
    clear sess_All_controls_same_area;
    
    sessions = repmat(Sess,1,size(sess_controls.ctrl_idx,2));
    controls = [double(sessions)', double(sess_controls.ctrl_idx)']; % session, ctrl idx
    ctrl_list = [ctrl_list; controls];
    
end

ctrl_list = [ctrl_list, (1:size(ctrl_list,1))']; % Session, control idx, order index i
dlmwrite(strcat(dir_out,sprintf('/control_list_same_area%s.txt',filename_ctrl)),ctrl_list,'delimiter','\t'); % session, control idx, order index i



%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS OTHER AREAS
%%%%%%%%%%%%%%%%%%%%%%%%

ctrl_list = [];

for i = sess_list  % For each session with at least one modulator

    Sess = sess_info{1}(i); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Controls_other_areas',Sess));
    load(strcat(dir_Sess,name_struct_input_OA)); % RS LFP split into 1 sec window and artifacts removed
    sess_controls = sess_All_controls_other_areas;
%     sess_controls = sess_control_lfp;
    clear sess_All_controls_other_areas;
    
    sessions = repmat(Sess,1,size(sess_controls.ctrl_idx,2));
    controls = [double(sessions)', double(sess_controls.ctrl_idx)']; % session, modulator idx, decod accuracy
    ctrl_list = [ctrl_list; controls];
    
end

ctrl_list = [ctrl_list, (1:size(ctrl_list,1))']; % Session, control idx, order index i
dlmwrite(strcat(dir_out,sprintf('/control_list_other_areas%s.txt',filename_ctrl)),ctrl_list,'delimiter','\t'); % session, control idx, order index i






