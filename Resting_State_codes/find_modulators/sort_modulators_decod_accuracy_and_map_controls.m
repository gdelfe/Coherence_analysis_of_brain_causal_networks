%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code sorts the modulators according to their decoding accuracy
% For all the controls, same area and different area, it lists them
% together with the session label and the i indx
% 
% @ Gino Del Ferraro, NYU, March 2021



clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
monkey = 'Archie';
dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));
dir_out = strcat(dir_main,sprintf('%s/Resting_state/%s/Modulators_controls',monkey,freq_band));

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)


name_structure = '/modulators_decod_accuracy.mat';
% 
% % -- define list of sessions
% if strcmp(monkey,'Maverick')
%     list_sess = 1:19;
%     list_sess(17) = [];
% else
%     list_sess = 1:length(sess_info{3});
% end

list_sess = 1:length(sess_info{3});
modulators = [];
cnt_el = 1;
for i= list_sess  % For each session with at least one modulator
    
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
        mod_list(idx) = [];
        accuracy(idx) = [];
    end 
    
    sessions = repmat(Sess,1,size(mod_list,2));
    mod_sess = [double(sessions)', double(mod_list)', accuracy] % session, modulator idx, decod accuracy
    modulators = [modulators; mod_sess];
     
end

modulators = [modulators, double(1:size(modulators,1))']
modulators = sortrows(modulators,3,'descend');
dlmwrite(strcat(dir_out,'/modulators_sorted_decod_accuracy_AUC.txt'),modulators,'delimiter','\t'); % session, modulator idx, decod accuracy, order index i

keyboard 
%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS SAME AREA
%%%%%%%%%%%%%%%%%%%%%%%%

name_struct_input = '/session_controls_same_area_info.mat';
ctrl_list = [];

for i = list_sess  % For each session with at least one modulator

    Sess = sess_info{1}(i); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Controls_same_area',Sess));
    load(strcat(dir_Sess,name_struct_input)); % load data structure info 
    sess_controls = sess_All_controls_same_area;
    clear sess_All_controls_same_area;
    
    sessions = repmat(Sess,1,size(sess_controls.ctrl_idx,2));
    controls = [double(sessions)', double(sess_controls.ctrl_idx)']; % session, ctrl idx
    ctrl_list = [ctrl_list; controls];
    
end

ctrl_list = [ctrl_list, (1:size(ctrl_list,1))']; % Session, control idx, order index i
dlmwrite(strcat(dir_out,'/control_list_same_area.txt'),ctrl_list,'delimiter','\t'); % session, control idx, order index i



%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS OTHER AREAS
%%%%%%%%%%%%%%%%%%%%%%%%

name_struct_input = '/session_controls_other_areas_info.mat';
ctrl_list = [];

for i = list_sess  % For each session with at least one modulator

    Sess = sess_info{1}(i); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Controls_other_areas',Sess));
    load(strcat(dir_Sess,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    sess_controls = sess_All_controls_other_areas;
    clear sess_All_controls_other_areas;
    
    sessions = repmat(Sess,1,size(sess_controls.ctrl_idx,2));
    controls = [double(sessions)', double(sess_controls.ctrl_idx)']; % session, modulator idx, decod accuracy
    ctrl_list = [ctrl_list; controls];
    
end

ctrl_list = [ctrl_list, (1:size(ctrl_list,1))']; % Session, control idx, order index i
dlmwrite(strcat(dir_out,'/control_list_other_areas.txt'),ctrl_list,'delimiter','\t'); % session, control idx, order index i

