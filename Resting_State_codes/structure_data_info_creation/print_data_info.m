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

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

list_sess = 1:length(sess_info{3});

% -- define list of sessions
if strcmp(freq_band,'beta_band')
    if strcmp(monkey,'Maverick')
        list_sess = 1:19;
        list_sess(17) = [];

    end
end

M1 = 0; CN = 0; dPFC = 0; OFC = 0; ACC = 0;
cnt = 0;

for i= list_sess  
    
    Sess = sess_info{1}(i); % Session number
%     display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    dir_ctrl_SA = strcat(dir_RS,sprintf('/Sess_%d/Controls_same_area',Sess));
    dir_ctrl_OA = strcat(dir_RS,sprintf('/Sess_%d/Controls_other_areas',Sess));

    
    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    sess_data
        
    load(strcat(dir_ctrl_SA,'/session_controls_same_area_info.mat')); % --- dataG: all data info and LFP
    ctrl_SA = sess_All_controls_same_area
% 
    load(strcat(dir_ctrl_OA,'/session_controls_other_areas_info.mat')); % --- dataG: all data info and LFP
    ctrl_OA = sess_All_controls_other_areas
    
    
    
%     idx_M1 = strfind(sess_data.mod_areas(1:end), 'M1');
%     M1 = M1 + length(find(not(cellfun('isempty', idx_M1))));
%     
%     idx_CN = strfind(sess_data.mod_areas(1:end), 'CN');
%     CN = CN + length(find(not(cellfun('isempty', idx_CN))));
%     
%     idx_dPFC = strfind(sess_data.mod_areas(1:end), 'dPFC');
%     dPFC = dPFC + length(find(not(cellfun('isempty', idx_dPFC))));
%     
%     idx_OFC = strfind(sess_data.mod_areas(1:end), 'OFC');
%     OFC = OFC + length(find(not(cellfun('isempty', idx_OFC))));
%     
%     idx_ACC = strfind(sess_data.mod_areas(1:end), 'ACC');
%     ACC = ACC + length(find(not(cellfun('isempty', idx_ACC))));
    
%     idx_ACC = strfind(sess_data.mod_areas(1:end), 'ACC');
%     ACC = ACC + length(find(not(cellfun('isempty', idx_ACC))));



end 



