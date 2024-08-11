


clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/vol/bd5/People/Gino/Gino_codes')
dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)
% sub_dir = ''
filename = '.mat'; % -- filename for sess_data_info.mat 


for i= 1:size(sess_info{1},1) 



    Sess = sess_info{1}(i); % Session number
    dir_Sess_mod = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Sess_mod,'/session_data_info.mat')); % --- dataG: all data info and LFP
    keyboard
end