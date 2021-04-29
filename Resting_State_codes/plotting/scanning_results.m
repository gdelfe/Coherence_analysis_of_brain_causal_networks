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
monkey = 'Archie';
recording = 'rec001';


dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));
dir_out = strcat(dir_main,sprintf('%s/Resting_state/%s/Modulators_controls',monkey,freq_band));

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

name_structure = '/session_data_info.mat';

list_sess = 1:length(sess_info{3});
modulators = [];
cnt_el = 1;
for i= list_sess  

    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    load(strcat(dir_Sess,name_structure))
    
end





