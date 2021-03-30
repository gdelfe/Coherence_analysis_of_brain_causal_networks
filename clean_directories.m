%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAN DIRECTORIES FROM FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

list_sess = 1:19;
list_sess(17) = [];

for i = list_sess %list_sess %size(sess_info{1},1) % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session numberdir_Sess
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));
    dir_Controls = strcat(dir_RS,sprintf('/Sess_%d/Controls_same_area',Sess));

    
    cd(dir_Controls)
    !mv sess_controls_same_area_lfp_001.mat session_controls_same_area_lfp_rec001.mat
    %     !mv sess_data_lfp.mat Modulators/session_data_lfp.mat

    
end






