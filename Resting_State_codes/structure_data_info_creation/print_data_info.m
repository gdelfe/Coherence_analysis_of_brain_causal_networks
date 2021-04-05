clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'beta_band';
monkey = 'Maverick';
dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

% -- define list of sessions
if strcmp(monkey,'Maverick')
    list_sess = 1:19;
    list_sess(17) = [];
else
    list_sess = 1:length(sess_info{3});
end

for i= list_sess  
    
    Sess = sess_info{1}(i); % Session number
%     display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    sess_data
    
end 



