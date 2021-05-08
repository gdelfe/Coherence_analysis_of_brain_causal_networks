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
monkey = 'Archie';
dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));
% dir_Stim = strcat(dir_main,sprintf('%s/Stim_data',monkey));


fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% -- define list of sessions
% if strcmp(monkey,'Maverick')
%     list_sess = 1:19;
%     list_sess(17) = [];
% else
%     list_sess = 1:length(sess_info{3});
% end

% sess_list = importdata(strcat(dir_Stim,'/Sessions_list.txt'));
% cd(dir_Stim)

for i = 1:size(sess_info{1},1) % For each session with at least one modulator
    
    
%     close all
    Sess = sess_info{1}(i); % Session numberdir_Sess
%     display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    cd(dir_Sess)
%     dir_Controls = strcat(dir_RS,sprintf('/Sess_%d/Controls_other_areas',Sess));
%     dir_Sess = strcat(dir_Stim,sprintf('/Sess_%d/Theta_band',i));
%     if isempty(find(sess_info{1}== i))
%         system(sprintf('rm -rf Sess_%d',i))
%     end 
     
    !mv sess_data_lfp_coherence_fk_200_W_5_movie.mat movie
%     !rm -rf last_recording
%     !rm -rf movie
%     !mv sess_all_controls_other_areas_lfp_001.mat session_controls_other_areas_lfp_rec001.mat
%     !mv sess_all_controls_other_areas_lfp.mat session_controls_other_areas_lfp_movie.mat

%     system(sprintf('mkdir %s/Theta_band/Sess_%d',dir_Stim,i))
%     system(sprintf('mv Data_with_theta_band.mat %s/Theta_band/Sess_%d',dir_Stim,i))
%     cd('..')
%     !rm -rf Theta_band
    %     !mv sess_data_lfp.mat Modulators/session_data_lfp.mat

    
end






