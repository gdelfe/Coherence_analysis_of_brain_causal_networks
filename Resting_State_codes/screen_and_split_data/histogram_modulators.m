
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code loads the entire RS time series for the MODULATORS and check for artifacts in the lfp
% It splits the RS time series into 1 sec windows and remove all the windows
% which contains artifacts (th = 4*std(lfp))
%
% INPUT : sess_data_info.mat
%         structure containing all the info about the session, i.e. idx modulators, sender, etc...
%
% OUTPUT: 1. sess_data_lfp.mat:
%               Strucure of data for each session containing:
%               a. lfp before removing artifacts
%               b. lfp after artifacts were removed
%
%    @ Gino Del Ferraro, December 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band_name = 'beta band';
freq_band = 'beta_band';
monkey = 'Maverick';
dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)
mod_numb = zeros(1,size(sess_info{1},1));

for i= 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all

    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    dir_Sess_mod = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Sess_mod,'/session_data_info.mat')); % --- dataG: all data info and LFP
    
    % -- load list electrodes, sender, receiver
    electrode = sess_data.RecordPair; % ---- all electrode pairs
    receiver = sess_data.receiver_pair;  % ---- receiver pair
    sender = sess_data.sender_pair; % ---- sender pair
    sess_data

    mod_Ch = sess_data.mod_idx; % control modulators
    mod_numb(i) = size(mod_Ch,2);

    
end



[max_mod,idx] = max(mod_numb);

% -- FIGURE: histogram of number of modulators 
fig = figure;
histogram(mod_numb,max_mod,'FaceAlpha',.6); grid on
legend('Numb of modulators')
title(sprintf('%s - %s, Numb of edges',monkey,freq_band_name),'FontSize',12)
ylabel('# of edges')
xlabel('# of modulators')

fig_name = strcat(dir_RS,sprintf('/%s_%s_histogram_modulators.png',monkey,freq_band));
saveas(fig,fig_name);














