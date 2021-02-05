
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs the permutation test for the coherence: it shuffles 
% the trail of one of the two channels when computing the MR, SR, MS
% coherence
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
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
step = 110;

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

name_struct_input = '/sess_data_lfp.mat';

iter = 1000; % number of iteration for the permutation test
cnt_sr = 1; % counter sender-receiver coherencies
cnt_el = 1; % counter for how many modulators excluding the receivers modulators
list_sess = 1:19;
list_sess(17) = []; % -- Session 17 and 20 are full of artifacts

max = 0;
for i = list_sess %1:size(sess_info{1},1)-1  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));
    
    load(strcat(dir_Sess,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
     mod_Ch = sess_data_lfp.mod_idx; % -- modulators (not controls!) index
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
    
    
    % %%%%%%% ALL Electrodes LFP %%%%%%%%%%%%%%%%%%%%%
    lfp_E_all = sess_data_lfp.lfp_E;
    
    cnt_m = 1;
    for Ch = mod_Ch % for all the modulators in the session
        
        if max < size(sess_data_lfp.outliers_tot(cnt_m).idx,2)
            max = size(sess_data_lfp.outliers_tot(cnt_m).idx,2)    
        end
        cnt_m = cnt_m + 1;
    end
    
    
    
end

max
% % MAX number of outliers is max = 52

    
    
    
    