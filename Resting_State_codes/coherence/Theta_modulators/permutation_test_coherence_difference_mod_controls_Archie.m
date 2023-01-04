
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the null distribution for the test the Null Hypothesis
% that there is no difference between the MS coherence and the Control-S
% coherence (same for the receiver). In other words, the null that the
% coherence difference MS-CS and MR-CR is not different from zero.
%
% This is used in Fig 5 and 6 in the theta paper
%
% Gino Del Ferraro, Jan 2023, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - PATH --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

name_modulators = '/coh_spec_m_fk_200_W_5_movie.mat';
name_SA = '/coh_spec_m_Controls_same_area_fk_200_W_5_movie.mat'; % -- name file to load 
name_OA = '/coh_spec_m_Controls_other_areas_fk_200_W_5_movie.mat'; % -- name file to load 

filename = '.mat'; % -- filename for sess_data_info.mat 

freq_band = 'theta_band';
monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_out = strcat(dir_RS_Theta,'/Modulators_Controls_avg_results');
iterations = 10000;

% Define directory path for each session
dir_input = strcat(dir_RS_Theta,'/Modulators_Controls_avg_results/movie_all_sessions');

% Load coherence files for modulators and controls
load(strcat(dir_input,name_modulators)); % RS LFP split into 1 sec window and artifacts removed
modulators = mod;
load(strcat(dir_input,name_SA)); % RS LFP split into 1 sec window and artifacts removed
ctrl_SA = mod;
load(strcat(dir_input,name_OA)); % RS LFP split into 1 sec window and artifacts removed
ctrl_OA = mod;
clear mod
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert structures into matrices

% Modulators
mod_mat = cell2mat(struct2cell(modulators)); % transform struct to mat for modulators
coh_ms = sq(mod_mat(1,:,:))'; % 1st field, c_ms
coh_mr = sq(mod_mat(2,:,:))'; %  2nd field, c_mr

% controls SA
ctrl_SA_mat = cell2mat(struct2cell(ctrl_SA)); % transform struct to mat for modulators
ctrl_SA_ms = sq(ctrl_SA_mat(1,:,:))'; % 1st field, c_ms
ctrl_SA_mr = sq(ctrl_SA_mat(2,:,:))'; %  2nd field, c_mr

% controls OA 
ctrl_OA_mat = cell2mat(struct2cell(ctrl_OA)); % transform struct to mat for modulators
ctrl_OA_ms = sq(ctrl_OA_mat(1,:,:))'; % 1st field, c_ms
ctrl_OA_mr = sq(ctrl_OA_mat(2,:,:))'; %  2nd field, c_mr

L_m = size(modulators,2);
L_SA = size(ctrl_SA,2);
L_OA = size(ctrl_OA,2);

% matrices to store the null distributions 
diff_SA_ms_all = zeros(iterations,409);
diff_OA_ms_all = zeros(iterations,409);
diff_SA_mr_all = zeros(iterations,409);
diff_OA_mr_all = zeros(iterations,409);

for iter=1:iterations
    
    % %%%%%%%%%%%%%%%%%%%
    % MS Coherence Null distribution 

    % Pick as many controls as modulators and pick them randomly  
    perm_SA = randperm(L_SA);
    perm_OA = randperm(L_OA);
    perm_SA = perm_SA(1:L_m); % same length (statistics) of the coherence MS
    perm_OA = perm_OA(1:L_m);
    
    % randomly pick X coherences in the ctrl SA and OA, where X=dim of coherences in  modulators
    SA_coh_ms = ctrl_SA_ms(perm_SA,:);
    OA_coh_ms = ctrl_OA_ms(perm_OA,:);
    
    % concatenate modulator and controls coherences
    mod_SA_ms = [coh_ms;SA_coh_ms];
    mod_OA_ms = [coh_ms;OA_coh_ms];
    
    % random permutation for null distribution 
    p_SA = randperm(size(mod_SA_ms,1));
    p_OA = randperm(size(mod_OA_ms,1));
    
    % permute the concatenated coherences in order to create two pseudo sets 
    mod_SA_ms = mod_SA_ms(p_SA,:);
    mod_OA_ms = mod_SA_ms(p_OA,:);
    
    % Pseudo coherence difference modulator-controls_SA
    pseudo_mod_SA_ms_1 = mod_SA_ms(1:round(size(mod_SA_ms,1)/2),:); % pseudo data set 1
    pseudo_mod_SA_ms_2 = mod_SA_ms(round(size(mod_SA_ms,1)/2)+1:end,:); % pseudo data set 2
    diff_SA_ms = mean(abs(pseudo_mod_SA_ms_1) - abs(pseudo_mod_SA_ms_2)); % pseudo coherence difference 
    
    pseudo_mod_OA_ms_1 = mod_OA_ms(1:round(size(mod_OA_ms,1)/2),:); % pseudo data set 1
    pseudo_mod_OA_ms_2 = mod_OA_ms(round(size(mod_OA_ms,1)/2)+1:end,:); % pseudo data set 2
    diff_OA_ms = mean(abs(pseudo_mod_OA_ms_1) - abs(pseudo_mod_OA_ms_2)); % pseudo coherence difference
    
    % MS pseudo differences 
    diff_SA_ms_all(iter,:) = diff_SA_ms;
    diff_OA_ms_all(iter,:) = diff_OA_ms;
    
    
  
    
    
    % %%%%%%%%%%%%%%%%%%%
    % MR Coherence Null distribution 

    % Pick as many controls as modulators and pick them randomly  
    perm_SA = randperm(L_SA);
    perm_OA = randperm(L_OA);
    perm_SA = perm_SA(1:L_m); % same length (statistics) of the coherence MS
    perm_OA = perm_OA(1:L_m);
    
    % randomly pick X coherences in the ctrl SA and OA, where X=dim of coherences in  modulators
    SA_coh_mr = ctrl_SA_mr(perm_SA,:);
    OA_coh_mr = ctrl_OA_mr(perm_OA,:);
    
    % concatenate modulator and controls coherences
    mod_SA_mr = [coh_mr; SA_coh_mr];
    mod_OA_mr = [coh_mr; OA_coh_mr];
    
    % random permutation for null distribution 
    p_SA = randperm(size(mod_SA_mr,1));
    p_OA = randperm(size(mod_OA_mr,1));
    
    % permute the concatenated coherences in order to create two pseudo sets 
    mod_SA_mr = mod_SA_mr(p_SA,:);
    mod_OA_mr = mod_SA_mr(p_OA,:);
    
    % Pseudo coherence difference modulator-controls_SA
    pseudo_mod_SA_mr_1 = mod_SA_mr(1:round(size(mod_SA_mr,1)/2),:); % pseudo data set 1
    pseudo_mod_SA_mr_2 = mod_SA_mr(round(size(mod_SA_mr,1)/2)+1:end,:); % pseudo data set 2
    diff_SA_mr = mean(abs(pseudo_mod_SA_mr_1) - abs(pseudo_mod_SA_mr_2)); % pseudo coherence difference 
    
    pseudo_mod_OA_mr_1 = mod_OA_mr(1:round(size(mod_OA_mr,1)/2),:); % pseudo data set 1
    pseudo_mod_OA_mr_2 = mod_OA_mr(round(size(mod_OA_mr,1)/2)+1:end,:); % pseudo data set 2
    diff_OA_mr = mean(abs(pseudo_mod_OA_mr_1) - abs(pseudo_mod_OA_mr_2)); % pseudo coherence difference
    
    % MS pseudo differences 
    diff_SA_mr_all(iter,:) = diff_SA_mr;
    diff_OA_mr_all(iter,:) = diff_OA_mr;
    
    
end

save(strcat(dir_out,'/pseudo_coherence_diff_MS_mod_ctrl_SA.mat'),'diff_SA_ms_all');
save(strcat(dir_out,'/pseudo_coherence_diff_MS_mod_ctrl_OA.mat'),'diff_OA_ms_all');
save(strcat(dir_out,'/pseudo_coherence_diff_MR_mod_ctrl_SA.mat'),'diff_SA_mr_all');
save(strcat(dir_out,'/pseudo_coherence_diff_MR_mod_ctrl_OA.mat'),'diff_OA_mr_all');







   