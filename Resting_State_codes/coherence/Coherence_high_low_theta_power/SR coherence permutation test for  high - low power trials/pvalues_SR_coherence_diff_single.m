clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - ADD PATH  --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
freq_band = 'theta_band';

dir_permutations_mav = strcat(dir_main,'/Maverick/Resting_State/high_low_theta/permutation_test/permutation_single_distribution');
dir_permutations_arc = strcat(dir_main,'/Archie/Resting_State/high_low_theta/permutation_test/permutation_single_distribution');
dir_HL_theta_Mav = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';
dir_HL_theta_Arc = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Archie/Resting_State/high_low_theta';
dir_both_monkeys = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys/high_low_theta';

% %%%%%%%%% Maverick High power %%%%%%
c_sr_mod_mav_H = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_high.mat')); % SR coherence - high theta - trial sorted by modulator power  
c_sr_mod_SA_mav_H = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_high_controls_SA.mat')); % SR coherence - high theta -  trials sorted by controls SA
c_sr_mod_OA_mav_H = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_high_controls_OA.mat')); % SR coherence - high theta - trials sorted by controls OA 

% %%%%%%%%% Maverick Low power %%%%%%
c_sr_mod_mav_L = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_low.mat')); % SR coherence - high theta - trial sorted by modulator power  
c_sr_mod_SA_mav_L = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_low_controls_SA.mat')); % SR coherence - high theta -  trials sorted by controls SA
c_sr_mod_OA_mav_L = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_low_controls_OA.mat')); % SR coherence - high theta - trials sorted by controls OA 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%  Archie High power %%%%%%
c_sr_mod_arc_H = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_high.mat')); % SR coherence - high theta - trial sorted by modulator power  
c_sr_mod_SA_arc_H = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_high_controls_SA.mat')); % SR coherence - high theta - trials sorted by controls SA
c_sr_mod_OA_arc_H = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_high_controls_OA.mat')); % SR coherence - high theta - trials sorted by controls OA

% %%%%%%%%% Archie Low power %%%%%%
c_sr_mod_arc_L = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_low.mat')); % SR coherence - high theta - trial sorted by modulator power  
c_sr_mod_SA_arc_L = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_low_controls_SA.mat')); % SR coherence - high theta -  trials sorted by controls SA
c_sr_mod_OA_arc_L = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_low_controls_OA.mat')); % SR coherence - high theta - trials sorted by controls OA 



% %%%% Concatenate results for Maverick and Archie 
% high power
c_sr_H = [c_sr_mod_mav_H.coh_all_c_sr_high; c_sr_mod_arc_H.coh_all_c_sr_high]; % concatenate data for the two monkeys - trial sorted by modulator power 
c_sr_SA_H = [c_sr_mod_SA_mav_H.coh_all_c_sr_high; c_sr_mod_SA_arc_H.coh_all_c_sr_high]; % concatenate data for the two monkeys - trial sorted by modulator power 
c_sr_OA_H = [c_sr_mod_OA_mav_H.coh_all_c_sr_high; c_sr_mod_OA_arc_H.coh_all_c_sr_high]; % concatenate data for the two monkeys - trial sorted by modulator power 
% low power
c_sr_L = [c_sr_mod_mav_L.coh_all_c_sr_low; c_sr_mod_arc_L.coh_all_c_sr_low]; % concatenate data for the two monkeys - trial sorted by modulator power 
c_sr_SA_L = [c_sr_mod_SA_mav_L.coh_all_c_sr_low; c_sr_mod_SA_arc_L.coh_all_c_sr_low]; % concatenate data for the two monkeys - trial sorted by modulator power 
c_sr_OA_L = [c_sr_mod_OA_mav_L.coh_all_c_sr_low; c_sr_mod_OA_arc_L.coh_all_c_sr_low]; % concatenate data for the two monkeys - trial sorted by modulator power 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST STATISTICS (modulators, SA, OA) ------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% high power 
mean_sr_H = mean(abs(c_sr_H));    % modulators 
mean_sr_SA_H = mean(abs(c_sr_SA_H));  % controls SA
mean_sr_OA_H = mean(abs(c_sr_OA_H)); % controls OA 

% low power
mean_sr_L = mean(abs(c_sr_L));
mean_sr_SA_L = mean(abs(c_sr_SA_L));
mean_sr_OA_L = mean(abs(c_sr_OA_L));

% DIFFERENCES in the theta range 4-10 Hz

test_diff_mod = mean(mean_sr_H(9:20)) - mean(mean_sr_L(9:20))
test_diff_SA = mean(mean_sr_SA_H(9:20)) - mean(mean_sr_SA_L(9:20))
test_diff_OA = mean(mean_sr_OA_H(9:20)) - mean(mean_sr_OA_L(9:20))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fk = 200;
f = linspace(1,fk,409); % frequency values (range)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD NULL DISTRIBUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Coherence differences permutation test distributions Maverick
load(strcat(dir_permutations_mav,'/coh_diff_permutation_mav.mat'))
diff_mod_mav = diff;
load(strcat(dir_permutations_mav,'/coh_diff_perm_ctlr_SA_mav.mat'))
diff_SA_mav = diff;
load(strcat(dir_permutations_mav,'/coh_diff_perm_ctlr_OA_mav.mat'))
diff_OA_mav = diff;

% Coherence differences permutation test distributions Archie
load(strcat(dir_permutations_arc,'/coh_diff_permutation_archie.mat'))
diff_mod_arc = diff;
load(strcat(dir_permutations_arc,'/coh_diff_perm_ctlr_SA_archie.mat'))
diff_SA_arc = diff;
load(strcat(dir_permutations_arc,'/coh_diff_perm_ctlr_OA_archie.mat'))
diff_OA_arc = diff;


f = linspace(0,200,409);

% Concatenate null distributions for Maverick and Archie
diff_distr_mod = [mean(diff_mod_mav(:,9:22),2); mean(diff_mod_arc(:,9:22),2)];
diff_distr_SA = [mean(diff_SA_mav(:,9:22),2); mean(diff_SA_arc(:,9:22),2)];
diff_distr_OA = [mean(diff_OA_mav(:,9:22),2); mean(diff_OA_arc(:,9:22),2)];

% compute p-values for modulators, controls SA and OA
pvalue_mod = nnz(diff_distr_mod > test_diff_mod)/length(diff_distr_mod)
pvalue_SA = nnz(diff_distr_SA > test_diff_SA)/length(diff_distr_SA)
pvalue_OA = nnz(diff_distr_OA > test_diff_OA)/length(diff_distr_OA)



figure;
histogram(diff_distr_mod,100,'FaceAlpha',0.5,'Normalization','pdf')







