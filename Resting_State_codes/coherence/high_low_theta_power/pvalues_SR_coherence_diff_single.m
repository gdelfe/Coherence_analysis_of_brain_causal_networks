clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - ADD PATH  --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
freq_band = 'theta_band';

dir_permutations_mav = strcat(dir_main,'/Maverick/Resting_State/high_low_theta/permutation_test/permutation_single_distribution');
dir_high_low_theta_mav = strcat(dir_main,'/Maverick/Resting_State/high_low_theta');
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


% TEST STATISTICS (modulators, SA, OA) ------------------------
% high power 
mean_sr_H = mean(abs(c_sr_H));
mean_sr_SA_H = mean(abs(c_sr_SA_H));
mean_sr_OA_H = mean(abs(c_sr_OA_H));

% low power
mean_sr_L = mean(abs(c_sr_L));
mean_sr_SA_L = mean(abs(c_sr_SA_L));
mean_sr_OA_L = mean(abs(c_sr_OA_L));

% DIFFERENCES in the theta range 4-10 Hz

test_diff_mod = mean(mean_sr_H(9:20)) - mean(mean_sr_L(9:20))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fk = 200;
f = linspace(1,fk,409); % frequency values (range)
















% Coherence differences permutation test distributions 
load(strcat(dir_permutations_mav,'/coh_diff_permutation_mav.mat'))
diff_mod_mav = diff;
load(strcat(dir_permutations_mav,'/coh_diff_perm_ctlr_SA_mav.mat'))
diff_SA_mav = diff;
load(strcat(dir_permutations_mav,'/coh_diff_perm_ctlr_OA_mav.mat'))
diff_OA_mav = diff;

% Test-statistics SR coherence difference - modulators
load(strcat(dir_high_low_theta_mav,'/coh_all_sess_sr_high.mat'))
load(strcat(dir_high_low_theta_mav,'/coh_all_sess_sr_low.mat'))

% Test-statistics SR coherence difference - controls SA
load(strcat(dir_high_low_theta_mav,'/coh_all_sess_controls_sender_SA-receiver.mat'))
% control-sender receiver coherence
coh_cr_SA_high = ctrl_send_SA_coh.cr_high;
coh_cr_SA_low = ctrl_send_SA_coh.cr_low;

% Test-statistics SR coherence difference - controls OA
load(strcat(dir_high_low_theta_mav,'/coh_all_sess_controls_sender_OA-receiver.mat'))
% control-sender receiver coherence
coh_cr_OA_high = ctrl_send_OA_coh.cr_high;
coh_cr_OA_low = ctrl_send_OA_coh.cr_low;

f = linspace(0,200,409);

% sender and receiver 
mean_sr_high = mean(abs(coh_all_c_sr_high),1);
mean_sr_low = mean(abs(coh_all_c_sr_low),1);

% controls SA sender and receveir 
mean_cr_SA_high = mean(abs(coh_cr_SA_high),1);
mean_cr_SA_low = mean(abs(coh_cr_SA_low),1);


% controls OA sender and receveir 
mean_cr_OA_high = mean(abs(coh_cr_OA_high),1);
mean_cr_OA_low = mean(abs(coh_cr_OA_low),1);


test_diff_mod = mean(mean_sr_high(9:22)) - mean(mean_sr_low(9:22)); % test statistics - mean difference modulator in theta range



figure;
histogram(mean(diff_mod_mav(:,9:22),2),100,'FaceAlpha',0.5,'Normalization','pdf')


