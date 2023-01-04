%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes
%
% Gino Del Ferraro, Jan 2023, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - PATH --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST STATISTICS - BOTH MONKEYS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_band = 'theta_band';
recording_both = 'last_rec-rec001_002';
decoding = 'AUC';
N = 100;
fk = 200;

dir_both_monkeys = strcat(dir_main,sprintf('both_monkeys/%s/modulators_vs_controls/%s/%s',freq_band,recording_both,decoding));

load(strcat(dir_both_monkeys,sprintf('/modulators_N_%d.mat',N))); % structure mod
load(strcat(dir_both_monkeys,sprintf('/controls_same_area_N_%d.mat',N)));
load(strcat(dir_both_monkeys,sprintf('/controls_other_areas_N_%d.mat',N)));

% %%%%%%%%%%%%%%%%%%%%%%%%%
% NULL DISTRIBUTIONS DATA 
% %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% MAVERICK - 

monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_in_Mav = strcat(dir_RS_Theta,'/Modulators_Controls_avg_results');

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

c_ms_perm = [];
c_mr_perm = [];

for i = 1:size(sess_info{1},1)
    Sess = sess_info{1}(i); % Session number
    dir_Mod_permute = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/permutations',Sess));
    load(strcat(dir_Mod_permute,'/coherence_MS_MR_permuted_single.mat'));
    
    for m = 1:size(coh_perm.mod,2)
        c_ms_perm = [c_ms_perm; abs(coh_perm.mod(m).c_ms_perm) ];
        c_mr_perm = [c_mr_perm; abs(coh_perm.mod(m).c_mr_perm) ];
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%
% ARCHIE - 

monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_in_Arc = strcat(dir_RS_Theta,'/Modulators_Controls_avg_results');

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


for i = 1:size(sess_info{1},1)
    Sess = sess_info{1}(i); % Session number
    dir_Mod_permute = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/permutations',Sess));
    load(strcat(dir_Mod_permute,'/coherence_MS_MR_permuted_single.mat'));
    
    for m = 1:size(coh_perm.mod,2)
        c_ms_perm = [c_ms_perm; abs(coh_perm.mod(m).c_ms_perm) ];
        c_mr_perm = [c_mr_perm; abs(coh_perm.mod(m).c_mr_perm) ];
    end
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PVALUES for difference wrt zero

% P-values calculation 
pval_ms = nnz(mean(modulators.mean_coh_ms(9:21)) < c_ms_perm(:,9:21))/size(c_ms_perm,1); % MS-coherence in theta band
pval_mr = nnz(mean(modulators.mean_coh_mr(9:21)) < c_mr_perm(:,9:21))/size(c_mr_perm,1); % MR-coherence in theta band 

% define pval = 0 such as 1/n_iterations
if pval_ms == 0 ; pval_ms = 1/size(c_ms_perm,1); end
if pval_mr == 0 ; pval_mr = 1/size(c_ms_perm,1); end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PVALUES for the difference

% Test statistics -- MS
diff_mod_SA_ms = modulators.mean_coh_ms - ctrl_SA.mean_coh_ms;
diff_mod_OA_ms = modulators.mean_coh_ms - ctrl_OA.mean_coh_ms;

% Null distributions -- MS
load(strcat(dir_in_Mav,'/pseudo_coherence_diff_MS_mod_ctrl_SA.mat'));
mav_mod_SA_ms = diff_SA_ms_all;
load(strcat(dir_in_Mav,'/pseudo_coherence_diff_MS_mod_ctrl_OA.mat'));
mav_mod_OA_ms = diff_OA_ms_all;
load(strcat(dir_in_Arc,'/pseudo_coherence_diff_MS_mod_ctrl_SA.mat'));
arc_mod_SA_ms = diff_SA_ms_all;
load(strcat(dir_in_Arc,'/pseudo_coherence_diff_MS_mod_ctrl_OA.mat'));
arc_mod_OA_ms = diff_OA_ms_all;

diff_mod_SA_ms_ALL = [mav_mod_SA_ms; arc_mod_SA_ms];
diff_mod_OA_ms_ALL = [mav_mod_OA_ms; arc_mod_OA_ms];

pval_diff_mod_SA_ms = nnz(mean(diff_mod_SA_ms(9:21)) < diff_mod_SA_ms_ALL(:,9:21) | -abs(mean(diff_mod_SA_ms(9:21))) > diff_mod_SA_ms_ALL(:,9:21))/size(diff_mod_SA_ms_ALL,1); % MS-coherence in theta band
pval_diff_mod_OA_ms = nnz(mean(diff_mod_OA_ms(9:21)) < diff_mod_OA_ms_ALL(:,9:21) | -abs(mean(diff_mod_OA_ms(9:21))) > diff_mod_OA_ms_ALL(:,9:21))/size(diff_mod_OA_ms_ALL,1); % MS-coherence in theta band




% Test statistics -- MR
diff_mod_SA_mr = modulators.mean_coh_mr - ctrl_SA.mean_coh_mr;
diff_mod_OA_mr = modulators.mean_coh_mr - ctrl_OA.mean_coh_mr;

% Null distributions -- MS
load(strcat(dir_in_Mav,'/pseudo_coherence_diff_MR_mod_ctrl_SA.mat'));
mav_mod_SA_mr = diff_SA_mr_all;
load(strcat(dir_in_Mav,'/pseudo_coherence_diff_MR_mod_ctrl_OA.mat'));
mav_mod_OA_mr = diff_OA_mr_all;
load(strcat(dir_in_Arc,'/pseudo_coherence_diff_MR_mod_ctrl_SA.mat'));
arc_mod_SA_mr = diff_SA_mr_all;
load(strcat(dir_in_Arc,'/pseudo_coherence_diff_MR_mod_ctrl_OA.mat'));
arc_mod_OA_mr = diff_OA_mr_all;

diff_mod_SA_mr_ALL = [mav_mod_SA_mr; arc_mod_SA_mr];
diff_mod_OA_mr_ALL = [mav_mod_OA_mr; arc_mod_OA_mr];

pval_diff_mod_SA_mr = nnz(abs(mean(diff_mod_SA_mr(9:21))) < diff_mod_SA_mr_ALL(:,9:21) | -abs(mean(diff_mod_SA_mr(9:21))) > diff_mod_SA_mr_ALL(:,9:21)  )/size(diff_mod_SA_mr_ALL,1); % MS-coherence in theta band
pval_diff_mod_OA_mr = nnz(mean(diff_mod_OA_mr(9:21)) < diff_mod_OA_mr_ALL(:,9:21) | -abs(mean(diff_mod_OA_mr(9:21))) > diff_mod_OA_mr_ALL(:,9:21))/size(diff_mod_OA_mr_ALL,1); % MS-coherence in theta band

figure;
histogram(diff_SA_mr_all(:,9),50)

