
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THETA PICK: 7.8 Hz, BETA PICK: 21 Hz

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)



addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/STIM_codes');
dir_Stim = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';
dir_Perm = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data/Permutation_test';


fk = 200; W = 5; iter = 1000;

% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_Stim,sprintf('/coh_spec_mr_fk_%d_W_%d.mat',fk,W))); % structure stim
modulators = mean_coh_and_spec_STIM(stim);

% %%%%%%%%% PERMUTED DATA  %%%%%%
load(strcat(dir_Perm,sprintf('/coh_spec_mr_permuted_fk_%d_W_%d.mat',fk,W))); % structure stim_perm
coh = stim_perm;
f = linspace(1,fk,size(stim(1).c_mr,2)); % frequency values (range)

% Plot to check --- permuted coherence 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all
el = 41; permS = 1 ; permR = 1;
plot(abs(coh(el).perm(permR).c_mr))
xlabel('freq (Hz)')
ylabel('coherence')
legend('Permuted coherence')
grid on


% Theta pick 7.3 Hz -> f(14)
% Beta pick 21.9 Hz -> f(44)
t_idx = 14; b_idx = 44;

theta_MR = zeros(size(coh,2),iter);
beta_MR = zeros(size(coh,2),iter);

% select only theta frequence and beta frequence at the peaks for each
% electrode 
% Modulators-Receivers 
for ch = 1:size(coh,2)
    for perm = 1:iter
        
        theta_MR(ch,perm) = coh(ch).perm(perm).c_mr(t_idx); % 7.8 Hz
        beta_MR(ch,perm) = coh(ch).perm(perm).c_mr(b_idx);  % 21 Hz
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THETA AND BETA TEST OF INDEPENDENCE --- MODULATORS-RECEIVERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pvalue calculation for each electrodes 
pth = 0.05;
theta_NZ = 0; % counts for theta non-zero
beta_NZ = 0; % counts for beta non-zero
for ch = 1:size(coh,2)
    
   stats.MR.Ch(ch).theta_MR_pval =  nnz(abs(theta_MR(ch,:)) > abs(stim(ch).c_mr(t_idx)))/iter;
   stats.MR.Ch(ch).beta_MR_pval =  nnz(abs(beta_MR(ch,:)) > abs(stim(ch).c_mr(b_idx)))/iter;
   
   % How many times theta and beta are significantly larger than zero (with
   % pval threshold = 0.05)
   theta_NZ = theta_NZ + int8(stats.MR.Ch(ch).theta_MR_pval <= pth); % if pval of the channel is smaller than pval threshold add 1
   beta_NZ = beta_NZ + int8(stats.MR.Ch(ch).beta_MR_pval <= pth); 
   
end


% counts of theta being significantly zero and beta being significantly zero
theta_Z = size(coh,2) - theta_NZ;
beta_Z = size(coh,2) - beta_NZ;

Obs = double([theta_NZ, theta_Z; beta_NZ, beta_Z]); 

stats.MR.Obs = Obs; % save into structure 

N = sum(sum(Obs)) % Sample size
E_theta_NZ = sum(Obs(1,:)) * sum(Obs(:,1))/N
E_theta_Z = sum(Obs(1,:)) * sum(Obs(:,2))/N
E_beta_NZ = sum(Obs(2,:)) * sum(Obs(:,1))/N
E_beta_Z = sum(Obs(2,:)) * sum(Obs(:,2))/N

Exp = [E_theta_NZ,E_theta_Z; E_beta_NZ,E_beta_Z];

stats.MR.Exp = Exp; % save into structure 

ch = 41;
fig = figure;
histogram(abs(theta_MR(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
hold on; histogram(abs(beta_MR(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
legend('p-val theta','p-val beta')
title('p-val distribution','FontSize',12)


O = reshape(Obs,[1,4]);
E = reshape(Exp,[1,4]);0.0
% Compute Chi-Squared 
chi_MR = 0;
for i = 1:4
    power((O(i) - E(i)),2)/E(i)
    chi_MR = chi_MR + power((O(i) - E(i)),2)/E(i);
    
end
chi_MR

stats.MR.ChiSquared = chi_MR; % save into structure 




