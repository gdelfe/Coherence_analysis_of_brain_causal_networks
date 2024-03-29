
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theta/Beta chi-squared independence test      
% with average values in an interval rather than a single value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Perm = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state/Permutation_test2';


fk = 200; W = 5; iter = 1000;

fk = 200; W = 5;
% %%%%%%%%% MODULATORS  %%%%%%
load(strcat(dir_RS,sprintf('/coh_spec_m_fk_%d_W_%d.mat',fk,W))); % structure mod
load(strcat(dir_RS,sprintf('/coh_spec_sr_fk_%d_W_%d.mat',fk,W))); % structure stim
stim_mod = stim;
mod_mod = mod;
modulators = mean_coh_and_spec_RS(mod,stim);

% %%%%%%%%% PERMUTED DATA  %%%%%%
load(strcat(dir_Perm,sprintf('/coh_spec_m_fk_%d_W_%d.mat',fk,W))); % structure mod
load(strcat(dir_Perm,sprintf('/coh_sr_permuted_fk_%d_W_%d.mat',fk,W))); % structure stim
f = linspace(1,fk,size(coh(1).perm(1).c_mr,2)); % frequency values (range)

% Plot to check --- permuted coherence 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all
el = 41; permS = 1 ; permR = 1;
plot(abs(coh(el).perm(permS).c_ms))
hold on 
plot(abs(coh(el).perm(permR).c_mr))
xlabel('freq (Hz)')
ylabel('coherence')
legend('Permuted coherence')


% Theta pick 7.8 Hz -> f(15)
% Beta pick 21 Hz -> f(42)
% M-S first peak 8.8 Hz, 2nd peak is the same as MR, i.e. 21 Hz
t_idx = 15; b_idx = 42; s_idx = 17;

theta_MR = zeros(size(coh,2),iter);
beta_MR = zeros(size(coh,2),iter);

theta_MS = zeros(size(coh_sr,2),iter);
beta_MS = zeros(size(coh_sr,2),iter);

% --- select theta and beta within a range around the peaks
% Modulators-Receivers 
for ch = 1:size(coh,2)
    for perm = 1:iter
        
        theta_MR(ch,perm) = sum(abs(coh(ch).perm(perm).c_mr(13:17)))/5; % 6.8 - 8.8 Hz
        beta_MR(ch,perm) = sum(abs(coh(ch).perm(perm).c_mr(40:44)))/5;  % 20-22 Hz
        
    end
end

% Modulators-Senders 
for ch = 1:size(coh_sr,2)
    for perm = 1:iter
        
        theta_MS(ch,perm) = sum(abs(coh_sr(ch).perm(perm).c_sr(15:19)))/5; % 7.8 - 9.8 Hz
        beta_MS(ch,perm) = sum(abs(coh_sr(ch).perm(perm).c_sr(40:44)))/5;  % 20 - 22 Hz
        
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
    
   stats_avg.MR.Ch(ch).theta_MR_pval =  nnz(theta_MR(ch,:) > sum(abs(mod(ch).c_mr(13:17)))/5)/iter;
   stats_avg.MR.Ch(ch).beta_MR_pval =  nnz(beta_MR(ch,:) > sum(abs(mod(ch).c_mr(40:44)))/5)/iter;
   
   % How many times theta and beta are significantly larger than zero (with
   % pval threshold = 0.05)
   theta_NZ = theta_NZ + int8(stats_avg.MR.Ch(ch).theta_MR_pval <= pth); % if pval of the channel is smaller than pval threshold add 1
   beta_NZ = beta_NZ + int8(stats_avg.MR.Ch(ch).beta_MR_pval <= pth); 
   
end

% counts of theta being significantly zero and beta being significantly zero
theta_Z = size(coh,2) - theta_NZ;
beta_Z = size(coh,2) - beta_NZ;

Obs = double([theta_NZ, theta_Z; beta_NZ, beta_Z])

stats_avg.MR.Obs = Obs; % save into structure 

N = sum(sum(Obs)) % Sample size
E_theta_NZ = sum(Obs(1,:)) * sum(Obs(:,1))/N;
E_theta_Z = sum(Obs(1,:)) * sum(Obs(:,2))/N;
E_beta_NZ = sum(Obs(2,:)) * sum(Obs(:,1))/N;
E_beta_Z = sum(Obs(2,:)) * sum(Obs(:,2))/N;

Exp = [E_theta_NZ,E_theta_Z; E_beta_NZ,E_beta_Z]

stats_avg.MR.Exp = Exp; % save into structure 

ch = 1;
fig = figure;
histogram(abs(theta_MR(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
hold on; histogram(abs(beta_MR(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
legend('p-val theta','p-val beta')
title('p-val distribution','FontSize',12)


O = reshape(Obs,[1,4]);
E = reshape(Exp,[1,4]);

% Compute Chi-Squared 
chi_MR = 0;
for i = 1:4
    power((O(i) - E(i)),2)/E(i);
    chi_MR = chi_MR + power((O(i) - E(i)),2)/E(i);
    
end
chi_MR

stats_avg.MR.ChiSquared = chi_MR; % save into structure 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THETA AND BETA TEST OF INDEPENDENCE --- MODULATORS-SENDERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pvalue calculation for each electrodes 
pth = 0.05;
theta_NZ = 0; % counts for theta non-zero
beta_NZ = 0; % counts for beta non-zero
for ch = 1:size(coh_sr,2)
    
   stats_avg.MS.Ch(ch).theta_pval =  nnz(abs(theta_MS(ch,:)) > sum(abs(stim(ch).c_sr(15:19)))/5)/iter;
   stats_avg.MS.Ch(ch).beta_pval =  nnz(abs(beta_MS(ch,:)) > sum(abs(stim(ch).c_sr(40:44)))/5)/iter;
   
   % How many times theta and beta are significantly larger than zero (with
   % pval threshold = 0.05)
   theta_NZ = theta_NZ + int8(stats_avg.MS.Ch(ch).theta_pval <= pth); % if pval of the channel is smaller than pval threshold add 1
   beta_NZ = beta_NZ + int8(stats_avg.MS.Ch(ch).beta_pval <= pth); 
   
end

% counts of theta being significantly zero and beta being significantly zero
theta_Z = size(coh_sr,2) - theta_NZ;
beta_Z = size(coh_sr,2) - beta_NZ;

Obs = double([theta_NZ, theta_Z; beta_NZ, beta_Z])

stats_avg.MS.Obs = Obs; % save into structure 


N = sum(sum(Obs)) % Sample size
E_theta_NZ = sum(Obs(1,:)) * sum(Obs(:,1))/N;
E_theta_Z = sum(Obs(1,:)) * sum(Obs(:,2))/N;
E_beta_NZ = sum(Obs(2,:)) * sum(Obs(:,1))/N;
E_beta_Z = sum(Obs(2,:)) * sum(Obs(:,2))/N;

Exp = [E_theta_NZ,E_theta_Z; E_beta_NZ,E_beta_Z]

stats_avg.MS.Exp = Exp; % save into structure 

ch = 2;
fig = figure;
histogram(abs(theta_MS(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
hold on; histogram(abs(beta_MR(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
legend('p-val theta','p-val beta')
title('p-val distribution','FontSize',12)


O = reshape(Obs,[1,4]);
E = reshape(Exp,[1,4]);
% Compute Chi-Squared 
chi_MS = 0;
for i = 1:4
    power((O(i) - E(i)),2)/E(i);
    chi_MS = chi_MS + power((O(i) - E(i)),2)/E(i);
    
end
chi_MS


stats_avg.MS.ChiSquared = chi_MS; % save into structure 







