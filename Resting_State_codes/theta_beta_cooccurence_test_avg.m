
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theta/Beta chi-squared co-occurrence test      
% with average values in an interval rather than a single value
% This code contructs the contingency matrix
%           theta Yes / theta No
%  Beta Yes
%  Beta No
%
% where yes and no mean 'different from zero' and 'not different from zero'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Perm = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state/Permutation_test';


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



% MR peaks: Theta peaks 7.8 Hz -> f(15); Beta pick 21 Hz -> f(42)
% MS peaks: Theta peaks 7.8 Hz -> f(15); Beta pick 21 Hz -> f(42)
% SR peaks: Theta peaks 8.3 Hz -> f(16); Beta pick 20.5 Hz -> f(41)

% -- modulators-receivers
theta_MR = zeros(size(coh,2),iter);
beta_MR = zeros(size(coh,2),iter);
% -- modulators-senders
theta_MS = zeros(size(coh,2),iter);
beta_MS = zeros(size(coh,2),iter);
% -- senders-receivers 
theta_SR = zeros(size(coh_sr,2),iter);
beta_SR = zeros(size(coh_sr,2),iter);


% --- select theta and beta within a range around the peaks
% Modulators-Receivers 
for ch = 1:size(coh,2)
    for perm = 1:iter
        
        theta_MR(ch,perm) = sum(abs(coh(ch).perm(perm).c_mr(13:17)))/5; % 6.8 - 8.8 Hz
        beta_MR(ch,perm) = sum(abs(coh(ch).perm(perm).c_mr(40:44)))/5;  % 20-22 Hz
        
    end
end

% Modulators-Senders 
for ch = 1:size(coh,2)
    for perm = 1:iter
        
        theta_MS(ch,perm) = sum(abs(coh(ch).perm(perm).c_ms(13:17)))/5; % 7.8 - 9.8 Hz
        beta_MS(ch,perm) = sum(abs(coh(ch).perm(perm).c_ms(40:44)))/5;  % 20 - 22 Hz
        
    end
end

% Senders-Receivers 
for ch = 1:size(coh_sr,2)
    for perm = 1:iter
        
        theta_SR(ch,perm) = sum(abs(coh_sr(ch).perm(perm).c_sr(14:18))); % 8.8 Hz
        beta_SR(ch,perm) = sum(abs(coh_sr(ch).perm(perm).c_sr(39:43)));  % 21 Hz
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THETA AND BETA TEST OF INDEPENDENCE --- MODULATORS-RECEIVERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pvalue calculation for each electrodes 
pth = 0.05;
theta_Y_beta_Y = 0; % counts for theta yes beta yes
theta_Y_beta_N = 0; % counts for theta yes beta no
theta_N_beta_Y = 0; % counts for theta no beta yes
theta_N_beta_N = 0; % counts for theta no beta no

for ch = 1:size(coh,2)
    
   stats_avg.MR.Ch(ch).theta_MR_pval =  nnz(theta_MR(ch,:) > sum(abs(mod(ch).c_mr(13:17)))/5)/iter;
   stats_avg.MR.Ch(ch).beta_MR_pval =  nnz(beta_MR(ch,:) > sum(abs(mod(ch).c_mr(40:44)))/5)/iter;
   
   % How many times theta and beta are significantly larger than zero (with
   % pval threshold = 0.05)
   theta_Y_beta_Y = theta_Y_beta_Y + int8(stats_avg.MR.Ch(ch).theta_MR_pval <= pth & stats_avg.MR.Ch(ch).beta_MR_pval <= pth); % theta Y and beta Y: if both pvalues are smaller than threshold
   theta_Y_beta_N = theta_Y_beta_N + int8(stats_avg.MR.Ch(ch).theta_MR_pval <= pth & stats_avg.MR.Ch(ch).beta_MR_pval > pth); % theta Y beta N
   theta_N_beta_Y = theta_N_beta_Y + int8(stats_avg.MR.Ch(ch).theta_MR_pval > pth & stats_avg.MR.Ch(ch).beta_MR_pval <= pth); % theta N beta Y
   theta_N_beta_N = theta_N_beta_N + int8(stats_avg.MR.Ch(ch).theta_MR_pval > pth & stats_avg.MR.Ch(ch).beta_MR_pval > pth);  % theta N beta N 
   
end


Obs = double([theta_Y_beta_Y, theta_N_beta_Y; theta_Y_beta_N, theta_N_beta_N])

stats_avg.MR.Obs = Obs; % save into structure 

N = sum(sum(Obs)) % Sample size
E_theta_Y_beta_Y = sum(Obs(1,:)) * sum(Obs(:,1))/N;
E_theta_N_beta_Y = sum(Obs(1,:)) * sum(Obs(:,2))/N;
E_theta_Y_beta_N = sum(Obs(2,:)) * sum(Obs(:,1))/N;
E_theta_N_beta_N = sum(Obs(2,:)) * sum(Obs(:,2))/N;

Exp = [E_theta_Y_beta_Y,E_theta_N_beta_Y; E_theta_Y_beta_N,E_theta_N_beta_N]

stats_avg.MR.Exp = Exp; % save into structure 

ch = 41;
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
theta_Y_beta_Y = 0; % counts for theta yes beta yes
theta_Y_beta_N = 0; % counts for theta yes beta no
theta_N_beta_Y = 0; % counts for theta no beta yes
theta_N_beta_N = 0; % counts for theta no beta no

for ch = 1:size(coh,2)
    
   stats_avg.MS.Ch(ch).theta_MS_pval =  nnz(abs(theta_MS(ch,:)) > sum(abs(mod(ch).c_ms(13:17)))/5)/iter;
   stats_avg.MS.Ch(ch).beta_MS_pval =  nnz(abs(beta_MS(ch,:)) > sum(abs(mod(ch).c_ms(40:44)))/5)/iter;
   
   % How many times theta and beta are significantly larger than zero (with
   % pval threshold = 0.05)
   theta_Y_beta_Y = theta_Y_beta_Y + int8(stats_avg.MS.Ch(ch).theta_MS_pval <= pth & stats_avg.MS.Ch(ch).beta_MS_pval <= pth); % theta Y and beta Y: if both pvalues are smaller than threshold
   theta_Y_beta_N = theta_Y_beta_N + int8(stats_avg.MS.Ch(ch).theta_MS_pval <= pth & stats_avg.MS.Ch(ch).beta_MS_pval > pth); % theta Y beta N
   theta_N_beta_Y = theta_N_beta_Y + int8(stats_avg.MS.Ch(ch).theta_MS_pval > pth & stats_avg.MS.Ch(ch).beta_MS_pval <= pth); % theta N beta Y
   theta_N_beta_N = theta_N_beta_N + int8(stats_avg.MS.Ch(ch).theta_MS_pval > pth & stats_avg.MS.Ch(ch).beta_MS_pval > pth);  % theta N beta N 
   
end


Obs = double([theta_Y_beta_Y, theta_N_beta_Y; theta_Y_beta_N, theta_N_beta_N])

stats_avg.MS.Obs = Obs; % save into structure 

N = sum(sum(Obs)) % Sample size
E_theta_Y_beta_Y = sum(Obs(1,:)) * sum(Obs(:,1))/N;
E_theta_N_beta_Y = sum(Obs(1,:)) * sum(Obs(:,2))/N;
E_theta_Y_beta_N = sum(Obs(2,:)) * sum(Obs(:,1))/N;
E_theta_N_beta_N = sum(Obs(2,:)) * sum(Obs(:,2))/N;

Exp = [E_theta_Y_beta_Y,E_theta_N_beta_Y; E_theta_Y_beta_N,E_theta_N_beta_N]

stats_avg.MS.Exp = Exp; % save into structure 

ch = 11;
fig = figure;
histogram(abs(theta_MS(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
hold on; histogram(abs(beta_MS(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THETA AND BETA TEST OF INDEPENDENCE --- SENDERS-RECEIVERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pvalue calculation for each electrodes 
pth = 0.05;
theta_Y_beta_Y = 0; % counts for theta yes beta yes
theta_Y_beta_N = 0; % counts for theta yes beta no
theta_N_beta_Y = 0; % counts for theta no beta yes
theta_N_beta_N = 0; % counts for theta no beta no

for ch = 1:size(coh_sr,2)
    
   stats_avg.SR.Ch(ch).theta_SR_pval =  nnz(abs(theta_SR(ch,:)) > sum(abs(stim(ch).c_sr(14:18)))/5)/iter;
   stats_avg.SR.Ch(ch).beta_SR_pval =  nnz(abs(beta_SR(ch,:)) > sum(abs(stim(ch).c_sr(39:43)))/5)/iter;
   
   % How many times theta and beta are significantly larger than zero (with
   % pval threshold = 0.05)
   theta_Y_beta_Y = theta_Y_beta_Y + int8(stats_avg.SR.Ch(ch).theta_SR_pval <= pth & stats_avg.SR.Ch(ch).beta_SR_pval <= pth); % theta Y and beta Y: if both pvalues are smaller than threshold
   theta_Y_beta_N = theta_Y_beta_N + int8(stats_avg.SR.Ch(ch).theta_SR_pval <= pth & stats_avg.SR.Ch(ch).beta_SR_pval > pth); % theta Y beta N
   theta_N_beta_Y = theta_N_beta_Y + int8(stats_avg.SR.Ch(ch).theta_SR_pval > pth & stats_avg.SR.Ch(ch).beta_SR_pval <= pth); % theta N beta Y
   theta_N_beta_N = theta_N_beta_N + int8(stats_avg.SR.Ch(ch).theta_SR_pval > pth & stats_avg.SR.Ch(ch).beta_SR_pval > pth);  % theta N beta N 
   
end


Obs = double([theta_Y_beta_Y, theta_N_beta_Y; theta_Y_beta_N, theta_N_beta_N])

stats_avg.SR.Obs = Obs; % save into structure 

N = sum(sum(Obs)) % Sample size
E_theta_Y_beta_Y = sum(Obs(1,:)) * sum(Obs(:,1))/N;
E_theta_N_beta_Y = sum(Obs(1,:)) * sum(Obs(:,2))/N;
E_theta_Y_beta_N = sum(Obs(2,:)) * sum(Obs(:,1))/N;
E_theta_N_beta_N = sum(Obs(2,:)) * sum(Obs(:,2))/N;

Exp = [E_theta_Y_beta_Y,E_theta_N_beta_Y; E_theta_Y_beta_N,E_theta_N_beta_N]

stats_avg.SR.Exp = Exp; % save into structure 

ch = 11;
fig = figure;
histogram(abs(theta_SR(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
hold on; histogram(abs(beta_MS(ch,:)),20,'Normalization','probability','FaceAlpha',.6); grid on
legend('p-val theta','p-val beta')
title('p-val distribution','FontSize',12)


O = reshape(Obs,[1,4]);
E = reshape(Exp,[1,4]);
% Compute Chi-Squared 
chi_SR = 0;
for i = 1:4
    power((O(i) - E(i)),2)/E(i);
    chi_SR = chi_SR + power((O(i) - E(i)),2)/E(i);
    
end
chi_SR

stats_avg.SR.ChiSquared = chi_SR; % save into structure 





