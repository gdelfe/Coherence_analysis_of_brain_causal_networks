
close all; clear all;

dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';

lfp_S = importdata(strcat(dir_main,'/lfp_S.txt'));
lfp_R = importdata(strcat(dir_main,'/lfp_R.txt'));

display(['Computing coherence-gram...'])
% ---- parameters for the coherence-gram
tapers = [4 4];
N = tapers(1);
dn = 0.01;
fs = 1000;


% --- coherence-gram --- UNCOMMENT if you want to recalculate

% [c_sr,tf,f,spec_r,spec_s] = tfcoh_GINO(lfp_S,lfp_R,tapers,1e3,dn,60,2,[],[],1);
% 
% data.c_sr = c_sr;
% data.f = f;
% data.tf = tf;
% data.spec_r = spec_r;
% data.spec_s = spec_s;
% data.H = Hix;
% data.L = Lix;


% save(strcat(dir_main,'/data_structure.mat'),'data');

load(strcat(dir_main,'/data_structure.mat'));
f = data.f;
tf = data.tf;
c_sr = data.c_sr;
H = data.H;
L = data.L;


% -- Figure: coherence spectrum
fig = figure; tvimage(abs(c_sr(:,:))); colorbar; % coherence spectrum

xticks = floor(linspace(1,length(tf),5));
xticklabels = tf(xticks);
xtickformat('%d')
% yticks = floor(linspace(1,length(f),10));
yticks = 1:50:length(f);

yticklabels = floor(f(yticks));
ytickformat('%.2f')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
title('Coherence-gram S-R','FontSize',12);
xlabel('time (sec)');
ylabel('freq (Hz)')
% ylim([0,120])
set(gcf, 'Position',  [100, 600, 1000, 600])

fname = strcat(dir_main,'/coherence-gram_SR.jpg');
saveas(fig,fname);




% ----- FIGURE: Difference: coh(high) - coh(low)
fig = figure;
plot(f,mean(abs(c_sr(H,:))))
hold on 
plot(f,mean(abs(c_sr(L,:))))
title('Coherence H vs L')
legend('choerence H','coherence L')
ylabel('coherence');
xlabel('freq (Hz)')
grid on
xlim([0 50])

fname = strcat(dir_main,'/coherence_H_L.jpg');
saveas(fig,fname);

% -- Compute coherence differences 
result_MA = mean(abs(c_sr(H,:))) - mean(abs(c_sr(L,:))); 

% ----- FIGURE: Difference: coh(high) - coh(low)
fig = figure;
plot(f,result_MA)
title('Coh diff: H val - L val')
legend('diff choerence')
ylabel('coherence');
xlabel('freq (Hz)')
grid on
xlim([0 50])

fname = strcat(dir_main,'/coherence_difference.jpg');
saveas(fig,fname);

        
% ------------------------  %%
% -- PERMUTATION TEST ----  %%
% ------------------------- %%

display(['Running permutation test...'])
iter = 1000; % number of iterations
diff_MA = zeros(iter,size(c_sr,2));
for j = 1:iter
    
    perm = randperm(length(tf));
    perm_H = perm(1:length(H));
    perm_L = perm(length(H)+1:length(L)+length(H)); 
    
    
    c_H_meanAbs = mean(abs(c_sr(perm_H,:)),1);
    c_L_meanAbs = mean(abs(c_sr(perm_L,:)),1);
    
    diff_meanAbs = c_H_meanAbs - c_L_meanAbs; % difference between H and L coherence for the permuted data
    
    diff_MA(j,:) = diff_meanAbs; % -- store permuted coherence difference in a matrix 
    
end
display(['Permutation test DONE'])



% ----- FIGURE: Pemuted Difference: coh(high) - coh(low)
fig = figure;
plot(f,diff_MA(10:15,:))
title('Coh diff permuted: H val - L val')
legend('permuted diff choerence')
ylabel('coherence');
xlabel('freq (Hz)')
grid on
xlim([0 50])

fname = strcat(dir_main,'/coherence_perm_10.jpg');
saveas(fig,fname);




% --- FIGURE: mean of the differences across samples vs frequency
fig = figure;
plot(f,mean(diff_MA));
legend('diff coherence mean - permutated data')
title('mean of the diff coh(high) - coh(low), perm data', 'FontSize',10);
grid on


% --- FIGURE: std of the differences across samples vs frequency
fig = figure;
plot(f,std(diff_MA));
title('Std of the diff coh(high) - coh(low), perm data', 'FontSize',10);
legend('std coherence perm data')
grid on

% Histogram of diff coherence at a given frequency 
freq = 60;
fig = figure; histogram(diff_MA(:,60),30,'Normalization','probability','FaceAlpha',.6); grid on
legend('diff coherence')
title('Difference at a given frequency','FontSize',12)




% ------------------------  %%
% -----   Z-SCORES     ----  %%
% ------------------------- %%

zscore_MA = (result_MA - mean(diff_MA))./std(diff_MA);

% -- FIGURE: zscores vs frequency
fig = figure;
plot(f,zscore_MA)
hold on
title('zscore vs frequency','FontSize',11)
legend('zscore MA')
ylabel('z-score')
xlabel('frequency')
xlim([0,30])
grid on


fname = strcat(dir_main,'/zscore_vs_freq.jpg');
saveas(fig,fname);

% ------------------------  %%
% -----   P-VALUES     ---  %%
% ------------------------- %%

pval_MA = 2*normcdf(-abs(zscore_MA));


%-- FIGURE: p-value vs frequency, assuming gaussian
fig = figure;
semilogy(f,pval_MA)
grid on
title('p-value vs frequency, assuming gaussian','FontSize',10)
legend('p-val','Location','southeast')
xlim([0,30])

fname = strcat(dir_main,'/pval_vs_freq.jpg');
saveas(fig,fname);

% -- FIGURE: histogram of p-values across frequencies
fig = figure;
histogram(pval_MA,20,'Normalization','probability','FaceAlpha',.6); grid on
legend('p-val MA')
title('p-val distribution','FontSize',12)

fname = strcat(dir_main,'/pval_histo.jpg');
saveas(fig,fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLUSTER CORRECTION  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p_th = 5e-02; % p-value threshold
z_th =norminv((0.5*p_th)) % get zscore threshold from p-val treshold


f_lim = 246; % f bin limit: -> consider only frequencies up to 30 Hz

[zmax_MA_Obs,id_zmax_MA_Obs,blocks_F_MA] = cluster_correction(result_MA,diff_MA,f,f_lim,z_th); % Mean Abs


% ------------------------ %
% --- PERMUTED VALUES ---- %
% ------------------------ %

zmaxMA_Perm = zeros(1,iter);
zscore_perm = zeros(iter,size(f,2));

for i=1:iter
    [zmaxMA_Perm(i),id_zmax] = cluster_correction(diff_MA(i,:),diff_MA,f,f_lim,z_th);
    zscore_perm(i,:) = (diff_MA(i,:) - mean(diff_MA))./std(diff_MA);
end


% -- FIGURE: permuted zscores vs frequency
fig = figure;
plot(f,zscore_perm(2,:))
hold on
title('permuted zscore vs frequency','FontSize',11)
legend('zscore MA')
ylabel('z-score')
xlabel('frequency')
xlim([0,30])
grid on

fname = strcat(dir_main,'/perm_zscore_2.jpg');
saveas(fig,fname);

% Histogram of the z-score max
fig = figure;
histogram(zmaxMA_Perm,30,'Normalization','probability','FaceAlpha',.6); grid on
legend('z-score cluster size MA')
title('z-score cluster size','FontSize',12)


fname = strcat(dir_main,'/zscore_max_histo.jpg');
saveas(fig,fname);

% ------------------------  %%
% -----   P-VALUES     ---  %%
% ------------------------- %%

% -- compute the cluster p-value from the histogram of zscore max
pClust_MA = nnz(zmaxMA_Perm > abs(zmax_MA_Obs) | zmaxMA_Perm < -abs(zmaxMA_Perm))/iter;


if pClust_MA == 0
    pClust_MA = 1./iter;
end

% -- store the cluster p-value for each channel
p_channel_MA(Ch) = pClust_MA; % ---- do not exectute this line 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTIPLE COMPARISONS  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------- %
% --- FDR Benjamini-Hochberg test  %
% -------------------------------- %

dir_pvalue = strcat(dir_RS,'Sess_2/p_th_0.005/');
% p_channel_MA = importdata = strcat(dir_pvalue,'p_channel_MA_step_110.txt')


% -- sort p-value in ascending order
[sort_p_MA,idx_p_MA] = sort(p_channel_MA);

alpha = 0.05; % confidence level for FDR
m = size(p_channel_MA,2); % total number of p-values

k_th = (1:m)*alpha/m; % FDR threshol p-values

idx_sign_MA = sort_p_MA <= k_th; % threshold p-values below FDR threshold
th_idx_MA = strfind(idx_sign_MA,[1 0]); % index of the highest significant p-value
Ch_sign_MA = idx_p_MA(1:th_idx_MA); % channels that show significance after FDR correction, with alpha confidence interval

% -- FIGURE: FDR plot
fig = figure;
plot(1:m,k_th)
hold on
plot(1:m,sort_p_MA)
grid on
ylabel('p-values')
xlabel('k')
legend('k*alpha/m','p-values MA','Location','northwest')
title('FDR Benjamini-Hochberg test, alpha=0.05','FontSize',11)






