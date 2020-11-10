%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code runs FDR multiple comparison across coherent modulators in the
% same Session
%
% INPUTS: Cluster corrected p values across frequency for each channel
% OUTPUT: FDR plots, list of channel that are significant under the FDR
% correction, corresponding pvalues
%
% Gino Del Ferraro, August 2020, Pesaran Lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all;

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%
addpath('/mnt/pesaranlab/People/Gino/DL-modulators/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State';
step = 110;
pval = 0.001 % pval cluster correction --- change the path of the folder in the for loop, Sess/p_th_X

for pval = [0.05, 0.01,0.001]
for alpha = [0.05, 0.01, 0.005, 0.001]  % alpha = 0.05; % confidence level for FDR

close all

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

for s = 1:size(sess_info{1},1)  % For all the sessions with at least one modulator 
    
    Sess = sess_info{1}(s); % Session number
    dir_Sess_pval = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d/p_th_%.3f',Sess,pval);
    
    p_channel_AM = importdata(strcat(dir_Sess_pval,sprintf('/p_channel_AM_step_%d.txt',step))); % cluser corrected pvalues for each channel 
    p_channel_MA = importdata(strcat(dir_Sess_pval,sprintf('/p_channel_MA_step_%d.txt',step)));
    
    p_channel_AM = p_channel_AM';
    p_channel_MA = p_channel_MA';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MULTIPLE COMPARISONS  %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -------------------------------- %
    % --- FDR Benjamini-Hochberg test  %
    % -------------------------------- %
    
    
    % -- sort p-value in ascending order
    [sort_p_AM,idx_p_AM] = sort(p_channel_AM);
    [sort_p_MA,idx_p_MA] = sort(p_channel_MA);
    
    m = size(p_channel_AM,2); % total number of p-values
    
     
    
    k_th = (1:m)*alpha/m; % FDR threshol p-values
        
    idx_sign_AM = sort_p_AM <= k_th; % get the index of the p-values below FDR threshold
    th_idx_AM = strfind(idx_sign_AM,[1 0]); % index of the highest significant p-value
    if isempty(th_idx_AM) % if all the FDR p-values are significant ->
        th_idx_AM  =  length(idx_p_AM); % -> take all of them
    end
    Ch_sign_AM = idx_p_AM(1:th_idx_AM); % channels that show significance after FDR correction, with alpha confidence interval
    if nnz(~idx_sign_AM) == length(idx_sign_AM) % if none of the pvalue is less than FDR thresh, set the coherence mod to empty
    Ch_sign_AM = [];
    end 
    
    idx_sign_MA = sort_p_MA <= k_th; % threshold p-values below FDR threshold
    th_idx_MA = strfind(idx_sign_MA,[1 0]); % index of the highest significant p-value
    if isempty(th_idx_MA) % if all the FDR p-values are significant ->
        th_idx_MA  =  length(idx_p_MA); % -> take all of them
    end
    Ch_sign_MA = idx_p_MA(1:th_idx_MA); % channels that show significance after FDR correction, with alpha confidence interval
    if nnz(~idx_sign_MA) == length(idx_sign_MA) % if none of the pvalue is less than FDR thresh, set the coherence mod to empty
        Ch_sign_MA = [];
    end
    
    
    dlmwrite(strcat(dir_Sess_pval,sprintf('/Coherent_modulators_AM_step_%d_pval_%.3f_alpha_%.3f.txt',step,pval,alpha)),Ch_sign_AM'); % coherent modulators
    dlmwrite(strcat(dir_Sess_pval,sprintf('/Coherent_modulators_MA_step_%d_pval_%.3f_alpha_%.3f.txt',step,pval,alpha)),Ch_sign_MA');
    
    
    % -- FIGURE: FDR plot
    fig = figure;
    plot(1:m,k_th)
    hold on
    plot(1:m,sort_p_AM)
    hold on
    plot(1:m,sort_p_MA)
    grid on
    ylabel('p-values')
    xlabel('k')
    legend('k*alpha/m','p-values AM','p-values MA','Location','northwest')
    title(sprintf('FDR BH test, Sess = %d, pval = %.3f alpha= %.3f',Sess,pval,alpha),'FontSize',11)
    
    fname = strcat(dir_Sess_pval,sprintf('/FDR_plot_step_%d_pval_%.3f_alpha_%.3f.jpg',step,pval,alpha));
    saveas(fig,fname);
    
end
end
end
