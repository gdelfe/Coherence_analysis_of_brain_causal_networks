function [zmax,id_zmax,blocks_F] = cluster_correction(coh_diff,coh_diff_perm,f,f_lim,z_th)

% CLUSTER_CORRECTION - Perform cluster correction for zscore across frequency
%  [ZMAX,ID_ZMAX, BLOCKS_F] = cluster_correction(ZSCORE,COH_DIFF,F_CLUST,F_LIM,Z_TH)
%
% Inputs: 
%     COH_DIFF = coherence(High) - coherence(Low) values vs frequency, where High and Low corresponds
%             to high and low values of the modulator score, 1D array
%     F_CLUST = frequency values array
%     F_LIM = max value for the frequency bin to use for the calculation
%     Z_TH = zscore threshold 
%     
% Outputs:
%     ZMAX = max weighted sum value of the zscores across cluster, i.e. zscore sum of the cluster 
%         with the highest sum
%     ID_MAX = id corresponding to the cluster with the highest sum
%     BLOCKS_F = range of frequencies at which z-score is clustered
%
%   Written by: Gino Del Ferraro, July 2020
%               NYU, Pesaran Lab



zscore = (coh_diff(1,:) - mean(coh_diff_perm))./std(coh_diff_perm); % compute zscore for the coherence difference 

f(1) = 1e-5;
z_clust = zscore(1:f_lim);
f_clust = f(1:f_lim);
z_clust(z_clust > -abs(z_th) & z_clust < abs(z_th)) = 0; % set to zero all the non-significant z-values, two tailed approach
f_clust(z_clust > -abs(z_th) & z_clust < abs(z_th)) = 0; % get the frequencies of the significant z-values

if nnz(z_clust ~= 0) % if some cluster has survived the thresholding
    % Group continous non-zero p-values in separate clusters
    clear blocks zClustSum
    wrap       = [0, z_clust, 0] ;
    temp       = diff( wrap ~= 0 ) ;
    blockStart = find( temp == 1 ) + 1 ;
    blockEnd   = find( temp == -1 ) ;
    blocks     = arrayfun( @(bId) wrap(blockStart(bId):blockEnd(bId)), ...
        1:numel(blockStart), 'UniformOutput', false ) ; % each block contains a series of non zero consecutive p-values
    
    % Group continous non-zero p-values in separate clusters
    clear blocks_F
    wrap_F       = [0, f_clust, 0] ;
    temp_F       = diff( wrap_F ~= 0 ) ;
    blockStart_F = find( temp_F == 1 ) + 1 ;
    blockEnd_F   = find( temp_F == -1 ) ;
    blocks_F     = arrayfun( @(bId_F) wrap_F(blockStart_F(bId_F):blockEnd_F(bId_F)), ...
        1:numel(blockStart_F), 'UniformOutput', false ) ; % each block contains a series of non zero consecutive p-values
    
    
    
    % Compute the sum of the z-scores in each cluster
    for i=1:length(blocks)
        zClustSum{i} = sum(abs(blocks{i}));
    end
    
    % Compute the max of the sum of the z-scores in each cluster
    [zmax,id_zmax] = max([zClustSum{:}]);
    
else
    zmax = 0.; % if no cluster survives the threshold, set max zscore cluster = 0.1
    id_zmax = 1;
    blocks_F = {0};
%     display(['No cluster has survived the thresholding!'])
end
    
end