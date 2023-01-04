
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the permutation test for the MS and MR coherence. 
% It creates the null distribution for such test by computing null
% coherence: the latter is computed by randomly permuting the lfp time
% points of the modulator and recomputing the coherence between the
% lfp-permuted and the sender (receiver) non-permuted lfp. Results are
% stored session by session. N iterations for each session = 1000
%
%    @ Gino Del Ferraro, Jan 2023, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

iterations = 1000;
% ---- parameters for the coherence-gram
tot_time = 150001;
nt = tot_time;
fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;
% --- coherence

for i = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    mod_Ch = sess_data_lfp.mod_idx; % -- modulators (not controls!) index
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
    
    % %%%%%%% ALL Electrodes LFP %%%%%%%%%%%%%%%%%%%%%
    lfp_E_all = sess_data_lfp.lfp_E;
    
    cnt_m = 1;
    for Ch = mod_Ch % for all the modulators in the session
        
        close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- COHERENCE- Modulator - Sender/Receiver -- %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if Ch ~= sess_data_lfp.receiver_idx % if the electrode is not the receiver itself
            
            % -- remove outliers from modulator, sender, and receiver
            
            % -- MODULATOR - SENDER
            % -- modulaor lfp
            lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
            outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;  % -- get the modulator's outliers
            % -- Sender  LFP
            lfp_S = sess_data_lfp.lfp_S;
            outliers_S = sess_data_lfp.outliers_S;
            % -- Receiver  LFP
            lfp_R = sess_data_lfp.lfp_R;
            outliers_R = sess_data_lfp.outliers_R;
            
            % -- outliers
            outliers_ESR = [outliers_E, outliers_S, outliers_R];
            outliers_ESR = unique(outliers_ESR);
            % -- remove outliers from sender and modulator
            lfp_S(outliers_ESR,:) = [];
            lfp_R(outliers_ESR,:) = [];
            lfp_E(outliers_ESR,:) = [];
            
            dim = size(lfp_E,2);
            
            c_mr_permute = [];
            c_ms_permute = [];
            display(['Computing permuted coherence for modulator ',num2str(cnt_m),' out of ',num2str(length(mod_Ch))])

            for iter = 1:iterations
                
                perm = randperm(dim);
                lfp_E_perm = lfp_E(:,perm); % permute time points in lfp_E time series
                
                % -- coherence for modulator-sender, modulator-receiver
                [c_ms,f] = coherency(lfp_E_perm,lfp_S,[N W],fs,fk,pad,0.05,1,1);
                [c_mr,f] = coherency(lfp_E_perm,lfp_R,[N W],fs,fk,pad,0.05,1,1);
                
                c_ms_permute = [c_ms_permute; c_ms];
                c_mr_permute = [c_mr_permute; c_mr];
                
            end
            
            coh_perm.mod(cnt_m).c_ms_perm = c_ms_permute;
            coh_perm.mod(cnt_m).c_mr_perm = c_mr_permute;
            cnt_m = cnt_m + 1;
            
            dir_Mod_permute = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/permutations',Sess));
            if ~exist(dir_Mod_permute, 'dir')
                mkdir(dir_Mod_permute)
            end
            
            save(strcat(dir_Mod_permute,'/coherence_MS_MR_permuted_single.mat'),'coh_perm');
            clear coh_perm
            
        end % if electrode is not modulator
    end % for all modulators
end % for all 
        
        
