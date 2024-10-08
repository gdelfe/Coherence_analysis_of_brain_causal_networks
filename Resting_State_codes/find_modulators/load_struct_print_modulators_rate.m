%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code loads the structures for a given p-value and alpha FDR
% threshold and computes the rate of coherent and causal modulators for
% that specific pval and alpha
%
% INPUT: - file that contains the list of sessions with causal modulators
%        - structure files AM and MA in .mat for a given pval and alpha FDR 
% OUTPUT: A summary structure which contains info about the rates of
%         coherent and causal modulators for a given pval and alpha across 
%         all the sessions
% 
%    @ Gino Del Ferraro, August 2020, Pesaran Lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%
addpath('/vol/bd5/People/Gino/Gino_codes')
dir_RS = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_State';

step = 110;

pval = 0.05; % pvalue threshold cluster correction --- change also the name of the path below 
alpha = 0.01; % FDR threshold 

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


dir_RS_pval = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/p_th_%.3f',pval)


% -- load structures files 
newAM = load(strcat(dir_RS,sprintf('/session_AM_pval_%.3f_alpha_%.2f.mat',pval,alpha)));
session_AM = newAM.session_AM;

newMA = load(strcat(dir_RS,sprintf('/session_MA_pval_%.3f_alpha_%.2f.mat',pval,alpha)));
session_MA = newMA.session_MA;
% 

% -- print structures on stdout 
format short
for s=1:size(sess_info{1},1)
    session_AM(s);
    session_MA(s);
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      COMPUTE COUNT AND RATE OF CAUSAL MODULATOR ACROSS SESSIONS     %
%          and how many coherent modulators are causal modulators     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tot_mod_AM = 0; tot_mod_MA = 0; % causal modulators in stim 
cnt_AM = 0; cnt_MA = 0; % correctly identified causal mod 
rate_cohe_AM = 0; rate_cohe_MA = 0;  % rate correctly identified mod 

cnt_mod = 0; % numb of modulators in total 

for s=1:size(sess_info{1},1) % for all the sessions with modulators 
  
    cnt_AM = cnt_AM + nnz(session_AM(s).causal_mod); % count of how many causal modulators identified correctly
    tot_mod_AM = tot_mod_AM + length(session_AM(s).causal_mod); % count of how many  causal modulator from stim exp
    rate_cohe_AM = rate_cohe_AM + session_AM(s).rate_cohe_tot; % rate of how many correcly identified modulators
    
    cnt_MA = cnt_MA + nnz(session_MA(s).causal_mod);
    tot_mod_MA = tot_mod_MA + length(session_MA(s).causal_mod);
    rate_cohe_MA = rate_cohe_MA + session_MA(s).rate_cohe_tot;
    
    cnt_mod = cnt_mod + length(session_AM(s).mod_idx); % number of modulators in a session
    
end

rate_AM = round(cnt_AM/tot_mod_AM,2); % rate AM: number of causal modulator decoded correctly/numb of causal modulator in stim exp
rate_MA = round(cnt_MA/tot_mod_MA,2); % rate MA: number of causal modulator decoded correctly/numb of causal modulator in stim exp

rate_cohe_AM = rate_cohe_AM/size(sess_info{1},1); % average rate of coherent modulator across section 
rate_cohe_MA = rate_cohe_MA/size(sess_info{1},1);

UsedSess = importdata(strcat(dir_RS,'/Sessions_list.txt'));
tot_ch = 0;
for Sess = UsedSess
    dir_Sess = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d',Sess);
    electrode = importdata(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairsend
    tot_ch = tot_ch + size(electrode,1) ;
end

rate_mod = cnt_mod/tot_ch; % rate of how many modulator across all sections 


summary.rate_CausCorrectDecoded_AM = rate_AM; % rate of correclty decoded causal modulators AM
summary.rate_CausCorrectDecoded_MA = rate_MA; % rate of correclty decoded causal modulators MA
summary.avg_rate_coherent_AM = rate_cohe_AM; % avg rate of coherent modulators across sections AM
summary.avg_rate_coherent_MA = rate_cohe_MA; % avg rate of coherent modulators across sections MA
summary.rate_tot_causal = rate_mod; % rate of tot causal modulators across sections 

summary % print 

save(strcat(dir_RS_pval,sprintf('/summary_pval_%.3f_alpha_%.3f.mat',pval,alpha)),'summary');

% new = load(strcat(dir_RS,sprintf('/summary_pval_%.3f_alpha_%.2f.mat',pval,alpha)))
% summary = new.summary;






