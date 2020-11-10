%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates structures to store all the info/data about the
% coherence and causal modulator across session. Two different structures
% are used to store the data, AM for abs mean, and MA for mean abs
% 
% At the end of the code, the rate of the causal modulator is also
% computed, together with how many causal mod are also coherent mod 
%
% INPUT: file that contains the list of sessions with causal modulators
% OUTPUT: .mat files for the two structures AM, MA
%         file causal_modulators_rate.txt where the rate are saved
% 
%    @ Gino Del Ferraro, August 2020, Pesaran Lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%
addpath('/mnt/pesaranlab/People/Gino/DL-modulators/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State';

step = 110;

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

pval = 0.05; % pvalue threshold cluster correction --- change also the name of the path below 
alpha = 0.001; % FDR threshold 

dir_RS_pval = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/p_th_%.3f',pval)

for pval = [0.01 0.001]
    for alpha = [0.05 0.01 0.005 0.001]
        
for s = 1:size(sess_info{1},1)
    
       
    Sess = sess_info{1}(s); % Session number
    dir_Sess = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d',Sess);
    dir_pval = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d/p_th_%.3f',Sess,pval);
    
    % -- load list electrodes, sender, receiver
    electrode = importdata(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairs
    receiver = importdata(strcat(dir_Sess,sprintf('/receiver_Sess_%d.txt',Sess)));  % ---- receiver pair
    sender = importdata(strcat(dir_Sess,sprintf('/sender_Sess_%d.txt',Sess))); % ---- sender pair
    
    modulator = importdata(strcat(dir_Sess,'/Modulators_idx.txt')); % import the idx of the modulator(s)
    
    AM_channel = importdata(strcat(dir_pval,sprintf('/Coherent_modulators_AM_step_%d_pval_%.3f_alpha_%.3f.txt',step,pval,alpha)));
    MA_channel = importdata(strcat(dir_pval,sprintf('/Coherent_modulators_MA_step_%d_pval_%.3f_alpha_%.3f.txt',step,pval,alpha)));
         
    session_AM(s).session_idx = [s, Sess];
    session_MA(s).session_idx = [s, Sess];
    for i=1:length(modulator) % for all the causal modulators
        
        % --- Average Mean data
        session_AM(s).causal_mod(i) = 1; % set coherence modulator = 1
        if isempty(find(AM_channel == modulator(i))) % if the coherence mod is not a causal mod, set it to zero
            session_AM(s).causal_mod(i) = 0;
        end
        
        session_AM(s).mod_idx(i) = modulator(i); % modulator index
        session_AM(s).mod(i,:) = electrode(modulator(i),:); % modulator pair
        
        % --- Mean Average data
        session_MA(s).causal_mod(i) = 1; % set coherence modulator = 1
        if isempty(find(MA_channel == modulator(i))) % if the coherence mod is not a causal mod, set it to zero
            session_MA(s).causal_mod(i) = 0;
        end
        
        session_MA(s).mod_idx(i) = modulator(i); % modulator index
        session_MA(s).mod(i,:) = electrode(modulator(i),:); % modulator pair
    end
    
    % --- Average Mean data
    session_AM(s).cohe_mod = AM_channel';  % coherence modulator list
    session_AM(s).tot_ch = size(electrode,1); % tot number of channels
    session_AM(s).numb_sign_ch = size(AM_channel,1); % number of coherent modulators  
    session_AM(s).rate_cohe_tot = round(size(AM_channel,1)/size(electrode,1),2); % rate of coherent modulators
    session_AM(s).receiver = receiver; % receiver pair
    session_AM(s).sender = sender; % sender pair
    
    % --- Mean Average data
    session_MA(s).cohe_mod = MA_channel';  % coherence modulator list
    session_MA(s).tot_ch = size(electrode,1); % tot number of channels
    session_MA(s).numb_sign_ch = size(MA_channel,1); % number of coherent modulators
    session_MA(s).rate_cohe_tot = round(size(MA_channel,1)/size(electrode,1),2);  % rate of coherent modulators
    session_MA(s).receiver = receiver; % receiver pair
    session_MA(s).sender = sender; % sender pair
    
end

% -- Write structures files
save(strcat(dir_RS,sprintf('/session_AM_pval_%.3f_alpha_%.2f.mat',pval,alpha)),'session_AM');
save(strcat(dir_RS,sprintf('/session_MA_pval_%.3f_alpha_%.2f.mat',pval,alpha)),'session_MA');




pval = 0.05
alpha = 0.01

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



    end
end 







