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
        
for s = 1:size(sess_info{1},1) % for all the sesssions with causal modulators
    
       
    Sess = sess_info{1}(s); % Session number
    dir_Sess = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d',Sess);
    dir_pval = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d/p_th_%.3f',Sess,pval);
    
    % -- load list electrodes, sender, receiver for that session 
    electrode = importdata(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairs
    receiver = importdata(strcat(dir_Sess,sprintf('/receiver_Sess_%d.txt',Sess)));  % ---- receiver pair
    sender = importdata(strcat(dir_Sess,sprintf('/sender_Sess_%d.txt',Sess))); % ---- sender pair
    
    modulator = importdata(strcat(dir_Sess,'/Modulators_idx.txt')); % import the idx of the modulator(s)
    
    AM_channel = importdata(strcat(dir_pval,sprintf('/Coherent_modulators_AM_step_%d_pval_%.3f_alpha_%.3f.txt',step,pval,alpha))); % list of coherent modulators index with abs(mean()) analysis
    MA_channel = importdata(strcat(dir_pval,sprintf('/Coherent_modulators_MA_step_%d_pval_%.3f_alpha_%.3f.txt',step,pval,alpha))); % list of coherent modulators index with mean(abs()) analysis
         
    session_AM(s).session_idx = [s, Sess]; % session i, and index
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



    end
end 







