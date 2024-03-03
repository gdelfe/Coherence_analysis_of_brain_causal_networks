
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the granger-causality test for sender-receiver (and viceversa),
%  sender-modulator (and viceversa), and receiver-modulator (and
%  viceversa). Perform the test for each trial. For each pair of electrode,
%  compute the rate of GC test passed (null-hypothesis not rejected) across
%  trials. Store this rate for each session, for all the pair SR, SM, RM
%  into a structure called 'gc_rate'.
%
%  The Lag value has been chosen to be 10, after trying different values
%  and based on biological motivations. 10 corresponds to 10 ms.
%
%  @ Gino Del Ferraro, March 2024, NYU, Pesaran Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

% addpath('T:/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat 
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
Lag = 10; % maxLag for the computation of GC test 

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

for i = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    dir_Mod_recording = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/%s',Sess,recording));

    if ~exist(dir_Mod_recording, 'dir')
        mkdir(dir_Mod_recording)
    end
    %     % ---  time parameter
    tot_time = 150001;
    
    outliers_SR = [sess_data_lfp.outliers_S, sess_data_lfp.outliers_R];
    outliers_SR = unique(outliers_SR)  % -- remove repeated entries in outliers
    
    % %%%%%%%%%%% Sender and Receiver LFP %%%%%%%%%%%%%%%%%%%
    lfp_S = sess_data_lfp.lfp_S;
    lfp_R = sess_data_lfp.lfp_R;
    
    
    % -- remove outliers from sender and receiver
    lfp_S(outliers_SR,:) = [];
    lfp_R(outliers_SR,:) = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- GRANGER CAUSALITY - Sender-Receiver -- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(['Computng sender-receiver GC test for session ',num2str(i)])
    tic
    gcRate_SR = granger_causality_X_Y(lfp_S, lfp_R, Lag); % sender-receeiver test
    gcRate_RS = granger_causality_X_Y(lfp_R, lfp_S, Lag); % receiver-sender test
    toc 

    % store GC for sender-receiver and receiver-sender in structure 
    gc_struct(i).SR = gcRate_SR;
    gc_struct(i).RS = gcRate_RS;
  
    
    mod_Ch = sess_data_lfp.mod_idx; % -- modulators index

    % --- Remove channels with artifacts for Maverick
    if strcmp(monkey,"Maverick")
        if Sess == 19
            mod_Ch(mod_Ch == 60) = []; % remove channel
        end
        if Sess == 41
            mod_Ch(mod_Ch == 8) = []; % remove channel
        end
    end
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
    
    % %%%%%%% ALL Electrodes LFP %%%%%%%%%%%%%%%%%%%%%
    lfp_E_all = sess_data_lfp.lfp_E;
    
    cnt_m = 1;
    for Ch = mod_Ch % for all the modulators in the session
        
        close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- GRANGER CAUSALITY --- Modulator -Sender(-Receiver) -- %
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
            % -- outliers
            outliers_ES = [outliers_E, outliers_S];
            outliers_ES = unique(outliers_ES);
            % -- remove outliers from sender and modulator 
            lfp_S(outliers_ES,:) = [];
            lfp_E(outliers_ES,:) = [];
            
            disp(['Computng sender-modulator GC test. Modulator ',num2str(cnt_m),' Out of ',num2str(length(mod_Ch))])
            tic % compute granger-test for SM and MS
            gcRate_SM = granger_causality_X_Y(lfp_S, lfp_E, Lag); % sender-modulator test
            gcRate_MS = granger_causality_X_Y(lfp_E, lfp_S, Lag); % modulator-sender test
            toc 
            
            % store GC for sender-modulator and modlator-sender in structure
            gc_struct(i).mod(cnt_m).SM = gcRate_SM; 
            gc_struct(i).mod(cnt_m).MS = gcRate_MS;

            
            % -- MODULATOR - RECEIVER
            % -- modulaor lfp
            lfp_E = sq(lfp_E_all(Ch,:,:));          % -- get lfp for only that channel
            outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;  % -- get the modulator's outliers
            % -- Receiver  LFP
            lfp_R = sess_data_lfp.lfp_R;
            outliers_R = sess_data_lfp.outliers_R;
            % -- outliers
            outliers_ER = [outliers_E, outliers_R];
            outliers_ER = unique(outliers_ER);
            % -- remove outliers from sender and modulator 
            lfp_R(outliers_ER,:) = [];
            lfp_E(outliers_ER,:) = [];
          
            disp(['Computng receiver-modulator GC test. Modulator ',num2str(cnt_m),' Out of ',num2str(length(mod_Ch))])
            tic % compute granger-test for RM and MR
            gcRate_RM = granger_causality_X_Y(lfp_R, lfp_E, Lag); % receiver-modulator test
            gcRate_MR = granger_causality_X_Y(lfp_E, lfp_R, Lag); % modulator-receiver test
            toc 
            
            % store GC for receiver-modulator and modlator-sender in structure
            gc_struct(i).mod(cnt_m).RM = gcRate_RM; 
            gc_struct(i).mod(cnt_m).MR = gcRate_MR;
            
        end                     
        cnt_m = cnt_m + 1; % counter for modulators within this session
    end
end


dir_GC_results = strcat(dir_RS_Theta,sprintf('/GC_results/'));
if ~exist(dir_GC_results, 'dir')
    mkdir(dir_GC_results)
end    


% Save coherence and spectrum data in structure format
save(strcat(dir_GC_results,sprintf('GC_rate_%d.mat',Lag)),'gc_struct');

% load(strcat(dir_GC_results,sprintf('GC_rate_%d.mat',Lag)))


