
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the STIM coherence between the causal modulators found by
% Shaoyu's and the receiver
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators and plots this coherence for the
% hits and the misses separately
%
% INPUT: file with session modulator info
%        .mat file with structure AM and MA information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

subjects = {'maverick','archie'};

%for
iSubject = 1% : length(subjects) % Loop on the animals
%     clearvars -except subjects iSubject
if strcmp(subjects{iSubject},'archie')
    archie_vSUBNETS220_rig3
else
    maverick_vSUBNETS220_rig3
end
PreStimSess = PreStimResponseAll_Database_NetworkEdge;

iter = 1000;

for iSess = 1 : numel(PreStimSess)
    %         pAccLLRperm = PreStimSess{iSess}{end-2};
    %         useSessIndx(iSess) = pAccLLRperm <= 0.05 & ~isnan(pAccLLRperm);
    
    pAccLLRperm = PreStimSess{iSess}{end-2};
    pFDR_logic = PreStimSess{iSess}{end-1};
    useSessIndx(iSess) = pFDR_logic & ~isnan(pAccLLRperm);
end

UsedSess = find(useSessIndx);

%%%%%%%%%%%%%%%%%%%
% - SET PATHS --- %
%%%%%%%%%%%%%%%%%%%


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
dir_Stim = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data';
step = 110;

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


cnt_el = 1;
list_sess = 1:19;
list_sess(17) = []; % -- Session 17 and 20 are full of artifacts

for i = list_sess % size(sess_info{1},1)  % For all the session with a modulator
    
    Sess = sess_info{1}(i); % Session number
    
    disp(['Session ' num2str(Sess) ' out of ' num2str(length(PreStimSess)) ' ...'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD LFP STIM ...                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bn_Pre = [-1005 -5]; % ms
    [Lfp_Pre, Data] = load_LFP_STIM(PreStimSess,Sess,bn_Pre,DATADIR);
    
    % %%%%%%%%%%%%%%%%% --->>> LFP data loaded at this point
    
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',   out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    
    % ---- bipolar referencing, pairs of electrodes
    dir_Sess = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Stim_data/Sess_%d',Sess);
    dir_RS_Sess = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state/Sess_%d',Sess);
    
    if ~exist(dir_Sess, 'dir')
        mkdir(dir_Sess)
    end
    
    load(strcat(dir_RS_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    
    % -- load list electrodes, sender, receiver
    electrode = sess_data.RecordPair; % ---- all electrode pairs
    receiver = sess_data.receiver_pair;  % ---- receiver pair
    sender = sess_data.sender_pair; % ---- sender pair
    
    % ---  freq parameter for the masking
    fmin = 10;
    fmax = 40;
    
    %% %%%%%% Assign LFP %%%%%%
    lfp_E = sq(Lfp_Pre(:,electrode(:,1),:) - Lfp_Pre(:,electrode(:,2),:)); % modulator lfp
    lfp_R = sq(Lfp_Pre(:,receiver(1),:) - Lfp_Pre(:,receiver(2),:)); % receiver lfp
    lfp_S = sq(Lfp_Pre(:,sender(1),:) - Lfp_Pre(:,sender(2),:)); % sender lfp
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Coherency LFP                  %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ---- parameters for the coherence-gram
    nt = 1000;
    fs = 1000;
    fk = 200;
    pad = 2;
    N = 1;
    W = 5;
    
    
    mod_Ch = sess_data.mod_idx; % -- modulators (not controls!) index
    
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
    
    dir_Modulators = strcat(dir_Sess,'/Modulators');
    if ~exist(dir_Modulators, 'dir')
        mkdir(dir_Modulators)
    end
    
        
    for Ch = mod_Ch % for all the modulators in the session
        
        close all
        
        if Ch ~= sess_data.receiver_idx % if the electrode is not the receiver itself
            
                                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -- Permutation test for the coherence MR, MR
             
            % -- coherence of modulator-receiver
            display(['Computing modulator-receiver coherence...'])
            for j=1:iter 
                
                perm = randperm(size(lfp_R,1));
%                 tic
                [c_mr,f,S_m,S_r] = coherency(sq(lfp_E(perm,Ch,:)),lfp_R,[N W],fs,fk,pad,0.05,1,11);
%                 toc                
                
                % -- structure assignements
                stim(cnt_el).perm(j).c_mr = c_mr ; % MR coherence
                stim(cnt_el).perm(j).s_m = S_m;  % M spectrum
                stim(cnt_el).perm(j).s_r = S_r; % R spectrum
            end
              

            
            cnt_el = cnt_el + 1;  % -- count # of electrodes (modulators)
            
        end
    end % --- end of all modulator channels
    
end

 
dir_Perm = strcat(dir_Stim,'/Permutation_test');
if ~exist(dir_Perm, 'dir')
    mkdir(dir_Perm)
end

save(strcat(dir_Perm,sprintf('/coh_spec_mr_permuted_fk_%d_W_%d.mat',fk,W)),'stim');

