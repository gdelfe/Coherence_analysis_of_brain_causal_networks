%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the STIM coherence between the causal modulators found by
% Shaoyu's and the receiver
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators and plots this coherence for the
% hits and the misses separately
%
% IMPORTANT: For each channel it uses a variable number of trials in order
% to have always STIM and RS with the same amount of trial for equal
% statistics. It achieves so by checking which between STIM and RS has less
% trials and choses that as number of trials for both 
% Random permutation of trial is applied when sampling in order to sample
% homogeneously 
%
% INPUT: file with session modulator info
%        .mat file with structure AM and MA information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all
set(0,'DefaultFigureVisible','off')
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

tot_trials_RS = 150; % number of total trial in the Resting state
n_trials_STIM = 110; % number of total trial in the STIM 
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
    
    % Resting state directory 
    dir_RS_Sess = sprintf('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state/Sess_%d',Sess);
    load(strcat(dir_RS_Sess,'/sess_data_lfp.mat')); % RS LFP data split into 1 sec window and artifacts annotated 
       
    % -- load list electrodes, sender, receiver
    electrode = sess_data.RecordPair; % ---- all electrode pairs
    receiver = sess_data.receiver_pair;  % ---- receiver pair
    sender = sess_data.sender_pair; % ---- sender pair
    
    % ---  freq parameter for the masking
    fmin = 10;
    fmax = 40;
    
    % %%%%%% Assign LFP %%%%%%
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
    
    dir_Modulators = strcat(dir_Sess,'/Modulators_n_trials_as_RS');
    if ~exist(dir_Modulators, 'dir')
        mkdir(dir_Modulators)
    end
    
    
    indx_list = [];
    cnt_m = 1;
    for Ch = mod_Ch % for all the modulators in the session
        
        close all
        
        if Ch ~= sess_data.receiver_idx % if the electrode is not the receiver itself
                        
            indx_list = [indx_list, cnt_el]; % store the cnt number --- needed for multiple plotting
            
            outliers_RS = size(sess_data_lfp.outliers_tot(cnt_m).idx,2); % number of artifacts in RS LFP
            n_trials_RS = (tot_trials_RS - outliers_RS);   % # of trials without artifacts RS 
            
            % --- Select STIM trials to match the number of RS trial for idential statistics
            n_trials = min(n_trials_RS,size(lfp_R,1)) % select the min numb of trial between RS and STIM 
                                      
            perm = randperm(size(lfp_R,1)); % Randomly permute the STIM trials to sample homogeneously
            trials = perm(1:n_trials); % get STIM trials: as many as RS'

            display(['-- Mod ',num2str(Ch),' -- N trials: ',num2str(size(trials,2))])
            
            % %%%% Reduce the number of trails in order to match the Resting State
            lfp_E = lfp_E(trials,:,:);
            lfp_R = lfp_R(trials,:);
            lfp_S = lfp_S(trials,:);
            
                       
            % -- labels Hits and Misses
            hitIndx = Data.spec.lfp.DetectedIndx{Ch}; % labels for the hits (which trial was a hit)
            missIndx = Data.spec.lfp.notDetectedIndx{Ch}; % labels for the misses (which trial was a miss)
            
            % -- coherence of modulator-receiver
            display(['Computing modulator-receiver coherence...'])
            tic
            [c_mr,f,S_m,S_r] = coherency(sq(lfp_E(:,Ch,:)),lfp_R,[N W],fs,fk,pad,0.05,1,11);
            toc
            % NOTE: since we are including only a fraction of the whole
            % trials. Hits and Miss trails need to be taken among these set
            % -- coherence of modulator-receiver HITS
            [C,tr_hits,i_hits] = intersect(trials,hitIndx); % -- get the trials that are both hitIndx and in trials 
            [c_mr_H,f_H,S_m_H,S_r_H] = coherency(sq(lfp_E(tr_hits,Ch,:)),lfp_R(tr_hits,:),[N W],fs,fk,pad,0.05,1,11);
            % -- coherence of modulator-receiver MISSES
            [C,tr_miss,i_miss] = intersect(trials,missIndx); % -- get the trials that are both missIndx and in trials 
            [c_mr_M,f_M,S_m_M,S_r_M] = coherency(sq(lfp_E(tr_miss,Ch,:)),lfp_R(tr_miss,:),[N W],fs,fk,pad,0.05,1,11);
            
            
            % --- FIGURE --------- %%
            % -- Coherence vs frequency --- %
            fig = figure;
            plot(f,abs(c_mr),'color',[0, 15, 26]/255)
            hold on
            plot(f,abs(c_mr_H),'color',[0, 153, 255]/255)
            hold on
            plot(f,abs(c_mr_M),'color',[255, 153, 51]/255)
            grid on
            title(sprintf('STIM: Abs coherence vs frequency, ch = %d, mod',Ch),'FontSize',10);
            legend('M-R coherence','M-R coherence HITS','M-R coherence MISSES')
            %         xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Modulators,sprintf('/coherency_ch_%d_N_%.2f_W_%d.jpg',Ch,N,W));
            saveas(fig,fname);
            
            
            % --- CALCULATION OF THE SPECTRUM USING DMTSPEC --- %%%%%%%%
            % -- ALL
            dmt_S_r = double(dmtspec(lfp_R,[1 3],fs,fk, pad, 0.05, 1));
            dmt_S_m = double(dmtspec(sq(lfp_E(:,Ch,:)),[1 3],fs,fk, pad, 0.05, 1));
            % -- HITS
            dmt_S_r_H = double(dmtspec(lfp_R(tr_hits,:),[1 3],fs,fk, pad, 0.05, 1));
            dmt_S_m_H = double(dmtspec(sq(lfp_E(tr_hits,Ch,:)),[1 3],fs,fk, pad, 0.05, 1));
            % -- MISSES 
            dmt_S_r_M = double(dmtspec(lfp_R(tr_miss,:),[1 3],fs,fk, pad, 0.05, 1));
            dmt_S_m_M = double(dmtspec(sq(lfp_E(tr_miss,Ch,:)),[1 3],fs,fk, pad, 0.05, 1));
            
           
            % --- FIGURE --------- %%
            % -- Spectrum vs frequency --- %
            fig = figure;
            hAx=axes;
            hAx.XScale='linear'
            hAx.YScale='log'
            plot(f,abs(S_m_H),'color',[0, 51, 204]/255)
            hold on
            plot(f,abs(S_r_H),'color',[0, 153, 255]/255)
            hold on
            plot(f,abs(S_m_M),'color',[255, 153, 51]/255)
            grid on
            plot(f,abs(S_r_M),'color',[255, 51, 0]/255)
            grid on
            title(sprintf('STIM: Abs coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
            legend('M Spectrum Hits','R Spectrum Hits','M Spectrum Misses','R Spectrum Misses')
            xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Modulators,sprintf('/spectrum_ch_%d_N_%.2f_W_%d.jpg',Ch,N,W));
            saveas(fig,fname);
            
            
             % --- FIGURE --------- %%
            % -- DMT Spectrum vs frequency --- %
            fig = figure;
            hAx=axes;
            hAx.XScale='linear'
            hAx.YScale='log'
            plot(f,abs(dmt_S_m_H),'color',[0, 51, 204]/255)
            hold on
            plot(f,abs(dmt_S_r_H),'color',[0, 153, 255]/255)
            hold on
            plot(f,abs(dmt_S_m_M),'color',[255, 153, 51]/255)
            grid on
            plot(f,abs(dmt_S_r_M),'color',[255, 51, 0]/255)
            grid on
            title(sprintf('STIM: Abs coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
            legend('M Spectrum Hits','R Spectrum Hits','M Spectrum Misses','R Spectrum Misses')
            xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Modulators,sprintf('/spectrum_dmt_ch_%d_N_%.2f_W_%d.jpg',Ch,N,W));
            saveas(fig,fname);
            
            
            % -- structure assignements
            stim(cnt_el).c_mr = c_mr ; % MR coherence
            stim(cnt_el).s_m = S_m;  % M spectrum
            stim(cnt_el).s_r = S_r; % R spectrum
            stim(cnt_el).dmt_s_m = dmt_S_m;  % M dmt spectrum
            stim(cnt_el).dmt_s_r = dmt_S_r; % R dmt spectrum
            
            stim(cnt_el).c_mr_H = c_mr_H ; % MR coherence hits
            stim(cnt_el).s_m_H = S_m_H;  % M spectrum hits
            stim(cnt_el).s_r_H = S_r_H; % R spectrum hits
            stim(cnt_el).dmt_s_m_H = dmt_S_m_H;  % M dmt spectrum
            stim(cnt_el).dmt_s_r_H = dmt_S_r_H; % R dmt spectrum
            
            stim(cnt_el).c_mr_M = c_mr_M ; % MR coherence misses
            stim(cnt_el).s_m_M = S_m_M;  % M spectrum misses
            stim(cnt_el).s_r_M = S_r_M; % R spectrum misses
            stim(cnt_el).dmt_s_m_M = dmt_S_m_M;  % M dmt spectrum misses
            stim(cnt_el).dmt_s_r_M = dmt_S_r_M; % R dmt spectrum misses
            
            stim(cnt_el).n_trials_RS_STIM = size(trials,2);
            
            cnt_el = cnt_el + 1;  % -- count # of electrodes (modulators)
            
        end
        cnt_m = cnt_m + 1;
    end % --- end of all modulator channels
    
end

save(strcat(dir_Stim,sprintf('/coh_spec_mr_sameRStrails_fk_%d_W_%d.mat',fk,W)),'stim');

