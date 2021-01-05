
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

name_struct_input = '/session_all_controls_other_areas_info.mat'

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
    
    load(strcat(dir_RS_Sess,name_struct_input)); % --- dataG: all data info and LFP
    sess_control = sess_All_controls_other_areas;
    clear sess_All_controls_other_areas;
    
    % -- load list electrodes, sender, receiver
    electrode = sess_control.RecordPair; % ---- all electrode pairs
    receiver = sess_control.receiver_pair;  % ---- receiver pair
    sender = sess_control.sender_pair; % ---- sender pair
    
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
    
    
    mod_Ch = sess_control.mod_idx; % -- modulators (not controls!) index
    ctrl_Ch = sess_control.ctrl_idx;
    
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- control Ch:  ',num2str(ctrl_Ch)])
    
    dir_Modulators = strcat(dir_Sess,'/Controls_other_areas');
    if ~exist(dir_Modulators, 'dir')
        mkdir(dir_Modulators)
    end
    
    
    indx_list = [];
    
    for Ch = ctrl_Ch % for all the modulators in the session
        
        close all
        
        
        if Ch ~= sess_control_lfp.receiver_idx % if the electrode is not the receiver itself
            
            
            indx_list = [indx_list, cnt_el]; % store the cnt number --- needed for multiple plotting
            
            % -- labels Hits and Misses
            hitIndx = Data.spec.lfp.DetectedIndx{Ch}; % labels for the hits (which trial was a hit)
            missIndx = Data.spec.lfp.notDetectedIndx{Ch}; % labels for the misses (which trial was a miss)
            
            % -- coherence of modulator-receiver
            display(['Computing modulator-receiver coherence...'])
            tic
            [c_mr,f,S_m,S_r] = coherency(sq(lfp_E(:,Ch,:)),lfp_R,[N W],fs,fk,pad,0.05,1,11);
            toc
            % -- coherence of modulator-receiver HITS
            [c_mr_H,f_H,S_m_H,S_r_H] = coherency(sq(lfp_E(hitIndx,Ch,:)),lfp_R(hitIndx,:),[N W],fs,fk,pad,0.05,1,11);
            % -- coherence of modulator-receiver MISSES
            [c_mr_M,f_M,S_m_M,S_r_M] = coherency(sq(lfp_E(missIndx,Ch,:)),lfp_R(missIndx,:),[N W],fs,fk,pad,0.05,1,11);
            
            
            % --- FIGURE --------- %%
            % -- Coherence vs frequency --- %
            fig = figure;
            plot(f,abs(c_mr),'color',[0, 15, 26]/255)
            hold on
            plot(f,abs(c_mr_H),'color',[0, 153, 255]/255)
            hold on
            plot(f,abs(c_mr_M),'color',[255, 153, 51]/255)
            grid on
            title(sprintf('STIM: Abs coherence vs frequency, ch = %d, controls other areas',Ch),'FontSize',10);
            legend('M-R coherence','M-R coherence HITS','M-R coherence MISSES')
            %         xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Modulators,sprintf('/coherency_ch_%d_N_%.2f_W_%d.jpg',Ch,N,W));
            saveas(fig,fname);
            
            
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
            title(sprintf('STIM: Abs coherence vs frequency, ch = %d, controls other areas',Ch),'FontSize',10);
            legend('M Spectrum Hits','R Spectrum Hits','M Spectrum Misses','R Spectrum Misses')
            xlim([0 60])
            set(gcf, 'Position',  [100, 600, 1000, 500])
            
            fname = strcat(dir_Modulators,sprintf('/spectrum_ch_%d_N_%.2f_W_%d.jpg',Ch,N,W));
            saveas(fig,fname);
            
            
            % -- structure assignements
            stim(cnt_el).c_mr = c_mr ; % MR coherence
            stim(cnt_el).s_m = S_m;  % M spectrum
            stim(cnt_el).s_r = S_r; % R spectrum
            
            stim(cnt_el).c_mr_H = c_mr_H ; % MR coherence hits
            stim(cnt_el).s_m_H = S_m_H;  % M spectrum hits
            stim(cnt_el).s_r_H = S_r_H; % R spectrum hits
            
            stim(cnt_el).c_mr_M = c_mr_M ; % MR coherence misses
            stim(cnt_el).s_m_M = S_m_M;  % M spectrum misses
            stim(cnt_el).s_r_M = S_r_M; % R spectrum misses
            
            cnt_el = cnt_el + 1;  % -- count # of electrodes (modulators)
        end
    end % --- end of all modulator channels
    
end


save(strcat(dir_Stim,sprintf('/coh_spec_mr_controls_other_areas_fk_%d_W_%d.mat',fk,W)),'stim');

