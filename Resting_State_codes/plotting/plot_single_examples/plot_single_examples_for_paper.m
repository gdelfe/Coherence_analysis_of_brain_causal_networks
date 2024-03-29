
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the MR and MS coherence across monkeys. It uses all the
% last_rec sessions for Maverick (because Maverick does not have bad
% sessions) and it uses a combination of rec001 and rec002 without bad
% sessions for Archie.
%
% The code computes MS and MR for a varying number of TOP modulators, and
% can be run also for all the modulators (N=100).
%
% It is flexible to both theta and beta frequency calculations
%
% OUTPUT: .png, .fig files and .mat files with the data for the coherence.
% All the output files are saved in the dir_both_monkeys path
%
% @Gino Del Ferraro, March 2021, NYU, Pesaran Lab


clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_fig = strcat(dir_main,'Figures_paper');

freq_band = 'theta_band';

% -- beta
% N_list = [10,20,30,40,50,100];
% N_mav_list = [4,8,12,16,20,41];
% N_arc_list = [5,10,15,20,25,51];

% -- theta
% N_list = [10, 20, 30,40,50,100];
% N_mav_list = [10, 20, 30,40,50,101];
% N_arc_list = [4,8,12,17,22,44];

%-- all of them
N_list = [100];
N_mav_list = [101];
N_arc_list = [44];



for i=1:length(N_list)
    
    close all
    N = N_list(i); % --- max number of modulators
    N_mav = N_mav_list(i);
    N_arc = N_arc_list(i);
    
    decode = 'AUC';
    decoding = 'AUC';
    titleN = sprintf('%d%% top modulators - last rec-rec001 002, %s',N,decode);
    
    % -- rec001-002
    arc_bad_sess = [8,22,30,31]; % -- theta band
    %         arc_bad_sess = [14,22,30,41]; % -- beta band
    
    % -- MOVIE
    %     arc_bad_sess = [8,22,30,31]; % -- theta band
    %         arc_bad_sess = [14,16,22,30,41]; % -- beta band
    
    recording_mav = 'last_recording';
    recording_arc = 'rec001_002_all_sessions';
    
    namef_mav = ''; % -- loading file for coherence averages
    namef_arc = '_rec001_002_all_sess';
    
    name = 'lat_rec-rec001_002_no_bad_sess';
    namefig = sprintf('%s.fig',name); % -- write out fig and png name
    namepng = sprintf('%s.png',name);
    
    recording_both = 'last_rec-rec001_002'; % -- output directory
    
    % -- decod accuracy
    %     fname_list_modulators_mav = '/modulators_sorted_decod_accuracy.txt';
    %     fname_list_modulators_arc = '/modulators_unsorted_decod_accuracy_all.txt';
    
    % -- AUC
    fname_list_modulators_mav = '/modulators_sorted_AUC.txt';
    fname_list_modulators_arc = '/modulators_unsorted_AUC_all.txt';
    
    
    dir_list_mav = strcat(dir_main,sprintf('Maverick/Resting_state/%s/Modulators_controls',freq_band));
    dir_list_arc = strcat(dir_main,sprintf('Archie/Resting_state/%s/Modulators_controls',freq_band));
    
    dir_Maverick = strcat(dir_main,sprintf('Maverick/Resting_State/%s/Modulators_Controls_avg_results/%s',freq_band));
    dir_Archie = strcat(dir_main,sprintf('Archie/Resting_State/%s/Modulators_Controls_avg_results/%s',freq_band));
    
    dir_Maverick_avg = strcat(dir_main,sprintf('Maverick/Resting_State/%s/Modulators_Controls_avg_results/%s',freq_band,recording_mav));
    dir_Archie_avg = strcat(dir_main,sprintf('Archie/Resting_State/%s/Modulators_Controls_avg_results/%s',freq_band,recording_arc));
    
    fk = 200; W = 5;
    dir_both_monkeys = strcat(dir_main,sprintf('both_monkeys/%s/modulators_vs_controls/%s/%s',freq_band,recording_both,decoding));
    if ~exist(dir_both_monkeys, 'dir')
        mkdir(dir_both_monkeys)
    end
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       MAVERICK              %
    
    mod_list = importdata(strcat(dir_list_mav,sprintf('%s',fname_list_modulators_mav)));
    fid = fopen(strcat(dir_Maverick,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
    sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
    fclose(fid);
    display(['Tot Maverick modulators is N = ',num2str(size(mod_list,1))])
    
    
    % -- select the first N index
    mod_idx = mod_list(1:N_mav,4);
    sess_numb_M = unique(mod_list(1:N_mav,1));
    
    % -- find the session index corresponding to the session with top modulators
    sess_idx_M = [];
    for i=1:length(sess_numb_M)
        sess_idx_M = [sess_idx_M, find(sess_info{1}==sess_numb_M(i))];
    end
    
    
    % %%%%%%%%% MODULATORS Maverick %%%%%%
    load(strcat(dir_Maverick_avg,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,namef_mav))); % structure mod
    load(strcat(dir_Maverick_avg,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,namef_mav))); % structure stim
    struct_mod_mav = mod;
    struct_stim_mav = stim;
    
    % -- get coherences for the specific modulators and specific sessions
    struct_mod_mav = struct_mod_mav(mod_idx);
    struct_stim_mav = struct_stim_mav(sess_idx_M);
    
    mod_list_mav = mod_list;
    
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       ARCHIE                %
    mod_list = importdata(strcat(dir_list_arc,sprintf('%s',fname_list_modulators_arc)));
    fid = fopen(strcat(dir_Archie,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
    sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
    fclose(fid);
    
    % -- remove modulators belonging to bad sessions
    idx_bad = [];
    for bad_s = arc_bad_sess
        idx_bad = [idx_bad; mod_list(mod_list(:,1)== bad_s,4)]; % -- find indexes of removed modulators
        mod_list(mod_list(:,1)==bad_s,:)=[];
    end
    
    mod_list = [mod_list(:,1:3), double(1:size(mod_list,1))']; % session, modulator idx, decod accuracy, order in coherence avg structure
    mod_list = sortrows(mod_list,3,'descend'); % sort in descending order
    display(['Tot Archie modulators after removal bad sessions is N = ',num2str(size(mod_list,1))])
    
    % -- select the first N index
    mod_idx = mod_list(1:N_arc,4);
    sess_numb_A = unique(mod_list(1:N_arc,1));
    
    % -- find the session index corresponding to the session with top modulators
    sess_idx_A = [];
    for i=1:length(sess_numb_A)
        sess_idx_A = [sess_idx_A, find(sess_info{1}==sess_numb_A(i))];
    end
    
    
    % %%%%%%%%% MODULATORS Archie %%%%%%
    load(strcat(dir_Archie_avg,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,namef_arc))); % structure mod
    load(strcat(dir_Archie_avg,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,namef_arc))); % structure stim
    struct_mod_arc = mod;
    struct_stim_arc = stim;
    
    % -- remove coherence channels associated to the bad modulators
    for idx = idx_bad
        struct_mod_arc(idx) = [];
    end
    
    struct_mod_arc = struct_mod_arc(mod_idx);
    struct_stim_arc = struct_stim_arc(sess_idx_A);
    
    
    % %%%%%%%%% Computes modulators for Maverick and Archie %%%%%%
    modulators = mean_coh_and_spec_across_monkeys(struct_mod_mav,struct_stim_mav,struct_mod_arc,struct_stim_arc);
    
    
    
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONTROLS SAME AREA MAVERICK              %
    
    load(strcat(dir_Maverick_avg,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s',fk,W,namef_mav)));
    load(strcat(dir_Maverick_avg,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s',fk,W,namef_mav)));
    mod_ctrl_SA_mav = mod;
    stim_ctrl_SA_mav = stim;
    
    ctrl_list = importdata(strcat(dir_list_mav,'/control_list_same_area.txt'));
    
    % -- get the index of the controls associated with sessions where the top modulators are found
    ctrl_idx_M = [];
    for i = sess_numb_M'
        ctrl_idx_M = [ctrl_idx_M; find(ctrl_list(:,1)== i)];
    end
    
    mod_ctrl_SA_mav = mod_ctrl_SA_mav(ctrl_idx_M);
    stim_ctrl_SA_mav = stim_ctrl_SA_mav(sess_idx_M);
    
    ctrl_SA_mav = ctrl_list;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONTROLS SAME AREA ARCHIE              %
    
    load(strcat(dir_Archie_avg,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s',fk,W,namef_arc)));
    load(strcat(dir_Archie_avg,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s',fk,W,namef_arc)));
    mod_ctrl_SA_arc = mod;
    stim_ctrl_SA_arc = stim;
    
    ctrl_list = importdata(strcat(dir_list_arc,'/control_list_same_area.txt'));
    
    % -- get the index of the controls associated with sessions where the top modulators are found
    ctrl_idx_A = [];
    for i = sess_numb_A'
        ctrl_idx_A = [ctrl_idx_A; find(ctrl_list(:,1)== i)];
    end
    
    
    mod_ctrl_SA_arc = mod_ctrl_SA_arc(ctrl_idx_A);
    stim_ctrl_SA_arc = stim_ctrl_SA_arc(sess_idx_A);
    
    
    ctrl_SA = mean_coh_and_spec_across_monkeys(mod_ctrl_SA_mav,stim_ctrl_SA_mav,mod_ctrl_SA_arc,stim_ctrl_SA_arc);
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONTROLS OTHER AREAS MAVERICK              %
    
    load(strcat(dir_Maverick_avg,sprintf('/coh_spec_m_Controls_other_areas_fk_%d_W_%d%s',fk,W,namef_mav)));
    load(strcat(dir_Maverick_avg,sprintf('/coh_spec_sr_Controls_other_areas_fk_%d_W_%d%s',fk,W,namef_mav)));
    mod_ctrl_OA_mav = mod;
    stim_ctrl_OA_mav = stim;
    
    ctrl_list = importdata(strcat(dir_list_mav,'/control_list_other_areas.txt'));
    
    % -- get the index of the controls associated with sessions where the top modulators are found
    ctrl_idx_M = [];
    for i = sess_numb_M'
        ctrl_idx_M = [ctrl_idx_M; find(ctrl_list(:,1)== i)];
    end
    
    mod_ctrl_OA_mav = mod_ctrl_OA_mav(ctrl_idx_M);
    stim_ctrl_OA_mav = stim_ctrl_OA_mav(sess_idx_M);
    
    ctrl_OA_mav = ctrl_list;

    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONTROLS OTHER AREAS ARCHIE              %
    
    load(strcat(dir_Archie_avg,sprintf('/coh_spec_m_Controls_other_areas_fk_%d_W_%d%s',fk,W,namef_arc)));
    load(strcat(dir_Archie_avg,sprintf('/coh_spec_sr_Controls_other_areas_fk_%d_W_%d%s',fk,W,namef_arc)));
    mod_ctrl_OA_arc = mod;
    stim_ctrl_OA_arc = stim;
    
    ctrl_list = importdata(strcat(dir_list_arc,'/control_list_other_areas.txt'));
    
    % -- get the index of the controls associated with sessions where the top modulators are found
    ctrl_idx_A = [];
    for i = sess_numb_A'
        ctrl_idx_A = [ctrl_idx_A; find(ctrl_list(:,1)== i)];
    end
    
    
    mod_ctrl_OA_arc = mod_ctrl_OA_arc(ctrl_idx_A);
    stim_ctrl_OA_arc = stim_ctrl_OA_arc(sess_idx_A);
    
    
    ctrl_OA = mean_coh_and_spec_across_monkeys(mod_ctrl_OA_mav,stim_ctrl_OA_mav,mod_ctrl_OA_arc,stim_ctrl_OA_arc);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f = linspace(1,fk,size(modulators.mean_coh_ms,2)); % frequency values (range)
    
    
    dir_RS = strcat(dir_main,sprintf('Maverick/Resting_state/theta_band'));
    
    fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
    sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
    fclose(fid);


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Figure: Coherence MR vs CR -- single electrodes
    
    id_fig = [1 2 3];
    id = [6 10 7]; % index modulator

    id = 13;
    
    id_list = [3 6 7 9 10 23]
    % chosen 6, 7, 10
    
    idr = [6 20 29]; % index modulator
    ids = [6 10 13];
    
    
    for id = 20    ytop = [0.88 0.35 0.8];
    
    %         close all
    
    Sess = mod_list_mav(id,1);
    %     display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    dir_ctrl_SA = strcat(dir_RS,sprintf('/Sess_%d/Controls_same_area',Sess));
    dir_ctrl_OA = strcat(dir_RS,sprintf('/Sess_%d/Controls_other_areas',Sess));
    
    
    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    sess_data
    
    mod_ch = mod_list_mav(id,2)
    
    load(strcat(dir_ctrl_SA,'/session_controls_same_area_info.mat')); % --- dataG: all data info and LFP
    ctrl_SA = sess_All_controls_same_area
    
    ctrl_SA_mav(ctrl_SA_mav(:,1) == Sess,:)
    
    ctrl_SA.ctrl_idx(find(strcmp(ctrl_SA.RecordPairMRIlabels(ctrl_SA.ctrl_idx,1),'CN')))
    temp = ctrl_SA_mav(ctrl_SA_mav(:,2,:) == 9,:) %%%%%%%%%%%%%%%%%
    SA_idx = temp(temp(:,1) == Sess,3);
    
    load(strcat(dir_ctrl_OA,'/session_controls_other_areas_info.mat')); % --- dataG: all data info and LFP
    ctrl_OA = sess_All_controls_other_areas
    
    ctrl_OA_mav(ctrl_SA_mav(:,1) == Sess,:)
     
    ctrl_OA.RecordPairMRIlabels(ctrl_OA.ctrl_idx,1)'
    ctrl_OA.ctrl_idx
    idx=28; %%%%%%%%%%%%%%%%%%%%%%%
    
    temp = ctrl_OA_mav(ctrl_OA_mav(:,2,:) == idx,:);
    OA_idx = temp(temp(:,1) == Sess,3)
    
    
    figMS = figure;
    
    plot(f,abs(struct_mod_mav(id).c_ms),'color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);  hold on
    plot(f,abs(mod_ctrl_SA_mav(SA_idx).c_ms),'color',[255, 51, 153]/255,'LineWidth',1.5); hold on
    plot(f,abs(mod_ctrl_OA_mav(OA_idx).c_ms),'color',[255, 128, 128]/255,'LineWidth',1.5);
    
    grid on
    set(gca,'FontSize',14)
    xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
    ylabel('Coherence','FontName','Arial','FontSize',15);
%     title(sprintf('ID = %d',id));
    set(gcf, 'Position',  [100, 600, 600, 400])
    xticks([ 10 20 30 40 50 60 70 80 90])
    xlim([1 95])
    ylim([0 0.46])
    
    
    cnt_fig = 13;
    fname = strcat(dir_fig,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d_id_%d.pdf',N,cnt_fig));
    saveas(figMS,fname)
    clear fname
    
    
    
    figMR = figure;
    
    plot(f,abs(struct_mod_mav(id).c_mr),'color',[28 199 139]/255,'LineWidth',1.5);  hold on
    plot(f,abs(mod_ctrl_SA_mav(SA_idx).c_mr),'color',[50 250 93]/255,'LineWidth',1.5); hold on
    plot(f,abs(mod_ctrl_OA_mav(OA_idx).c_mr),'color',[19 148 92]/255,'LineWidth',1.5);
    
    
    grid on
    set(gca,'FontSize',14)
    xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
    ylabel('Coherence','FontName','Arial','FontSize',15);
%     title(sprintf('ID = %d',id));
    set(gcf, 'Position',  [100, 600, 600, 400])
    xticks([ 10 20 30 40 50 60 70 80 90])
    xlim([1 95])
    ylim([0 0.65])
    
    
    cnt_fig = 13;
    fname = strcat(dir_fig,sprintf('/coherency_MR_mod_vs_ctrl_both_monkeys_N_%d_id_%d.pdf',N,cnt_fig));
    saveas(figMR,fname)
    clear fname
    
    
    
    
    end 
    
   
    
    fname = strcat(dir_fig,sprintf('/coherency_MR_mod_vs_ctrl_both_monkeys_N_%d_id_%d.pdf',N,cnt_fig));
    saveas(figMR,fname)
    clear fname
       
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Figure: Coherence MS vs CS -- single electrodes
    
    id_fig = [1 2 3];
 
    
    id_r = [1 3 10]; % index contrls
    id_c = [1 3 10]; % index contrls
    ytop = [0.65 0.46 0.46];
    
    close all 
    for i=1:3
%         close all
        fig = figure;

        plot(f,abs(struct_mod_mav(idr(i)).c_ms),'color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);  hold on
        plot(f,abs(mod_ctrl_SA_mav(id_c(i)).c_ms),'color',[255, 51, 153]/255,'LineWidth',1.5); hold on
        plot(f,abs(mod_ctrl_OA_mav(id_c(i)).c_ms),'color',[255, 128, 128]/255,'LineWidth',1.5);     
        grid on
        
        set(gca,'FontSize',14)
        xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
        if i == 1
            ylabel('Coherence','FontName','Arial','FontSize',15);
        end
        
        if i == 3
            legend('Modulator - Sender','Control-SA - Sender','Control-OA - Sender','FontSize',10,'FontName','Arial')
        end
        set(gcf, 'Position',  [100, 600, 600, 400])
        xticks([ 10 20 30 40 50 60 70 80 90])
        xlim([1 95])
        ylim([0 ytop(i)])
        
            
        fname = strcat(dir_fig,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d_id_%d.pdf',N,id_fig(i)));
        saveas(fig,fname)
    
    end
    

    
end
% keyboard





