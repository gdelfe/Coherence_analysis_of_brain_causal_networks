

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code loads the entire RS time series for the CONTROLS and check for artifacts in the lfp 
% It splits the RS time series into 1 sec windows and remove all the windows
% which contains artifacts (th = 4*std(lfp))
%
% Controls are for the Sender but they are still stored in the
% lfp_E variable
%
% INPUT : sess_control_info.mat
%         structure containing all the info about the session, i.e. idx modulators, sender, etc...
%
% OUTPUT: 1. sess_control_lfp.mat: 
%               Strucure of data for each session containing:
%               a. lfp before removing artifacts
%               b. lfp after artifacts were removed 
%
%    @ Gino Del Ferraro, December 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
step = 110;

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)
name_structure_data_info = '/session_Sender_controls_other_areas_info.mat';

for i=1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    addpath(sprintf('/vol/sas8/Maverick_RecStim_vSUBNETS220/%s/%s/',sess_info{2}{i},sess_info{3}{i})) % add path of the specific RS session

   % file = 'rec004.Frontal.lfp.dat'
    file = sprintf('rec%s.Frontal.lfp.dat',sess_info{3}{i})
    fid = fopen(file);
    format = 'float=>single';
    
    CH = 220; % tot number of channels
    FS = 1000; % sampling
    
    data = fread(fid,[CH,inf],format); % load the RS data
    % h = fread(fid,[CH,diff(bn)*FS./1e3],format);
    % ---- bipolar referencing, pairs of electrodes
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sender = strcat(dir_RS,sprintf('/Sess_%d/Sender_controls_other_areas',Sess));
    
    

    %load(strcat(dir_Sess,'/session_control_info.mat')); % --- dataG: all data info and LFP
    load(strcat(dir_Sender,name_structure_data_info)); % --- dataG: all data info and LFP
    sess_control = sess_Send_ctrl_other_areas;
    clear sess_Send_ctrl_other_areas;

    % -- load list electrodes, sender, receiver
    electrode = sess_control.RecordPair; % ---- all electrode pairs
    receiver = sess_control.receiver_pair;  % ---- receiver pair
    sender = sess_control.sender_pair; % ---- sender pair

    % ---  time parameter
    tot_time = 150001;
    % ---  freq parameter for the masking
    fmin = 10;
    fmax = 40;
    
    % ---- Lfp of the resting state for that specific pair of electrodes
    lfp_E_ns = data(electrode(:,1),:) - data(electrode(:,2),:); % all the electrodes
    lfp_S_ns = data(sender(1),:) - data(sender(2),:); % sender
    lfp_R_ns = data(receiver(1),:) - data(receiver(2),:); % receiver
    
    % include signal up to time where signal is not corrupted
    lfp_E_ns = lfp_E_ns(:,1:tot_time);
    lfp_S_ns = lfp_S_ns(:,1:tot_time);
    lfp_R_ns = lfp_R_ns(:,1:tot_time);
    
    
    % create matrices to store the split data: trial x time
    lfp_S = zeros(floor(tot_time/FS),1000);
    lfp_R = zeros(floor(tot_time/FS),1000);
    lfp_E = zeros(size(lfp_E_ns,1),floor(tot_time/FS),1000); % channel x trial x time 
    
    
% %   %-- sanity check LFP
%     figure;
%     plot(lfp_S_ns)
%     hold on
%     plot(lfp_R_ns)
%     hold on 
%     plot(lfp_E_ns(6,:,:))
%     legend('sender','receiver','electrode')


%     --- Split the Lenghty RS time series into 1000 ms windows
%     format: channel x win_indx xtime. For R and S size_channel = 1  

    % --- change this operation with reshape 
    delta = 1000;
    cnt = 1;
    for j = 0:delta:(tot_time - delta)
        lfp_S(cnt,:) = lfp_S_ns(1,j+1:j+delta);
        lfp_R(cnt,:) = lfp_R_ns(1,j+1:j+delta);
        
        lfp_E(:,cnt,:) =  lfp_E_ns(:,j+1:j+delta,:);
        cnt = cnt + 1;
    end
    
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %  FIGURES          %%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % -- plot sender / receiver LFP
%     fig = figure;
%     plot(lfp_S_ns);
%     hold on
%     plot(lfp_R_ns);
%     grid on 
%     title('Lfp sender and receiver ','FontSize',11);
%     xlabel('time (sec)');
%     ylabel('Lfp');
%     legend('Sender','Receiver','FontSize',11);
%     set(gcf, 'Position',  [100, 600, 1000, 600])
% 
%     fname = strcat(dir_Sess,'/lfp_Sender_Receiver.png');
%     saveas(fig,fname)
%     
    % -- plot controls LFP
%     
%     ctrl_Ch = sess_control.ctrl_idx; % control modulators 
%     figure;
%     ctrl_Ch = 51;
%     for Ch = ctrl_Ch
%        plot(lfp_E_ns(Ch,:));
%        hold on 
%     end
%     title('Lfp control(s) all controls same area','FontSize',11);
%     xlabel('time (sec)');
%     ylabel('Lfp');
%     grid on
%     set(gcf, 'Position',  [100, 600, 1000, 600])
%     
% 
%     fname = strcat(dir_Sess,'/lfp_all_Controls_same_area.png');
%     saveas(fig,fname)
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ASSIGN LFP SPLIT ALL ELECTRODES with artifacts included
    % ------------------------------------------------------------
    
    sess_control_lfp = sess_control; 
    
    % -- store LFP split 
    sess_control_lfp.lfp_S = lfp_S;
    sess_control_lfp.lfp_R = lfp_R;
    sess_control_lfp.lfp_E = lfp_E; % -- all electrodes 
    

    % %%%%%%%%% SESSIONS WITH HIGH STD %%%%%%%%%%%%%%%%%%%%%%%%
    % --- find sessions with high std (>150)
    if std(lfp_S_ns,[],2) > 150  badSess(i).std_S = std(lfp_S_ns,[],2);
        display(['Bad Sess -- ',num2str(Sess),' Sender std > 150 -- ']); end 
    
    if std(lfp_R_ns,[],2) > 150  badSess(i).std_R = std(lfp_R_ns,[],2); 
        display(['Sess -- ',num2str(Sess),' Receiver -- ',num2str(Ch),' std > 150']); end 

    ctrl_Ch = sess_control.ctrl_idx;
    cnt_M = 1;
    for Ch = ctrl_Ch
            if std(lfp_E_ns(Ch,:,:),[],2) > 150  badSess(i).std_E(cnt_M).std = std(lfp_E_ns(Ch,:,:),[],2); 
            display(['Sess -- ',num2str(Sess),' Electrode -- ',num2str(Ch),' std > 150'])
            end 
            cnt_M = cnt_M + 1;
    end 
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  FIND OUTLIERS, i.e. windows with artifacts, and remove them  %%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % -- outliers Sender 
    th_S = 4*std(lfp_S_ns,[],2); % -- threshold for LFP Sender 
    max_S_split = max(abs(lfp_S),[],2);  % -- max of LFP for each time window   
    sess_control_lfp.outliers_S = find(max_S_split > th_S)';
    
    % -- outliers Receiver 
    th_R = 4*std(lfp_R_ns,[],2); % -- threshold for LFP Sender
    max_R_split = max(abs(lfp_R),[],2);  % -- max of LFP for each time window
    sess_control_lfp.outliers_R = find(max_R_split > th_R)';

  
    % -- outliers all electrodes  
    th_E = 4*std(lfp_E_ns,[],2); % -- threshold for LFP all electrodes 
    max_E_split = max(abs(lfp_E),[],3);  % -- max of LFP for each time window, for each channel
    
    outliers = [];
    outliers = [outliers, sess_control_lfp.outliers_S]; % -- stuck up outliers sender
    outliers = [outliers, sess_control_lfp.outliers_R]; % -- stuck up outliers receiver 
    
    % -- modulators outliers
    cnt_m = 1;
    for Ch = ctrl_Ch % -- for each modulator find outliers
        sess_control_lfp.outliers_E(cnt_m).idx = find(max_E_split(Ch,:) > th_E(Ch)); % -- find outliers for this channel 
        sess_control_lfp.outliers_tot(cnt_m).idx = [outliers, sess_control_lfp.outliers_E(cnt_m).idx];    % -- stuck up outliers of receiver, sender, and modulators
        sess_control_lfp.outliers_tot(cnt_m).idx = unique(sess_control_lfp.outliers_tot(cnt_m).idx); % -- remove repeated entries in outliers M, S, R
        cnt_m = cnt_m + 1;
    end
    
    save(strcat(dir_Sender,'/sess_Send_ctrl_other_areas_lfp.mat'),'sess_control_lfp');

    
end


keyboard 

% save(strcat(dir_RS,'/all_sessions_split.mat'),'sess','-v7.3');


for i=1:20

    display(['Session --- ',num2str(i)])
    badSess(i).std_R
    badSess(i).std_S
    badSess(i).std_E
end

  




