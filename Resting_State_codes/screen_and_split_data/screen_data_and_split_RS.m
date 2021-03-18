
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code loads the entire RS time series for the MODULATORS and check for artifacts in the lfp 
% It splits the RS time series into 1 sec windows and remove all the windows
% which contains artifacts (th = 4*std(lfp))
%
% INPUT : sess_data_info.mat
%         structure containing all the info about the session, i.e. idx modulators, sender, etc...
%
% OUTPUT: 1. sess_data_lfp.mat: 
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
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'beta_band';
monkey = 'Archie'


dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Stim = strcat(dir_main,sprintf('%s/Stim_data/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

for i=1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
%     addpath(sprintf('/vol/sas8/Maverick_RecStim_vSUBNETS220/%s/%s/',sess_info{2}{i},sess_info{3}{i})) % add path of the specific RS session
    addpath(sprintf('/vol/sas5a/Archie_RecStim_vSUBNETS220_2nd/%s/%s/',sess_info{2}{i},sess_info{3}{i})) % add path of the specific RS session

    
   % file = 'rec004.Frontal.lfp.dat'
%     file = sprintf('rec%s.Frontal.lfp.dat',sess_info{3}{i}) % -- Maverick
    file = sprintf('rec%s.Frontal_1.lfp.dat',sess_info{3}{i}) % -- Archie
    fid = fopen(file);
    format = 'float=>single';
    
    CH = 220; % tot number of channels
    FS = 1000; % sampling
    
    data = fread(fid,[CH,inf],format); % load the RS data
    % h = fread(fid,[CH,diff(bn)*FS./1e3],format);
    % ---- bipolar referencing, pairs of electrodes
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));
    if ~exist(dir_Sess, 'dir')
        mkdir(dir_Sess)
    end
    

    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP

    % -- load list electrodes, sender, receiver
    electrode = sess_data.RecordPair; % ---- all electrode pairs
    receiver = sess_data.receiver_pair;  % ---- receiver pair
    sender = sess_data.sender_pair; % ---- sender pair

    % ---  time parameter
    tot_time = 150001;
    
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
    
    mod_Ch = sess_data.mod_idx; % control modulators 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  FIGURES          %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
%     % -- plot controls LFP
%     
%     figure;
%     for Ch = mod_Ch
%        plot(lfp_E_ns(Ch,:));
%        hold on 
%     end
%     title('Lfp control(s)','FontSize',11);
%     xlabel('time (sec)');
%     ylabel('Lfp');
%     grid on
%     set(gcf, 'Position',  [100, 600, 1000, 600])
% 
%     fname = strcat(dir_Sess,'/lfp_Controls.png');
%     saveas(fig,fname)
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ASSIGN LFP SPLIT ALL ELECTRODES with artifacts included
    % ------------------------------------------------------------
    
    sess_data_lfp = sess_data; 
    
    % -- store LFP split 
    sess_data_lfp.lfp_S = lfp_S;
    sess_data_lfp.lfp_R = lfp_R;
    sess_data_lfp.lfp_E = lfp_E; % -- all electrodes 
    

    % %%%%%%%%% SESSIONS WITH HIGH STD %%%%%%%%%%%%%%%%%%%%%%%%
    % --- find sessions with high std (>150)
    if std(lfp_S_ns,[],2) > 150  badSess(i).std_S = std(lfp_S_ns,[],2);
        display(['Sess -- ',num2str(Sess),' Sender -- ']); end 
    
    if std(lfp_R_ns,[],2) > 150  badSess(i).std_R = std(lfp_R_ns,[],2); 
        display(['Sess -- ',num2str(Sess),' Receiver -- ',num2str(Ch)]); end 

    cnt_M = 1;
    for Ch = mod_Ch
            if std(lfp_E_ns(Ch,:,:),[],2) > 150  badSess(i).std_E(cnt_M).std = std(lfp_E_ns(Ch,:,:),[],2); 
            display(['Sess -- ',num2str(Sess),' Modulator -- ',num2str(Ch)])
            end 
            cnt_M = cnt_M + 1;
    end 
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  FIND OUTLIERS, i.e. windows with artifacts, and remove them  %%%%%%%
    
    % -- outliers Sender 
    th_S = 4*std(lfp_S_ns,[],2); % -- threshold for LFP Sender 
    max_S_split = max(abs(lfp_S),[],2);  % -- max of LFP for each time window   
    sess_data_lfp.outliers_S = find(max_S_split > th_S)';
    
    % -- outliers Receiver 
    th_R = 4*std(lfp_R_ns,[],2); % -- threshold for LFP Sender
    max_R_split = max(abs(lfp_R),[],2);  % -- max of LFP for each time window
    sess_data_lfp.outliers_R = find(max_R_split > th_R)';

  
    % -- outliers all electrodes  
    th_E = 4*std(lfp_E_ns,[],2); % -- threshold for LFP all electrodes 
    max_E_split = max(abs(lfp_E),[],3);  % -- max of LFP for each time window, for each channel
    
    outliers = [];
    outliers = [outliers, sess_data_lfp.outliers_S]; % -- stuck up outliers sender
    outliers = [outliers, sess_data_lfp.outliers_R]; % -- stuck up outliers receiver 
    
    % -- control outliers
    cnt_m = 1;
    for Ch = mod_Ch % -- for each modulator find outliers
        sess_data_lfp.outliers_E(cnt_m).idx = find(max_E_split(Ch,:) > th_E(Ch)); % -- find outliers for this channel
        sess_data_lfp.outliers_tot(cnt_m).idx = [outliers, sess_data_lfp.outliers_E(cnt_m).idx];    % -- stuck up outliers of receiver, sender, and this control electrode
        sess_data_lfp.outliers_tot(cnt_m).idx = unique(sess_data_lfp.outliers_tot(cnt_m).idx); % -- remove repeated entries in outliers C, S, R
        cnt_m = cnt_m + 1;
    end

    sess_data_lfp
    save(strcat(dir_Sess,'/session_data_lfp.mat'),'sess_data_lfp');

    
end


keyboard 

% save(strcat(dir_RS,'/all_sessions_split.mat'),'sess','-v7.3');


for i=1:2
% 
% %     display(['Session --- ',num2str(i)])
% %     badSess(i).std_R
% %     badSess(i).std_S
% %     badSess(i).std_E
% %     sess_data_lfp(i)
%         
    Sess = sess_info{1}(i); % Session number
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d',Sess));
    load(strcat(dir_Sess,'/sess_data_lfp.mat')); % --- dataG: all data info and LFP
    sess_data_lfp
    
end
% 
%   





