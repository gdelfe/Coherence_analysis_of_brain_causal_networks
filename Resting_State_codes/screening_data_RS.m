
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code check the artifacts in the lfp 
%
%    @ Gino Del Ferraro, December 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_base = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Resting_state';
step = 110;

fid = fopen(strcat(dir_base,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

% -- load structure files
load(strcat(dir_base,'/session_AM.mat'))
load(strcat(dir_base,'/session_MA.mat'))


% -- print structures on stdout
%format short
for s=1:size(sess_info{1},1)
    session_AM(s)
    session_MA(s)
end


cnt_m = 1; % counter for the modulator-receiver/sender coherencies 
cnt_sr = 1; % counter sender-receiver coherencies 
cnt_el = 0; % counter for how many modulators excluding the receivers modulators 

for i=1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    % addpath('/vol/sas8/Maverick_RecStim_vSUBNETS220/160125/004')
    addpath(sprintf('/vol/sas8/Maverick_RecStim_vSUBNETS220/%s/%s/',sess_info{2}{i},sess_info{3}{i})) % add path of the specific RS session
    
    % file = 'rec004.Frontal.lfp.dat'
    file = sprintf('rec%s.Frontal.lfp.dat',sess_info{3}{i})
    fid = fopen(file);
    format = 'float=>single';
    
    CH = 220; % tot number of channels
    FS = 1000; % sampling
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    data = fread(fid,[CH,inf],format); % load the RS data
    % h = fread(fid,[CH,diff(bn)*FS./1e3],format);
    % ---- bipolar referencing, pairs of electrodes
    dir_Sess = strcat(dir_base,sprintf('/Sess_%d',Sess));
    if ~exist(dir_Sess, 'dir')
        mkdir(dir_Sess)
    end
    
    % -- load list electrodes, sender, receiver
    electrode = importdata(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairs
    receiver = importdata(strcat(dir_Sess,sprintf('/receiver_Sess_%d.txt',Sess)));  % ---- receiver pair
    sender = importdata(strcat(dir_Sess,sprintf('/sender_Sess_%d.txt',Sess))); % ---- sender pair
    
    receiver_idx = find(electrode == receiver(1));
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
    delta = 1000;
    cnt = 1;
    for j = 0:delta:(tot_time - delta)
        lfp_S(cnt,:) = lfp_S_ns(1,j+1:j+delta);
        lfp_R(cnt,:) = lfp_R_ns(1,j+1:j+delta);
        
        lfp_E(:,cnt,:) =  lfp_E_ns(:,j+1:j+delta,:);
        cnt = cnt + 1;
    end
    

    
    % -- plot sender / receiver LFP
    fig = figure;
    plot(lfp_S_ns);
    hold on
    plot(lfp_R_ns);
    grid on 
    title('Lfp sender and receiver ','FontSize',11);
    xlabel('time (sec)');
    ylabel('Lfp');
    legend('Sender','Receiver','FontSize',11);
    set(gcf, 'Position',  [100, 600, 1000, 600])

    fname = strcat(dir_Sess,'/lfp_Sender_Receiver.png');
    saveas(fig,fname)
    
    % -- plot modulators LFP
    mod_Ch = session_AM(i).mod_idx; % causal modulator channel
    figure;
    for Ch = mod_Ch
       plot(lfp_E_ns(Ch,:));
       hold on 
    end
    title('Lfp modulator(s)','FontSize',11);
    xlabel('time (sec)');
    ylabel('Lfp');
    grid on
    set(gcf, 'Position',  [100, 600, 1000, 600])

    fname = strcat(dir_Sess,'/lfp_Modulators.png');
    saveas(fig,fname)
    
    
    % --- Store session info and LFP
    sess(i).sess_idx = [i, Sess];
    % -- store LFP 
    sess(i).lfp_S = lfp_S;
    sess(i).lfp_R = lfp_R;
    sess(i).lfp_E = lfp_E; % -- all electrodes 
    
    sess(i).mod_idx = mod_Ch; % -- modulator idx
    sess(i).rec_idx = receiver_idx;
    
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
    
    
    % -- outliers Sender 
    th_S = 4*std(lfp_S_ns,[],2); % -- threshold for LFP Sender 
    max_S_split = max(abs(lfp_S),[],2);  % -- max of LFP for each time window   
    sess(i).outliers_S = find(max_S_split > th_S)';
    
    % -- outliers Receiver 
    th_R = 4*std(lfp_R_ns,[],2); % -- threshold for LFP Sender
    max_R_split = max(abs(lfp_R),[],2);  % -- max of LFP for each time window
    sess(i).outliers_R = find(max_R_split > th_R)';
    
    % %%%%%%%%%%% Sender and Receiver %%%%%%%%%%%%%%%%%%%
    lfp_S_clean = lfp_S;
    lfp_R_clean = lfp_R;
    
    % %%%%%%% ALL Electrodes %%%%%%%%%%%%%%%%%%%%%
    lfp_E_clean = lfp_E;

  
    % -- outliers all electrodes  
    th_E = 4*std(lfp_E_ns,[],2); % -- threshold for LFP all electrodes 
    max_E_split = max(abs(lfp_E),[],3);  % -- max of LFP for each time window, for each channel
    
    outliers = [];
    outliers = [outliers, sess(i).outliers_S]; % -- stuck up outliers sender
    outliers = [outliers, sess(i).outliers_R]; % -- stuck up outliers receiver 
    
    cnt_m = 1;
    for Ch = mod_Ch % -- for each modulator find outliers
        sess(i).outliers_E(cnt_m).idx = find(max_E_split(Ch,:) > th_E(Ch)); % -- find outliers for this channel 
        outliers = [outliers, sess(i).outliers_E(cnt_m).idx];    % -- stuck up outliers modulators
        cnt_m = cnt_m + 1;
    end 

        
    outliers = unique(outliers)  % -- remove repeated entries in outliers 
    
    % -- remove outliers from sender and receiver
    lfp_S_clean(outliers,:) = [];
    lfp_R_clean(outliers,:) = [];
    
    % --- save lfp onto strucure
    sess(i).lfp_S_clean = lfp_S_clean;
    sess(i).lfp_R_clean = lfp_R_clean;
    
    
    % -- remove outliers from modulators  
    cnt_m = 1;
    for Ch = mod_Ch % -- for all the electrodes
        lfp_E_clean = sq(lfp_E(Ch,:,:));          % -- get lfp for only that channel
        lfp_E_clean(outliers,:) = [];       % -- remove outliers for this channel
        sess(i).lfp_E_clean(cnt_m).lfp = lfp_E_clean;   % -- save to structure
        cnt_m = cnt_m + 1; 
    end
    
    current_session = sess(i);
    save(strcat(dir_Sess,'/session_split.mat'),'current_session');

    
end
    
save(strcat(dir_base,'/all_sessions_split.mat'),'sess','-v7.3');


for i=1:20

    display(['Session --- ',num2str(i)])
    badSess(i).std_R
    badSess(i).std_S
    badSess(i).std_E
end

  





