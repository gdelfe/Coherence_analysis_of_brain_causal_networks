
clear all; close all;

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%
addpath('/mnt/pesaranlab/People/Gino/DL-modulators/Gino_codes')
dir_RS = '/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State';

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% -- load structure files
newAM = load(strcat(dir_RS,'/session_AM.mat'))
session_AM = newAM.session_AM;

newMA = load(strcat(dir_RS,'/session_MA.mat'))
session_MA = newMA.session_MA;


for i=1:size(sess_info{1},1) % for all the sessions with modulator
    
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
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',   out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    data = fread(fid,[CH,inf],format); % load the RS data
    % h = fread(fid,[CH,diff(bn)*FS./1e3],format);
    % ---- bipolar referencing, pairs of electrodes
    dir_Sess = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d',Sess);
    if ~exist(dir_Sess, 'dir')
        mkdir(dir_Sess)
    end
    
    % -- load list electrodes, sender, receiver
    electrode = importdata(strcat(dir_Sess,sprintf('/recorded_pairs_modulators_Sess_%d.txt',Sess))); %Data.RecordPair;   % ---- all potential modulator pairs
    receiver = importdata(strcat(dir_Sess,sprintf('/receiver_Sess_%d.txt',Sess)));  % ---- receiver pair
    sender = importdata(strcat(dir_Sess,sprintf('/sender_Sess_%d.txt',Sess))); % ---- sender pair
    
    
    % ---  time parameter
    tot_time = 150000;
    % ---  freq parameter for the masking
    fmin = 10;
    fmax = 40;
    
    % ---- Lfp of the resting state for that specific pair of electrodes
    lfpRS = data(electrode(:,1),:) - data(electrode(:,2),:); % all potential modulators
    lfp_S = data(sender(1),:) - data(sender(2),:); % sender
    lfp_R = data(receiver(1),:) - data(receiver(2),:); % receiver
    
    % include signal up to time where signal is not corrupted
    lfpRS = lfpRS(:,1:tot_time);
    lfp_S = lfp_S(:,1:tot_time);
    lfp_R = lfp_R(:,1:tot_time);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Receiver Spectrogram
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % parameters used by Shaoyu
    k = 4;
    fk = [0 60];
    tapers = [0.5 5];
    dn = 0.005;
    fs = FS;
    pad = 2;
    
    % --- receiver spectrogram
    [spec_R, fR , tiR] = tfspec_GINO(lfp_R(:,:),tapers,fs,dn,fk,pad,0.05,0,1);
    fig = figure; tvimage(sq(log(spec_R(1:500,:)))); title(sprintf('Receiver - N = %.2f, W = %d, dn = %.4f, k = %d',tapers(1),tapers(2),dn,k)); colorbar;
    
    fname = strcat(dir_Sess,sprintf('/spectrogram_receiver_Sess_%d.jpg',Sess));
    saveas(fig,fname);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -----  Receiver Score
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    display([sprintf('Receiver score Session %d',Sess)])
    % -- parameters
    fmin = 10;
    fmax = 40;
    Delta_ms = 0; % time (ms) at which to start computing spectrogram. Default is 0 which means 0 ms, beginning of the time series.
    
    N = tapers(1);
    nt = 1000; % time length (ms) for the baseline in the stimulation experiment
    nwin = single(floor((nt-N*fs)./(dn*fs))); % numb of moving window (# points) in the resultig
    
    bin = (nt-N*fs)/nwin; % how many milliseconds correspon to a bsess_info{3}{i}in in the spectrogram
    Delta_bin = Delta_ms/bin;
    Nstep = (N*fs)/bin; % how many step correspond to a temporal shift of N
    % step = Nstep + nwin; % shift of N + nwin. To go from end of the window to mid of next window
    step = 110; % with steps < 100 there are overlapping windows for the computation of the beta power
    rscore = []; % receiver score
    for dt = (Delta_bin + nwin/2 +1):step:(size(spec_R,1)-(nwin+Nstep)/2) % from mid window + 1, to end of time - N/2 - nwin/2, with step of nwin + N
        rscore = [rscore, mean2(log(spec_R(dt:dt+nwin/2,fmin:fmax)))]; % mean of the log beta power in the range t = [500, 750] - keep in mind the spectrogram used a range [500,1000]
    end
    
    dt = (Delta_bin + nwin/2 +1):step:(size(spec_R,1)-(nwin+Nstep)/2); % BINS: from mid window + 1, to end of whole time - N/2 - nwin/2, with step of nwin + N
    ts = tiR(dt+(nwin+Nstep)/2); % Starting time (ms) for the computation of the cross-correlation sender-receiver
    
    dlmwrite(strcat(dir_Sess,'/receiver_score.txt'),[ts' rscore'],'delimiter','\t');
    
    
    mod_Ch = session_AM(i).mod_idx; % causal modulator channel
    cnt_m = 1;
    
    mscore = zeros(size(mod_Ch,2),271); % modulator score. 271 comes from the # of scores obtained from the for loop above 
    
    for Ch = mod_Ch % for all the modulators in the session
        
        
        % -- modulator spectrogram
        [specRS, fRS , tiRS] = tfspec_GINO(lfpRS(Ch,:),tapers,fs,dn,fk,pad,0.05,0,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----  Modulator Score
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        display([sprintf('Modulator score Ch %d',Ch)])
        
        ms = [];
        for dt = (Delta_bin + nwin/2 +1):step:(size(specRS,1)-(nwin+Nstep)/2) % from mid window + 1, to end of time - N/2 - nwin/2, with step of nwin + N
            ms = [ms mean2(log(specRS(dt:dt+nwin/2,fmin:fmax)))]; % mean of the log beta power in the range t = [500, 750] - keep in mind the spectrogram used a range [500,1000]
            %     ms = [ms, mean2(log(specRS(dt:dt+(nwin+Nstep)/2,fmin:fmax)))]; % mean of the log beta power in the range t = [500,1000] - keep in mind the spectrogram in this case used a range [500,1250]
        end
        
        mscore(cnt_m,:) = ms;
        
        % directory path to save files
        dir_Ch = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d/p_th_0.001/Ch_%d',Sess,Ch);
        dlmwrite(strcat(dir_Sess,sprintf('/modulator_score_ch_%d.txt',Ch)),[ts' mscore(cnt_m,:)'],'delimiter','\t');
        cnt_m = cnt_m + 1;
        
    end
    
    
    % -- receiver and modulator score --- %
    leg = cell(size(mod_Ch,2)+1,1); % dynamic legend 
    fig = figure;
    set(0,'DefaultLineLineWidth',2)
    plot(ts,rscore - mean(rscore)); hold on
    leg{1} = 'receiver score '
    for cnt_m = 1:size(mod_Ch,2)
        plot(ts,mscore(cnt_m,:) - mean(mscore(cnt_m,:))); hold on
        leg{cnt_m + 1} = sprintf('mod score ch %d',mod_Ch(cnt_m))
    end
    xlim([2000 50000])
    legend(leg,'FontSize',12)
    title(sprintf('receiver and modulator score Sess %d',Sess),'FontSize',12);
    grid on
    hold off
    set(gcf, 'Position',  [100, 600, 1200, 600])
    
    
    fname = strcat(dir_Sess,sprintf('/receiver_and_modulator_score_Sess_%d.png',Sess));
    saveas(fig,fname);
    
end
