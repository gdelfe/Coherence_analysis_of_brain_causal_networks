
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between the causal 
% modulators found by Shaoyu's and both the sender and the receiver
%
% It computes the mean and the std of such coherence, across all the
% channels that have causal modulators
%
% INPUT: file with session modulator info
%        .mat file with structure AM and MA information 
%
% OUTPUT: txt files with the values of the coherence MR, MS, SR and
% corresponding figures
%
%    @ Gino Del Ferraro, August 2020, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/DL-modulators/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/DL-modulators/Gino_codes/shadedErrorBar-master')
dir_base = '/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State';
step = 110;

fid = fopen(strcat(dir_base,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

set(0,'DefaultLineLineWidth',2)

% -- load structure files
newAM = load(strcat(dir_base,'/session_AM.mat'))
session_AM = newAM.session_AM;

newMA = load(strcat(dir_base,'/session_MA.mat'))
session_MA = newMA.session_MA;


% -- print structures on stdout
%format short
for s=1:size(sess_info{1},1)
    session_AM(s)
    session_MA(s)
end

% --- generate control structure, i.e. non-modulaors 
for s=1:size(sess_info{1},1) % for a given session 
    channels = 1:session_AM(s).tot_ch;  % create array with all the channels label
    ch = session_AM(s).mod_idx ;          % all the modulator
    channels([ch]) = [];        % remove the modulators from the list
    control(s).mod_idx = channels(1:size(session_AM(s).mod_idx,2)); % make control same size as # modulator    
end


% -- matrices to store coherence for each mod and then compute
% average
ftot = 1638;
coh_ms_AM = zeros(48,ftot); % to store abs mean coherence mod-sender. # of causal mod, # of freq bins
coh_mr_AM = zeros(48,ftot); % to store abs mean coherence mod-receiver.# of causal mod, # of freq bins

coh_ms_MA = zeros(48,ftot); % mean abs choerence # of causal mod, # of freq bins
coh_mr_MA = zeros(48,ftot); % mean abs coherence # of causal mod, # of freq bins

fk = 200;
coh_sr_MA = importdata(strcat(dir_base,sprintf('/coherence_sr_MA_fk_%d.txt',fk)));
coh_sr_AM = importdata(strcat(dir_base,sprintf('/coherence_sr_AM_fk_%d.txt',fk)));

spec_M = zeros(48,ftot); % spectrum modulators
spec_S = zeros(48,ftot); % spectrum senders
spec_R = zeros(48,ftot); % spectrum receiver 

cnt = 1;
cnt_sr = 1;

% ---- parameters for the coherence-gram
tapers = [4 4];
N = tapers(1);
nt = tot_time;
dn = 0.01;
fs = 1000;
nwin = single(floor((nt-N*fs)/(dn*fs)));


for i=1:size(sess_info{1},1)  % For all the session with a modulator
    
    
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --------- COHERENCE-GRAM --------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    mod_Ch = control(i).mod_idx; % causal modulator channel
    
    for Ch = mod_Ch % for all the modulators in the session
        
        close all        
        
        display(['Computing coherence-gram...'])

        [c_ms,tf,f,spec_r,spec_s] = tfcoh_GINO(lfpRS(Ch,:),lfp_S,tapers,1e3,dn,fk,2,[],[],1); % coherence modulator-sender
        [c_mr,tf,f,spec_r,spec_s] = tfcoh_GINO(lfpRS(Ch,:),lfp_R,tapers,1e3,dn,fk,2,[],[],1); % coherence modulator-receiver
        
        %     dlmwrite(strcat(dir_Sess,'/sr_coherogram.txt'),c_sr,'delimiter',' ');
        % -- Figure: coherence spectrum
        
        % --- FIGURE: MODULATOR-SENDER COHEROGRAM
        fig_ms = figure; tvimage(abs(c_ms(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:100:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('M-S Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Sess,sprintf('/control_MS_coherogram_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig_ms,fname);
        
        % --- FIGURE: MODULATOR-RECEIVER COHEROGRAM
        fig_mr = figure; tvimage(abs(c_mr(:,:))); colorbar; % coherence spectrum
        xticks = floor(linspace(1,length(tf),5));
        xticklabels = tf(xticks);
        xtickformat('%d')
        yticks = 1:100:length(f);
        yticklabels = floor(f(yticks));
        ytickformat('%.2f')
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
        title(sprintf('M-R Coherogram, Sess = %d, ch = %d',Sess,Ch),'FontSize',12);
        xlabel('time (sec)');
        ylabel('freq (Hz)')
        % ylim([0,120])
        set(gcf, 'Position',  [100, 600, 1000, 600])
        
        fname = strcat(dir_Sess,sprintf('/control_MR_coherogram_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig_mr,fname);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MEAN ABS COHERENCE            %%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        % --- FIGURE --------- %%
        % -- Coherence vs frequency --- %
        fig = figure;
        plot(f,coh_sr_MA(cnt_sr,:))
        hold on
        plot(f,mean(abs(c_ms(:,:)),1))
        hold on
        plot(f,mean(abs(c_mr(:,:)),1))
        grid on
        title(sprintf('Mean Abs coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
        legend('S-R coherence','M-S coherence','M-R coherence')
%         xlim([0 60])
        set(gcf, 'Position',  [100, 600, 1000, 500])
    

        fname = strcat(dir_Sess,sprintf('/control_coherence_vs_freq_MA_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig,fname);
        
        coh_ms_MA(cnt,:) = mean(abs(c_ms(:,:)),1) ; % assign coherence M-S value for this modulator
        coh_mr_MA(cnt,:) = mean(abs(c_mr(:,:)),1) ; % assign coherence M-R value for this modulator
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ABS MEAN COHERENCE            %%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        % --- FIGURE --------- %%
        % -- Coherence vs frequency --- %
        fig = figure;
         plot(f,coh_sr_AM(cnt_sr,:))
        hold on
        plot(f,abs(mean(c_ms(:,:),1)))
        hold on
        plot(f,abs(mean(c_mr(:,:),1)))
        grid on
        title(sprintf('Abs Mean coherence vs frequency, ch = %d, causal mod',Ch),'FontSize',10);
        legend('S-R coherence','M-S coherence','M-R coherence')
%         xlim([0 60])
        set(gcf, 'Position',  [100, 600, 1000, 500])
        
        fname = strcat(dir_Sess,sprintf('/control_coherence_vs_freq_AM_ch_%d_fk_%d.jpg',Ch,fk));
        saveas(fig,fname);
        
        coh_ms_AM(cnt,:) = abs(mean(c_ms(:,:),1)) ; % assign coherence M-S value for this modulator
        coh_mr_AM(cnt,:) = abs(mean(c_mr(:,:),1)) ; % assign coherence M-R value for this modulator
        
        cnt = cnt + 1;
        
    end
    cnt_sr = cnt_sr + 1;
    
    
end

keyboard 
% % -- write data
% dlmwrite(strcat(dir_base,sprintf('/coherence_ms_MA_fk_%d.txt',fk)),coh_ms_MA,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/coherence_mr_MA_fk_%d.txt',fk)),coh_mr_MA,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/coherence_sr_MA_fk_%d.txt',fk)),coh_sr_MA,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/coherence_ms_AM_fk_%d.txt',fk)),coh_ms_AM,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/coherence_mr_AM_fk_%d.txt',fk)),coh_mr_AM,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/coherence_sr_AM_fk_%d.txt',fk)),coh_sr_AM,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/frequency_range_fk_%d.txt',fk)),f,'delimiter',' ');

% % -- write data: CONTROLS 
% dlmwrite(strcat(dir_base,sprintf('/control_coherence_ms_MA_fk_%d.txt',fk)),coh_ms_MA,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/control_coherence_mr_MA_fk_%d.txt',fk)),coh_mr_MA,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/control_coherence_ms_AM_fk_%d.txt',fk)),coh_ms_AM,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/control_coherence_mr_AM_fk_%d.txt',fk)),coh_mr_AM,'delimiter',' ');
% dlmwrite(strcat(dir_base,sprintf('/control_frequency_range_fk_%d.txt',fk)),f,'delimiter',' ');

% --- load data 
coh_ms_MA = importdata(strcat(dir_base,sprintf('/coherence_ms_MA_fk_%d.txt',fk)));
coh_mr_MA = importdata(strcat(dir_base,sprintf('/coherence_mr_MA_fk_%d.txt',fk)));
coh_sr_MA = importdata(strcat(dir_base,sprintf('/coherence_sr_MA_fk_%d.txt',fk)));
coh_ms_AM = importdata(strcat(dir_base,sprintf('/coherence_ms_AM_fk_%d.txt',fk)));
coh_mr_AM = importdata(strcat(dir_base,sprintf('/coherence_mr_AM_fk_%d.txt',fk)));
coh_sr_AM = importdata(strcat(dir_base,sprintf('/coherence_sr_AM_fk_%d.txt',fk)));
f = importdata(strcat(dir_base,sprintf('/frequency_range_fk_%d.txt',fk)));

% -- load data: CONTROLS 
coh_ms_AM_control = importdata(strcat(dir_base,sprintf('/control_coherence_ms_AM_fk_%d.txt',fk)));
coh_mr_AM_control = importdata(strcat(dir_base,sprintf('/control_coherence_mr_AM_fk_%d.txt',fk)));
coh_ms_MA_control = importdata(strcat(dir_base,sprintf('/control_coherence_ms_MA_fk_%d.txt',fk)));
coh_mr_MA_control = importdata(strcat(dir_base,sprintf('/control_coherence_mr_MA_fk_%d.txt',fk)));

% --- mean coherences
mean_cho_ms_MA = mean(coh_ms_MA(:,:));  % modulator - sender
mean_cho_mr_MA = mean(coh_mr_MA(:,:));  % modulator - receiver
mean_cho_sr_MA = mean(coh_sr_MA(:,:));  % sender - receiver 
mean_cho_ms_AM = mean(coh_ms_AM(:,:));
mean_cho_mr_AM = mean(coh_mr_AM(:,:));
mean_cho_sr_AM = mean(coh_sr_AM(:,:));

% --- mean coherences: CONTROLS
mean_cho_ms_AM_control = mean(coh_ms_AM_control(:,:));  % modulator - sender
mean_cho_mr_AM_control = mean(coh_mr_AM_control(:,:));  % modulator - receiver 
mean_cho_ms_MA_control = mean(coh_ms_MA_control(:,:));  % modulator - sender
mean_cho_mr_MA_control = mean(coh_mr_MA_control(:,:));  % modulator - receiver 

% --- std coherences
std_cho_ms_MA = std(coh_ms_MA(:,:));  % modulator - sender
std_cho_mr_MA = std(coh_mr_MA(:,:));  % modulator - receiver
std_cho_sr_MA = std(coh_sr_MA(:,:));  % sender - receiver
std_cho_ms_AM = std(coh_ms_AM(:,:));
std_cho_mr_AM = std(coh_mr_AM(:,:));
std_cho_sr_AM = std(coh_sr_AM(:,:));

% --- std coherences: CONTROLS 
std_cho_ms_AM_control = std(coh_ms_AM_control(:,:));  % modulator - sender
std_cho_mr_AM_control = std(coh_mr_AM_control(:,:));  % modulator - receiver
std_cho_ms_MA_control = std(coh_ms_MA_control(:,:));  % modulator - sender
std_cho_mr_MA_control = std(coh_mr_MA_control(:,:));  % modulator - receiver

% --- Error bars
err_ms_MA = std_cho_ms_MA/sqrt(48);
err_mr_MA = std_cho_mr_MA/sqrt(48);
err_sr_MA = std_cho_sr_MA/sqrt(20);
err_ms_AM = std_cho_ms_AM/sqrt(48);
err_mr_AM = std_cho_mr_AM/sqrt(48);
err_sr_AM = std_cho_sr_AM/sqrt(20);

% --- Error bars: CONTROLS 
err_ms_AM_control = std_cho_ms_AM_control/sqrt(48);
err_mr_AM_control = std_cho_mr_AM_control/sqrt(48);
err_ms_MA_control = std_cho_ms_MA_control/sqrt(48);
err_mr_MA_control = std_cho_mr_MA_control/sqrt(48);


fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
% hold all
% errorbar(f,mean_cho_ms_MA,err_ms_MA); hold on

% shadedErrorBar(f,mean_cho_ms_MA,err_ms_MA,'lineProps','b','patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_mr_MA,err_mr_MA,'lineprops',{'color',[255, 83, 26]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_sr_MA,err_sr_MA,'lineprops',{'color',[230, 184 , 0]/255},'patchSaturation',0.4); hold on

% shadedErrorBar(f,mean_cho_ms_AM,err_ms_AM,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_mr_AM,err_mr_AM,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_sr_AM,err_sr_AM,'lineprops',{'color',[0 204 204]/255},'patchSaturation',0.4); hold on
% 
% % -- added new 
% shadedErrorBar(f,mean_cho_ms_AM,err_ms_AM,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% shadedErrorBar(f,mean_cho_ms_AM_control,err_ms_AM_control,'lineprops',{'color',[201/255 79/255 138/255]},'patchSaturation',0.4); hold on

shadedErrorBar(f,mean_cho_mr_AM,err_mr_AM,'lineprops',{'color',[115, 230, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_cho_mr_AM_control,err_mr_AM_control,'lineprops',{'color',[51, 102, 0]/255},'patchSaturation',0.4); hold on

grid on
% title('Mean Abs coherence of MS, MR, SR','FontSize',11);
title('MR Abs Mean coherence MODULATORS and CONTROLS','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
% legend('M-S mean abs','M-R mean abs','M-S abs mean','M-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
% legend('M-S mean abs','M-R mean abs','S-R mean abs','FontSize',10)
% legend('M-S abs mean','M-R abs mean','S-R abs mean','FontSize',10)
legend('M-R abs mean causal mod','M-R abs mean control','FontSize',10)
% xlim([0 60])
set(gcf, 'Position',  [100, 600, 1100, 600])

% fname = strcat(dir_base,sprintf('/coherence_MA_MR_MS_SR_fk_%d.png',fk));
fname = strcat(dir_base,sprintf('/coherence_AM_MR_modulators_and_controls_fk_%d.png',fk));
saveas(fig,fname)


fname = strcat(dir_base,sprintf('/control_coherence_AM_mean_causal_MR_SR_fk_%d.png',fk));
fname = strcat(dir_base,sprintf('/control_coherence_MA_mean_causal_MR_SR_fk_%d_zoom.png',fk));
saveas(fig,fname)


% Prepare data

y=randn(30,80)*5;
x=(1:size(y,2))-40;
yP = sin( linspace(-2*pi,2*pi,length(x)) )*20;
y = bsxfun(@plus,y,yP)+60;

% Make the plot
clf
figure;
errorbar(f,mean_cho_ms_MA,err_ms_MA); hold on
shadedErrorBar(f,mean(coh_ms_MA),err_ms_MA);

% Overlay the raw data
hold on
plot(f,coh_ms_MA,'-','color',[0.5,0.5,0.95])

grid on


