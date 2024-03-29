
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code considers the stimulation experiments. 
% It computes the MR coherence for high/low modulator theta power trials
% NOTE: chosing the trials length for the definition of high/low theta
% matters: [0,500] ms in the baseline do not produce any significant
% difference, whereas [0,1000]ms does.
%
%    @ Gino Del Ferraro, March 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')c
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - PATHS and NAMES --- %
%%%%%%%%%%%%%%%%%%%
addpath('/mnt/pesaranlab/Matlab/monkeys')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
dir_high_low_theta = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_Stim_Theta = strcat(dir_main,sprintf('/%s/Stim_data/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


% ---- parameters for the coherence-gram
tapers = [0.5 7];
k = 2*tapers(1)*tapers(2)-1;
N = tapers(1);
fs = 1000;
fk = 200;
dn = 0.01;
nt = 1000;
nwin = floor((nt-N*fs)./(dn*fs));          % calculate the number of windows

% time parameters: beginning and end of lfp signal used to determine the
% high/low modulator theta power for each trial 
t_i = 496;
t_f = 995;
t_tot = 500;
out_name = '/mod_rec_stim_gram.mat'
img_out = '_1_1000';
fig_caption = '[-1000, 0]ms';

% for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
% 
%     Sess = sess_info{1}(s); % Session number
%     display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
%     dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
%     
%     load(strcat(dir_Stim_Sess,'/sess_data_stim.mat'));
%     
%     mod_rec_stim_gram = sess_data_stim;
%     hit = mod_rec_stim_gram.hits;
%     miss = mod_rec_stim_gram.misses;
%      
%     % receiver idx and LFP
%     r = mod_rec_stim_gram.receiver_idx;
%     lfp_R = sq(mod_rec_stim_gram.lfp_E(:,r,:));
%         
%     cnt_m = 1;
%     for m = mod_rec_stim_gram.mod_idx
%     
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % % Coherency Modulator-Receiver high/low theta         
%         
%         % modulator's LFP
%         lfp_M = sq(mod_rec_stim_gram.lfp_E(:,m,:));
%         
%         % high/low power indexes 
% %         high_idx = mod_rec_stim.mod(cnt_m).high_idx;
% %         low_idx = mod_rec_stim.mod(cnt_m).low_idx;
%         
%         
%         % Compute the spectrum for each trial. Format: iTrial x times
%         W = 3;
%         [spec, f, err] = dmtspec(lfp_M(:,t_i:t_f),[t_tot/1e3,W],1e3,200);
%         
%         % Find low and high theta from the spectrum 
%         theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
%         theta_pow_mean = mean(theta_pow); % get the average theta power
%         theta_pow = theta_pow - theta_pow_mean; % rescale the theta power by removing the mean value
%         
%         [sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
%         cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution
%         
%         % low and high theta power indexes 
%         low_theta_pow = sort_theta(1:cut);
%         low_idx = trial_idx(1:cut);
%         
%         high_theta = sort_theta(end-cut+1:end);
%         high_idx = trial_idx(end-cut+1:end);
%         
%         mod_rec_stim_gram.mod(cnt_m).low_theta_pow = low_theta_pow;
%         mod_rec_stim_gram.mod(cnt_m).low_idx = low_idx;
%         mod_rec_stim_gram.mod(cnt_m).high_theta = high_theta;
%         mod_rec_stim_gram.mod(cnt_m).high_idx = high_idx;
%         
%         
%         W = 5;
%         % high/low theta power coherency
%       
%         
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % COHERENCE-GRAM
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         % --- coherence-gram 
%         [c_mr,f,spec_m,spec_r] = tfcoh(lfp_M,lfp_R,tapers,1e3,dn,200,2,[],[],1);
%         [c_mr_high,f, spec_m_high, spec_r_high] = tfcoh(lfp_M(high_idx,:),lfp_R(high_idx,:),tapers,1e3,dn,200,2,[],[],1); % high power trials
%         [c_mr_low,f, spec_m_low, spec_r_low] = tfcoh(lfp_M(low_idx,:),lfp_R(low_idx,:),tapers,1e3,dn,200,2,[],[],1); % low power trials 
%         
% %         % separate calculation of the spectrogram
% %         k = 4;
% %         fk = [0 60];
% %         tapers = [0.5 5];
% %         dn = 0.005;
% %         pad = 2;
% %        
% %         [specRS, fRS , tiRS] = tfspec(lfp_M,tapers,fs,dn,fk,pad,0.05,1);
%         
%         
% %         figure;
% %         tvimage(abs(c_mr)); colorbar;
% %         
% %         xticks = floor(linspace(1,size(c_mr,1),5));     
% %         xticklabels = floor(linspace(1,nt,5));
% %         xtickformat('%d')
% %         % yticks = floor(linspace(1,length(f),10));
% %         yticks = 1:20:length(f);
% %         
% %         yticklabels = floor(f(yticks));
% %         ytickformat('%.2f')
% %         set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
% %         title('Coherence-gram S-R','FontSize',12);
% %         xlabel('time (sec)');
% %         ylabel('freq (Hz)')
% %         ylim([0,120])
% %         set(gcf, 'Position',  [100, 600, 1000, 600])
%         
% 
% %         figure;
% %         tvimage(log(sq(specRS))); colorbar;        
% %         xticks = floor(linspace(1,size(specRS,1),5));     
% %         xticklabels = floor(linspace(1,nt,5));
% %         xtickformat('%d')
% %         yticks = floor(linspace(1,length(f),10));
% %         yticks = 1:20:length(f);
% %         
% %         yticklabels = floor(f(yticks));
% %         ytickformat('%.2f')
% %         set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
% %         title('Coherence-gram S-R','FontSize',12);
% %         xlabel('time (sec)');
% %         ylabel('freq (Hz)')
% %         ylim([0,120])
% %         set(gcf, 'Position',  [100, 600, 1000, 600])
% 
%         % hits/misses coherences
%         [c_mr_hit,f, spec_m_hit, spec_r_hit] = tfcoh(lfp_M(hit,:),lfp_R(hit,:),tapers,1e3,dn,200,2,[],[],1);
%         [c_mr_miss,f, spec_m_miss, spec_r_miss] = tfcoh(lfp_M(miss,:),lfp_R(miss,:),tapers,1e3,dn,200,2,[],[],1);
%         
%         % high/low theta power coherence-gram
%         mod_rec_stim_gram.mod(cnt_m).c_mr = c_mr;
%         mod_rec_stim_gram.mod(cnt_m).c_mr_high = c_mr_high;
%         mod_rec_stim_gram.mod(cnt_m).c_mr_low = c_mr_low;
%         % spectrum 
%         mod_rec_stim_gram.mod(cnt_m).s_m = spec_m;
%         mod_rec_stim_gram.mod(cnt_m).s_r = spec_r;
%         mod_rec_stim_gram.mod(cnt_m).s_m_high = spec_m_high;
%         mod_rec_stim_gram.mod(cnt_m).s_r_high = spec_r_high;
%         mod_rec_stim_gram.mod(cnt_m).s_m_low = spec_m_low;
%         mod_rec_stim_gram.mod(cnt_m).s_r_low = spec_r_low;
% 
%         
%         % hits/misses coherence-gram
%         mod_rec_stim_gram.mod(cnt_m).c_mr_hit = c_mr_hit;
%         mod_rec_stim_gram.mod(cnt_m).c_mr_miss = c_mr_miss;
%         % spectrum hit/miss
%         mod_rec_stim_gram.mod(cnt_m).s_m_hit = spec_m_hit;
%         mod_rec_stim_gram.mod(cnt_m).s_r_hit = spec_r_hit;
%         mod_rec_stim_gram.mod(cnt_m).s_m_miss = spec_m_miss;
%         mod_rec_stim_gram.mod(cnt_m).s_r_miss = spec_r_miss;
%         
%         mod_rec_stim_gram.freq = f;
%         
%         cnt_m = cnt_m +1;
%     end 
%     
%     save(strcat(dir_Stim_Sess,out_name),'mod_rec_stim_gram');
%     
%     clear mod_rec_stim_gram
%     
%     
% end 


keyboard;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd part, to be run independently
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_m_high = [];
s_m_low = [];
s_m_hit = [];
s_m_miss = [];

c_mr_high = []; 
c_mr_low = [];
c_mr_hit = [];
c_mr_miss = [];
for s = 1:size(sess_info{1},1) 

    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    
    load(strcat(dir_Stim_Sess,out_name));

    cnt_m = 1;
    for m = mod_rec_stim_gram.mod_idx
        
        if mod_rec_stim_gram.receiver_idx ~= m
        
            
            c_mr_high = cat(3, c_mr_high,abs(mod_rec_stim_gram.mod(cnt_m).c_mr_high));
            c_mr_low = cat(3, c_mr_low, abs(mod_rec_stim_gram.mod(cnt_m).c_mr_low));
            c_mr_hit = cat(3, c_mr_hit, abs(mod_rec_stim_gram.mod(cnt_m).c_mr_hit)); 
            c_mr_miss = cat(3, c_mr_miss, abs(mod_rec_stim_gram.mod(cnt_m).c_mr_miss));
            
            s_m_high = cat(3, s_m_high,abs(mod_rec_stim_gram.mod(cnt_m).s_m_high));
            s_m_low = cat(3, s_m_low, abs(mod_rec_stim_gram.mod(cnt_m).s_m_low));
            s_m_hit = cat(3, s_m_hit, abs(mod_rec_stim_gram.mod(cnt_m).s_m_hit)); 
            s_m_miss = cat(3, s_m_miss, abs(mod_rec_stim_gram.mod(cnt_m).s_m_miss));
            
            
            
        end
    
        cnt_m = cnt_m + 1;
    end 
    

end


keyboard;

% compute the mean of the coherence-gram across modulators-receiver pairs 
mean_coh_mr_high = mean(c_mr_high,3);
mean_coh_mr_low = mean(c_mr_low,3); 
mean_coh_mr_hit = mean(c_mr_hit,3);
mean_coh_mr_miss = mean(c_mr_miss,3);

mean_s_m_high = mean(s_m_high,3);
mean_s_m_low = mean(s_m_low,3); 
mean_s_m_hit = mean(s_m_hit,3);
mean_s_m_miss = mean(s_m_miss,3);


f = mod_rec_stim_gram.freq;

%%%%%%%%%%%%%%%%%%%%%%
% FIGURES %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% %%%% MR coherence-gram high power %%%%%%%%%%%%%%%%%%%%%%

fig = figure;
tvimage(mean_coh_mr_high); colorbar;

ax = gca;
clim = ax.CLim;
xticks = floor(linspace(1,size(mean_coh_mr_high,1),5));
xticklabels = floor(linspace(1,nt,5));
xtickformat('%d')
% yticks = floor(linspace(1,length(f),10));
yticks = 1:20:length(f);

yticklabels = floor(f(yticks));
ytickformat('%.2f')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
title(sprintf('Mean coherence-gram (across mod) for high theta pow trials - %s',fig_caption),'FontSize',10);
xlabel('time (sec)');
ylabel('freq (Hz)')
ylim([0,120])
set(gcf, 'Position',  [100, 600, 1000, 600])


fname = strcat(dir_Stim_Theta,sprintf('/coherence-gram_high_pow_MR%s.jpg',img_out));
saveas(fig,fname);



% %%%% MR coherence-gram low power %%%%%%%%%%%%%%%%%%%%%%


fig = figure;
tvimage(mean_coh_mr_low); colorbar;

ax = gca;
ax.CLim = clim;
xticks = floor(linspace(1,size(mean_coh_mr_low,1),5));
xticklabels = floor(linspace(1,nt,5));
xtickformat('%d')
% yticks = floor(linspace(1,length(f),10));
yticks = 1:20:length(f);

yticklabels = floor(f(yticks));
ytickformat('%.2f')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
title(sprintf('Mean coherence-gram (across mod) for low theta pow trials - %s',fig_caption),'FontSize',10);
xlabel('time (sec)');
ylabel('freq (Hz)')
ylim([0,120])
set(gcf, 'Position',  [100, 600, 1000, 600])


fname = strcat(dir_Stim_Theta,sprintf('/coherence-gram_low_pow_MR%s.jpg',img_out));
saveas(fig,fname);



% %%%% MR spectrogram high power %%%%%%%%%%%%%%%%%%%%%%


fig = figure;
tvimage(log(mean_s_m_high)); colorbar;

ax = gca;
clim = ax.CLim;
xticks = floor(linspace(1,size(mean_s_m_high,1),5));
xticklabels = floor(linspace(1,nt,5));
xtickformat('%d')
% yticks = floor(linspace(1,length(f),10));
yticks = 1:20:length(f);

yticklabels = floor(f(yticks));
ytickformat('%.2f')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
title(sprintf('Mean spectrogram (across mod) for high theta pow trials -  %s',fig_caption),'FontSize',10);
xlabel('time (sec)');
ylabel('freq (Hz)')
ylim([0,120])
set(gcf, 'Position',  [100, 600, 1000, 600])


fname = strcat(dir_Stim_Theta,sprintf('/spectrogram_high_pow_M%s.jpg',img_out));
saveas(fig,fname);


% %%%% MR spectrogram low power %%%%%%%%%%%%%%%%%%%%%%


fig = figure;
tvimage(log(mean_s_m_low)); colorbar;

ax = gca;
clim = ax.CLim;
xticks = floor(linspace(1,size(mean_s_m_low,1),5));
xticklabels = floor(linspace(1,nt,5));
xtickformat('%d')
% yticks = floor(linspace(1,length(f),10));
yticks = 1:20:length(f);

yticklabels = floor(f(yticks));
ytickformat('%.2f')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'YTick', yticks, 'YTickLabel', yticklabels)
title(sprintf('Mean spectrogram (across mod) for low theta pow trials -  %s',fig_caption),'FontSize',10);
xlabel('time (sec)');
ylabel('freq (Hz)')
ylim([0,120])
set(gcf, 'Position',  [100, 600, 1000, 600])


fname = strcat(dir_Stim_Theta,sprintf('/spectrogram_low_pow_M%s.jpg',img_out));
saveas(fig,fname);





