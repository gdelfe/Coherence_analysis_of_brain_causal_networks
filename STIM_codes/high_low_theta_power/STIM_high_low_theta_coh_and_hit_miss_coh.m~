
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

% set(0,'DefaultFigureVisible','off')
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
fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;

% time parameters: beginning and end of lfp signal used to determine the
% high/low modulator theta power for each trial 
t_i = 1;
t_f = 1000;

for s = 1:size(sess_info{1},1)  % For each session with at least one modulator

    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    
    load(strcat(dir_Stim_Sess,'/sess_data_stim.mat'));
    
    mod_rec_stim = sess_data_stim;
    hit = mod_rec_stim.hits;
    miss = mod_rec_stim.misses;
     
    % receiver idx and LFP
    r = mod_rec_stim.receiver_idx;
    lfp_R = sq(mod_rec_stim.lfp_E(:,r,:));
        
    cnt_m = 1;
    for m = mod_rec_stim.mod_idx
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Coherency Modulator-Receiver high/low theta         
        
        % modulator's LFP
        lfp_M = sq(mod_rec_stim.lfp_E(:,m,:));
        
        % high/low power indexes 
%         high_idx = mod_rec_stim.mod(cnt_m).high_idx;
%         low_idx = mod_rec_stim.mod(cnt_m).low_idx;
        
        
        % Compute the spectrum for each trial. Format: iTrial x times
        W = 3;
        [spec, f, err] = dmtspec(lfp_M(:,t_i:t_f),[1000/1e3,W],1e3,200);
        
        % Find low and high theta from the spectrum 
        theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
        theta_pow_mean = mean(theta_pow); % get the average theta power
        theta_pow = theta_pow - theta_pow_mean; % rescale the theta power by removing the mean value
        
        [sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
        cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution
        
        % low and high theta power indexes 
        low_theta_pow = sort_theta(1:cut);
        low_idx = trial_idx(1:cut);
        
        high_theta = sort_theta(end-cut+1:end);
        high_idx = trial_idx(end-cut+1:end);
        
        mod_rec_stim.mod(cnt_m).low_theta_pow = low_theta_pow;
        mod_rec_stim.mod(cnt_m).low_idx = low_idx;
        mod_rec_stim.mod(cnt_m).high_theta = high_theta;
        mod_rec_stim.mod(cnt_m).high_idx = high_idx;
        
        
        W = 5;
        % high/low theta power coherency
        [c_mr,f] = coherency(lfp_M,lfp_R,[N W],fs,fk,pad,0.05,1,1);
        [c_mr_high,f] = coherency(lfp_M(high_idx,:),lfp_R(high_idx,:),[N W],fs,fk,pad,0.05,1,1);
        [c_mr_low,f] = coherency(lfp_M(low_idx,:),lfp_R(low_idx,:),[N W],fs,fk,pad,0.05,1,1);
        
        % hits/misses coherences
        [c_mr_hit,f] = coherency(lfp_M(hit,:),lfp_R(hit,:),[N W],fs,fk,pad,0.05,1,1);
        [c_mr_miss,f] = coherency(lfp_M(miss,:),lfp_R(miss,:),[N W],fs,fk,pad,0.05,1,1);
        
        % high/low theta power coherency
        mod_rec_stim.mod(cnt_m).c_mr = c_mr;
        mod_rec_stim.mod(cnt_m).c_mr_high = c_mr_high;
        mod_rec_stim.mod(cnt_m).c_mr_low = c_mr_low;
        
        % hits/misses coherences
        mod_rec_stim.mod(cnt_m).c_mr_hit = c_mr_hit;
        mod_rec_stim.mod(cnt_m).c_mr_miss = c_mr_miss;
        
        mod_rec_stim.freq = f;
        
        cnt_m = cnt_m +1;
    end 
    
    save(strcat(dir_Stim_Sess,'/mod_rec_stim.mat'),'mod_rec_stim');
    
    clear mod_rec_stim
    
    
end 


keyboard;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd part, to be run independently
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c_mr_high = [];
c_mr_low = [];
c_mr_hit = [];
c_mr_miss = [];
for s = 1:size(sess_info{1},1) 

    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    
    load(strcat(dir_Stim_Sess,'/mod_rec_stim.mat'));
%     beg = beg.mod_rec_stim;
%     last = load(strcat(dir_Stim_Sess,'/mod_rec_stim_500_1000_ms_spec.mat'));
%     last = last.mod_rec_stim;
%     all = load(strcat(dir_Stim_Sess,'/mod_rec_stim.mat'));
%     all = all.mod_rec_stim;
% 
%     mod_rec_stim = beg;
    
    cnt_m = 1;
    for m = mod_rec_stim.mod_idx
        
        if mod_rec_stim.receiver_idx ~= m
        
            c_mr_high = [c_mr_high; mod_rec_stim.mod(cnt_m).c_mr_high];
            c_mr_low = [c_mr_low; mod_rec_stim.mod(cnt_m).c_mr_low];
            c_mr_hit = [c_mr_hit; mod_rec_stim.mod(cnt_m).c_mr_hit]; 
            c_mr_miss = [c_mr_miss; mod_rec_stim.mod(cnt_m).c_mr_miss];
            
%             b = beg.mod(cnt_m).low_idx;
%             l = last.mod(cnt_m).low_idx;
%             a = all.mod(cnt_m).low_idx;
%             common.sess(s).mod(cnt_m).low_idx_beg_last = setxor(b,l);
%             common.sess(s).mod(cnt_m).low_idx_beg_all = setxor(b,a);
%             common.sess(s).mod(cnt_m).low_idx_all_last = setxor(a,l);
%             
%             if ~isempty(setxor(b,l))
%                 display(['sess low',num2str(s)])
%             end
%             if ~isempty(setxor(b,a))
%                 display(['sess low',num2str(s)])
%             end
%             if ~isempty(setxor(a,l))
%                 display(['sess low',num2str(s)])
%             end
%             
%             b = beg.mod(cnt_m).high_idx;
%             l = last.mod(cnt_m).high_idx;
%             a = all.mod(cnt_m).high_idx;
%             common.sess(s).mod(cnt_m).high_idx_beg_last = setxor(b,l);
%             common.sess(s).mod(cnt_m).high_idx_beg_all = setxor(b,a);
%             common.sess(s).mod(cnt_m).high_idx_all_last = setxor(a,l);
%             
%             if ~isempty(setxor(b,l))
%                 display(['sess high',num2str(s)])
%             end
%             if ~isempty(setxor(b,a))
%                 display(['sess high',num2str(s)])
%             end
%             if ~isempty(setxor(a,l))
%                 display(['sess high',num2str(s)])
%             end
%             
        end
    
        cnt_m = cnt_m + 1;
    end 
    

end


keyboard;


mean_coh_mr_high = mean(abs(c_mr_high));
mean_coh_mr_low = mean(abs(c_mr_low)); 
mean_coh_mr_hit = mean(abs(c_mr_hit)); 
mean_coh_mr_miss = mean(abs(c_mr_miss)); 

err_coh_mr_high = std(abs(c_mr_high))/sqrt(size(c_mr_high,1)); 
err_coh_mr_low = std(abs(c_mr_low))/sqrt(size(c_mr_low,1));
err_coh_mr_hit = std(abs(c_mr_hit))/sqrt(size(c_mr_hit,1));
err_coh_mr_miss = std(abs(c_mr_miss))/sqrt(size(c_mr_miss,1));

f = mod_rec_stim.freq;

%%%%%%%%%%%%%%%%%%%%%%
% FIGURES %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% %%%% MR high power/low power %%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureVisible','on')

fig = figure;
shadedErrorBar(f,mean_coh_mr_high,err_coh_mr_high,'lineprops',{'color',[179, 71, 0]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_coh_mr_low,err_coh_mr_low,'lineprops',{'color',[255, 148, 77]/255 },'patchSaturation',0.5);


grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title('STIM Coherence MR high vs low theta power trial - high/low pow [-1000,0]ms','FontSize',10)
legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_Stim_Theta,sprintf('/MR_all_coherence_mean_stim_1000_ms.jpg',cnt_m));
saveas(fig,fname);



fig = figure;

shadedErrorBar(f,mean_coh_mr_hit,err_coh_mr_hit,'lineprops',{'color',[204, 0, 0]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_coh_mr_miss,err_coh_mr_miss,'lineprops',{'color',[255, 128, 128]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title('Coherence MR hits vs misses trials','FontSize',12)
legend('Hits','Misses','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_Stim_Theta,sprintf('/MR_hit_miss_coherence_mean_stim.jpg',cnt_m));
saveas(fig,fname);





