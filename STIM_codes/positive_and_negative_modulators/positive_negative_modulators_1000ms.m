
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code considers the stimulation experiments. For each trial, it
% checks whether the receiver's response is associated with high-(low)-theta
% power of the modulator's electrode. Modulators with high-theta power for  
% receiver response are considered "positive", modulators with low theta 
% power for receiver response are considered "negative". 
%
% NOTE: in this code we use LFP baseline length 1000 ms before stimulation 
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
dir_out = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_Stim_Theta = strcat(dir_main,sprintf('/%s/Stim_data/%s',monkey,freq_band));
dir_out = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys'

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

name_structure = '/modulators_decod_accuracy.mat';

% time parameters for 500 ms test
t_tot = 500
ti = 496; 
tf = 995;


% positive and negative modulators 
positive = 0;
negative = 0;
tot_m  = 0; % total number of modulators across sessions
tot_conf = [0 0; 0 0];

for s = 1:size(sess_info{1},1)  % For each session with at least one modulator

    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));


    load(strcat(dir_Stim_Sess,'/sess_data_stim.mat'));
    
    hit = sess_data_stim.hits;
    miss = sess_data_stim.misses;

    lfp_E = sess_data_stim.lfp_E;
    cnt_m = 1;
    for m = sess_data_stim.mod_idx
        
        display(['-- Modulator ',num2str(m),' --  ',num2str(cnt_m),', out of tot  ',num2str(length(sess_data_stim.mod_idx)),' '])

        lfp_m = sq(lfp_E(:,m,:)); % modulator lfp
        % Compute the spectrum for each trial. Format: iTrial x times
        W = 3; % frequency smoothing 
%         [spec, f,err] = dmtspec(lfp_m,[1000/1e3,W],1e3,200); % spectrum 1000 ms before onset 
        [spec, f,err] = dmtspec(lfp_m(:,ti:tf),[t_tot/1e3,W],1e3,200); % spectrum 500 ms before onset

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
        
        HH = intersect(hit,high_idx); % high power - hit trials, true positive
        ML = intersect(miss,low_idx); % low power - miss trials, true negative
        HL = intersect(hit,low_idx); % low power - miss trials, false negative
        MH = intersect(miss,high_idx); % high power - hit trials, false positive
        
        HH_rate = length(HH)/cut;
        ML_rate = length(ML)/cut;
        HL_rate = length(HL)/cut;
        MH_rate = length(MH)/cut;
        
        confusion = [HH_rate, HL_rate; MH_rate, ML_rate]; % confusion matrix 
        tot_conf = tot_conf + confusion;
        % store theta power, idx, and confusion matrix for each modulator 
        sess_data_stim.mod(cnt_m).theta_pow = theta_pow;
        sess_data_stim.mod(cnt_m).low_theta_pow = low_theta_pow;
        sess_data_stim.mod(cnt_m).low_idx = low_idx;
        sess_data_stim.mod(cnt_m).high_theta = high_theta;
        sess_data_stim.mod(cnt_m).high_idx = high_idx;
        sess_data_stim.mod(cnt_m).confusion = confusion;
        
        if sess_data_stim.receiver_idx ~= m % exclude modulators which are receivers 
            % count number of positive and negative modulators
            
            if HH_rate + ML_rate > HL_rate + MH_rate
                positive = positive + 1;
            else
                negative = negative + 1;
            end 
            tot_m = tot_m + 1;
        end
        
        cnt_m = cnt_m + 1;

    end %  for each modulator within session 
      
    save(strcat(dir_Sess,'/sess_data_stim.mat'),'sess_data_stim');
    clear sess_data_stim
    
end  % for each session

tot_conf = tot_conf./tot_m;

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COUNT number of positive/negative modulators 
% with receivers which are not modulators themselves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tot_m = 0;
% positive = 0;
% negative = 0;
% for s = 1:size(sess_info{1},1)
% 
%     Sess = sess_info{1}(s); % Session number
%     display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
%     dir_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
%     
%     load(strcat(dir_Sess,'/sess_data_stim.mat'));
%     
%     cnt_m = 1;
%     for m = sess_data_stim.mod_idx
%                 
%         if sess_data_stim.receiver_idx ~= m % exclude modulators which are receivers 
%             % count number of positive and negative modulators
%             
%             HH_rate = sess_data_stim.mod(cnt_m).confusion(1,1);
%             ML_rate = sess_data_stim.mod(cnt_m).confusion(2,2);
%             HL_rate = sess_data_stim.mod(cnt_m).confusion(1,2);
%             MH_rate = sess_data_stim.mod(cnt_m).confusion(2,1);
%             
%             if HH_rate + ML_rate > HL_rate + MH_rate
%                 positive = positive + 1;
%             else
%                 negative = negative + 1;
%             end 
%             tot_m = tot_m + 1;
%         end
%         
%         cnt_m = cnt_m + 1;
%     end 
%     
% 
% end 

keyboard 

score = [];
score_hit = [];
for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    
    load(strcat(dir_Sess,'/sess_data_stim.mat'));
%     load(strcat(dir_Sess,'/Data_with_theta_band.mat')); % load stim data for info about hits/misses
    

    % hit and miss trials -- they are the same for each modulator within the same session
%     hit = Data.spec.lfp.DetectedIndx{1};
%     miss = Data.spec.lfp.notDetectedIndx{1};
    
%     sess_data_stim.hits = hit;
%     sess_data_stim.misses = miss;
        

    cnt_m = 1;
    for m = sess_data_stim.mod_idx
        
        if sess_data_stim.receiver_idx ~= m % exclude modulators which are receivers
            % count number of positive and negative modulators
            
            HH_rate = sess_data_stim.mod(cnt_m).confusion(1,1);
            ML_rate = sess_data_stim.mod(cnt_m).confusion(2,2);
            HL_rate = sess_data_stim.mod(cnt_m).confusion(1,2);
            MH_rate = sess_data_stim.mod(cnt_m).confusion(2,1);
            
            conf_score = HH_rate + ML_rate - HL_rate - MH_rate;
            AUC_score = sess_data_stim.auc(cnt_m);
            dec_score = sess_data_stim.Decod_Accuracy(cnt_m);
            score = [score; conf_score, AUC_score, dec_score];
%             score_hit = [score_hit; HH_rate, length(hit)];
            
        end
    end % for each modulators
     
%     save(strcat(dir_Sess,'/sess_data_stim.mat'),'sess_data_stim');

end
   
keyboard 

scores.conf_score = score(:,1);
scores.AUC_score = score(:,2);
scores.dec_score = score(:,3);
scores.HH_rate = score_hit(:,1);
scores.n_hits = score_hit(:,2);

save(strcat(dir_out,'/modulators_scores_Archie.mat'),'scores');
% load(strcat(dir_high_low_theta,'/modulators_scores.mat'))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AUC vs Confusion score
fig = figure;
scatter(score(:,1),score(:,2),'filled','b');
grid on 
ylabel('AUC score')
xlabel('Confusion matrix score')
mdl = fitlm(score(:,1),score(:,2))
R2 = mdl.Rsquared.Ordinary
str = sprintf('R2 = %.2f',R2);
text(0,0.6,str)

fname = strcat(dir_high_low_theta,'/scatter_score_AUC.jpg');
saveas(fig,fname);


% Decoding accuracy vs Confusion score 
fig = figure;
scatter(score(:,1),score(:,3),'filled','r');
grid on 
ylabel('Decoding accuracy score')
xlabel('Confusion matrix score')
mdl = fitlm(score(:,1),score(:,3))
R2 = mdl.Rsquared.Ordinary
str = sprintf('R2 = %.2f',R2);
text(0,0.75,str)

fname = strcat(dir_high_low_theta,'/scatter_score_decod_acc.jpg');
saveas(fig,fname);


% HH_score vs numb of hits 
fig = figure;
scatter(score_hit(:,1),score_hit(:,2),'filled','k');
grid on 
xlabel('High pow - Hits rate')
ylabel('# of hits trials')
mdl = fitlm(score(:,1),score(:,3))
R2 = mdl.Rsquared.Ordinary
str = sprintf('R2 = %.2f',R2);
text(0.3,120,str)

fname = strcat(dir_high_low_theta,'/scatter_HH_rate_vs_hits.jpg');
saveas(fig,fname);



