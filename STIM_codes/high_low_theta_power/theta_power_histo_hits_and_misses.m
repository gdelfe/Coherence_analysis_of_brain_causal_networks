

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code considers the stimulation experiments. 
% 
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
dir_sort = strcat(dir_main,sprintf('/%s/Resting_state/%s/Modulators_controls',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

modulators = importdata(strcat(dir_sort,'/modulators_sorted_decod_accuracy_v2.txt')); % session, modulator idx, decod accuracy, order index i

in_name = '/mod_rec_stim_500_1000.mat';
% fig_name = '';
% title_name = 'Top 5 positive mod [500,1000]ms'

N = 5;

% time parameters: beginning and end of lfp signal used to determine the
% high/low modulator theta power for each trial 
t_i = 496;
t_f = 995;
t_tot = 500;


for n = 1:N
    
    Sess = modulators(n,1); % session
    dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));    
    load(strcat(dir_Stim_Sess,in_name));
    
    m = modulators(n,2);
    cnt_m = find(mod_rec_stim.mod_idx == m); % modulator
    
    if mod_rec_stim.receiver_idx ~= m
        
        % modulator's LFP
        lfp_M = sq(mod_rec_stim.lfp_E(:,m,:));
        
        W = 3;
        [spec, f, err] = dmtspec(lfp_M(:,t_i:t_f),[t_tot/1e3,W],1e3,200);
        theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
        
%         theta_pow_mean = mean(theta_pow); % get the average theta power
%         theta_pow = theta_pow - theta_pow_mean; % rescale the theta power by removing the mean value
        
        [sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
        cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution
        
        % low and high theta power indexes 
        low_theta_pow = sort_theta(1:cut);
        low_idx = trial_idx(1:cut);
        
        high_theta_pow = sort_theta(end-cut+1:end);
        high_idx = trial_idx(end-cut+1:end);
        
        hit = mod_rec_stim.hits;
        miss = mod_rec_stim.misses;
        
        a = mod_rec_stim.mod(cnt_m).confusion(1,1);
        b = mod_rec_stim.mod(cnt_m).confusion(1,2);
        c = mod_rec_stim.mod(cnt_m).confusion(2,1);
        d = mod_rec_stim.mod(cnt_m).confusion(2,2);
        D = a + d;
        ND = b + c;
        
        cut = length(high_idx);
        acc = (a + d)/2
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ROUGH ROC 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L = length(theta_pow);
        acc_list = [];
        for l = 1:L % varying the threshold for the ROC
            
            pot_miss = trial_idx(1:l); % everything at the left of the threshold
            pot_hit = trial_idx(l+1:end); % everything at the right of the threshold
            
            
            Tmiss = intersect(pot_miss,mod_rec_stim.misses);
            Thit = intersect(pot_hit,mod_rec_stim.hits);
            acc = (length(Tmiss) + length(Thit))/L;
            acc_list = [acc_list; [l, acc]];
        end
        [acc_max, idx_max] = max(acc_list(:,2)) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
        fig = figure;
        nbins = 30;
%         nbins_h = max( rount(length(theta_pow(hit)/nbins));
%         nbins_m = length(theta_pow(miss));
        h1 = histogram(theta_pow(hit),nbins,'FaceAlpha',0.5,'FaceColor','r');
        hold on 
        h2 = histogram(theta_pow(miss),nbins,'FaceAlpha',0.5,'FaceColor','b');
        grid on
        legend([h1,h2],{'hits','misses'},'AutoUpdate', 'off');
        title(sprintf('Theta power of modulator for hit and miss trials, top mod n. %d',n));
        xlabel('log(theta pow)');
        ylabel('count');
        set(gcf, 'Position',  [100, 600, 700, 500])
        ycoord = max(max(histcounts(theta_pow(hit),nbins)),max(histcounts(theta_pow(miss),nbins)));
        decod = mod_rec_stim.Decod_Accuracy(cnt_m);
        str = sprintf("Shaoyu dec acc = %.2f" + newline + "ROC acc = %.2f" +newline+ "CM acc = %.2f" + newline + "Conf. Diag. = %.2f" + newline + "Conf. non-Diag. = %.2f"+ newline + 'Conf. Mat. = ' + newline + '%.2f  %.2f'+ newline + '%.2f  %.2f',decod,acc_max,acc,D,ND,a,b,c,d);
        text(min(theta_pow),ycoord - 3,str);
        hold on;
        
        xline(low_theta_pow(end),'-',{'low threshold'}, 'LineWidth', 2, 'Color', 'b');
        xline(high_theta_pow(1),'-',{'high threshold'}, 'LineWidth', 2, 'Color', 'r');
        xline(sort_theta(idx_max),'-',{'ROC threshold'}, 'LineWidth', 2, 'Color', 'k');
        
        fname = strcat(dir_Stim_Theta,sprintf('/histo_hits_misses_top_m_%d.jpg',n));
        saveas(fig,fname);
      
    end
end


