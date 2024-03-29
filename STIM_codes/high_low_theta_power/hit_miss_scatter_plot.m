
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

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
dir_Stim_Theta = strcat(dir_main,sprintf('/%s/Stim_data/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


out_name = '/mod_rec_stim_500_1000.mat';
fig_name = '_500_1000ms';
% fig_caption = '[-1000,-500]ms';

c_mr_av = [];
theta_pow_HM = [];

for s = 1:size(sess_info{1},1) 

    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));
    
    load(strcat(dir_Stim_Sess,out_name));

    
    cnt_m = 1;
    for m = mod_rec_stim.mod_idx
        
        if mod_rec_stim.receiver_idx ~= m
            
%             miss_pow = log(mean(mod_rec_stim.mod(cnt_m).spec_m_miss((9:19))));
%             hit_pow = log(mean(mod_rec_stim.mod(cnt_m).spec_m_hit((9:19))));
            
            miss_pow = log(mean((sq(mod_rec_stim.Spec.psd.notDetected(:,5:9))),2));
            hit_pow = log(mean((sq(mod_rec_stim.Spec.psd.Detected(:,5:9))),2));
            
            theta_pow_HM = [theta_pow_HM; miss_pow, hit_pow];
            
          
%             c_mr_hit = mean(abs(mod_rec_stim.mod(cnt_m).c_mr_hit(9:19)));
%             c_mr_miss = mean(abs(mod_rec_stim.mod(cnt_m).c_mr_miss(9:19)));
            
%             c_mr_av = [c_mr_av; c_mr_miss, c_mr_hit];
             
        end
    
        cnt_m = cnt_m + 1;
    end 
    
end


fig = figure;
scatter(theta_pow_HM(:,1),theta_pow_HM(:,2),'filled');
grid on 
title('Hits vs miss mean theta power');
xlabel('misses log(mean(theta pow))');
ylabel('hits log(mean(theta pow))');

fname = strcat(dir_Stim_Theta,'/scatter_hits_vs_misses.jpg');
saveas(fig,fname);







