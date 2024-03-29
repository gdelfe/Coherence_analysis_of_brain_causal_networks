

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
fig_name = '_500_1000_ms_TOP10_mod';
title_name = 'Top 10% mod [-500,0]ms'

N = 10;

c_mr_high = [];
c_mr_low = [];
c_mr_hit = [];
c_mr_miss = [];
spec_hit = [];
spec_miss = [];

cnt = 0;
for n = 1:N
    
    Sess = modulators(n,1); % session
    dir_Stim_Sess = strcat(dir_Stim_Theta,sprintf('/Sess_%d',Sess));    
    load(strcat(dir_Stim_Sess,in_name));
    
    m = modulators(n,2);
    cnt_m = find(mod_rec_stim.mod_idx == m); % modulator
    
    if mod_rec_stim.receiver_idx ~= m
        
        a = mod_rec_stim.mod(cnt_m).confusion(1,1);
        b = mod_rec_stim.mod(cnt_m).confusion(1,2);
        c = mod_rec_stim.mod(cnt_m).confusion(2,1);  
        d = mod_rec_stim.mod(cnt_m).confusion(2,2);
        D = a + d;
        ND = b + c;
        
%         if 2*(D - ND)/(D + ND) > 0.3; 
            display(['---------------',num2str(n)]);
            c_mr_high = [c_mr_high; mod_rec_stim.mod(cnt_m).c_mr_high];
            c_mr_low = [c_mr_low; mod_rec_stim.mod(cnt_m).c_mr_low];
            c_mr_hit = [c_mr_hit; mod_rec_stim.mod(cnt_m).c_mr_hit];
            c_mr_miss = [c_mr_miss; mod_rec_stim.mod(cnt_m).c_mr_miss];
            
            spec_hit = [spec_hit; mod_rec_stim.mod(cnt_m).spec_m_hit];
            spec_miss = [spec_miss; mod_rec_stim.mod(cnt_m).spec_m_miss];
            
            cnt = cnt + 1;
            mod_rec_stim.mod(cnt_m).confusion
            mod_rec_stim.Decod_Accuracy(cnt_m)
            display(['decod acc modulator -- ', num2str(modulators(n,3))]);

%         end
        
    end
    
end 


% coherence
mean_coh_mr_high = mean(abs(c_mr_high));
mean_coh_mr_low = mean(abs(c_mr_low)); 
mean_coh_mr_hit = mean(abs(c_mr_hit)); 
mean_coh_mr_miss = mean(abs(c_mr_miss)); 

err_coh_mr_high = std(abs(c_mr_high))/sqrt(size(c_mr_high,1)); 
err_coh_mr_low = std(abs(c_mr_low))/sqrt(size(c_mr_low,1));
err_coh_mr_hit = std(abs(c_mr_hit))/sqrt(size(c_mr_hit,1));
err_coh_mr_miss = std(abs(c_mr_miss))/sqrt(size(c_mr_miss,1));


% Spectrum 
mean_s_m_hit = mean(log(spec_hit)); 
mean_s_m_miss = mean(log(spec_miss)); 
err_s_m_hit = std(log(spec_hit))/sqrt(size(spec_hit,1));
err_s_m_miss = std(log(spec_miss))/sqrt(size(spec_miss,1));

f = mod_rec_stim.freq;
fspec = linspace(0,200,size(spec_hit,2));



%%%%%%%%%%%%%%%%%%%%
% FIGURES 
%%%%%%%%%%%%%%%%%%%%

% COHERENCE %%%%%%%%%%%%%%%%

fig = figure;

shadedErrorBar(f,mean_coh_mr_hit,err_coh_mr_hit,'lineprops',{'color',[204, 0, 0]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_coh_mr_miss,err_coh_mr_miss,'lineprops',{'color',[255, 128, 128]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title(sprintf('Coherence MR hits vs misses trials - %s',title_name),'FontSize',12)
legend('Hits','Misses','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_Stim_Theta,sprintf('/MR_hit_miss_coherence%s.jpg',fig_name));
saveas(fig,fname);



% SPECTRUM hits vs misses 

fig = figure;

shadedErrorBar(fspec,mean_s_m_hit,err_s_m_hit,'lineprops',{'color',[204, 0, 0]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(fspec,mean_s_m_miss,err_s_m_miss,'lineprops',{'color',[255, 128, 128]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',13);
ylabel('Spectrum','FontName','Arial','FontSize',13);
title(sprintf('Mean Spetrum of the modulators  - %s',title_name),'FontSize',12)
legend('Hits','Misses','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
% xlim([1 95])
% ylim([0 0.25])
grid on

fname = strcat(dir_Stim_Theta,sprintf('/spectrum_mod_hit_vs_miss%s.jpg',fig_name));
saveas(fig,fname);



