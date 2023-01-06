


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_out = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys';

freq_band = 'theta_band';

iterations = 10000

% ARCHIE -----------------
monkey = 'Archie';
filename = '_rec001_002_all_sess'; % -- loading file name for coherence averages ******************
recording = 'rec001_002_all_sessions'; % -- folder where to load coherency files  *************

dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Controls = strcat(dir_RS,sprintf('/Modulators_Controls_avg_results/%s',recording));

load(strcat(dir_Controls,sprintf('/coh_spec_m_fk_200_W_5%s.mat',filename))); % structure mod
mod_arc = mod;
clear mod

% MAVERICK ----------------
monkey = 'Maverick';
filename = ''; % -- loading file name for coherence averages ******************
recording = 'last_recording'; % -- folder where to load coherency files  *************

dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_Controls = strcat(dir_RS,sprintf('/Modulators_Controls_avg_results/%s',recording));

load(strcat(dir_Controls,sprintf('/coh_spec_m_fk_200_W_5%s.mat',filename))); % structure mod
mod_mav = mod;
clear mod 

% concatenate Archie with Maverick
mod = struct;
mod.c_ms = {mod_mav.c_ms mod_arc.c_ms};
mod.c_mr = {mod_mav.c_mr mod_arc.c_mr};

% Structure to randomly select pseudo MS and MR
pseudo_coh = struct;
pseudo_coh = {mod.c_ms{:} mod.c_mr{:}};

diff = [];
L = length(mod.c_ms);
for iter = 1:iterations
    
    perm = randperm(L);
    pseudo_diff = abs(pseudo_coh{perm(1)})- abs(pseudo_coh{perm(2)});
    diff = [diff; pseudo_diff];
    
end 


save(strcat(dir_out,'/null_distribution_coh_MR_vs_MS_difference.mat'),'diff')


figure;
histogram(diff(:,10),50)

c_ms = []; c_mr = [];
for i=1:length(mod.c_ms)
    c_ms = [c_ms; mod.c_ms{i}];
    c_mr = [c_mr; mod.c_ms{i}];
end

mean_cms = mean(abs(c_ms));
mean_cmr = mean(abs(c_mr));

% INCOMPLETE


