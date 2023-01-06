
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the SR coherence difference high-low for trial sorted
% with:  1. Modulators' power; 2. Ctrl-SA power; 3. Ctrl-OA power
% as a function of the modulator ranking.
%
% Results are shown across monkeys
% @ Gino Del Ferraro, Jan 2023, NYU

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_out = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys/high_low_theta';

freq_band = 'Theta_band';

dir_mod_ctrl_Arc = strcat(dir_main,sprintf('Archie/Resting_state/%s/Modulators_controls',freq_band));
dir_mod_ctrl_Mav = strcat(dir_main,sprintf('Maverick/Resting_state/%s/Modulators_controls',freq_band));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD SR COHERENCES HIGH/LOW power for 1. modulator's power, 2. Ctrl-SA
% power, 3. Ctrl-OA power
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%
% MODULATORS
% %%%%%%%%%%%%%%%%%%%

% Load SR coherence for high/low modulator power, ordered by ranking
% --- Maverick -----------------
monkey = 'Maverick';
dir_high_low_theta_mav = strcat(dir_main,sprintf('/%s/Resting_State/high_low_theta',monkey));
mav_coh_H = load(strcat(dir_high_low_theta_mav,'/coh_all_sess_sr_high_mod_ranked.mat')) % SR high coherence
mav_coh_L = load(strcat(dir_high_low_theta_mav,'/coh_all_sess_sr_low_mod_ranked.mat')) % SR low coherence 

% list of ranked modulators
dir_sorted = strcat(dir_main,sprintf('/%s/Resting_state/%s/Modulators_controls',monkey,freq_band));
mod_mav = load(strcat(dir_sorted,'/modulators_sorted_AUC.txt')); % session, mod id, rank/score, indx

% dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
% fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
% sess_info_mav = textscan(fid,'%d%s%s'); % sess label, date, RS label
% fclose(fid);

% --- Archie ---------------
monkey = 'Archie';
dir_high_low_theta_arc = strcat(dir_main,sprintf('/%s/Resting_State/high_low_theta',monkey));
arc_coh_H = load(strcat(dir_high_low_theta_arc,'/coh_all_sess_sr_high_mod_ranked.mat')) % SR high coherence
arc_coh_L = load(strcat(dir_high_low_theta_arc,'/coh_all_sess_sr_low_mod_ranked.mat')) % SR low coherence

% list of ranked modulators
dir_sorted = strcat(dir_main,sprintf('/%s/Resting_state/%s/Modulators_controls',monkey,freq_band));
mod_arc = load(strcat(dir_sorted,'/modulators_sorted_AUC.txt')); % session, mod id, rank/score, indx

% dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));
% fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
% sess_info_arc = textscan(fid,'%d%s%s'); % sess label, date, RS label
% fclose(fid);

% ---------------------

% %%%%%%%%%%%%%%%%%%%
% SAME AREA CONTROLS
% %%%%%%%%%%%%%%%%%%%

% ARCHIE
monkey = 'Archie';
arc_ctrl_list_SA = importdata(strcat(dir_mod_ctrl_Arc,'/control_list_same_area.txt')); % session, modulator idx, order index i

dir_high_low_theta = strcat(dir_main,sprintf('/%s/Resting_State/high_low_theta',monkey));
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_high_controls_SA.mat')) % SR coherence high according to control's power
coh_high_SA_arc = coh_all_c_sr_high;
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_low_controls_SA.mat'))
coh_low_SA_arc = coh_all_c_sr_low;

% MAVERICK
monkey = 'Maverick';
mav_ctrl_list_SA = importdata(strcat(dir_mod_ctrl_Mav,'/control_list_same_area.txt')); % session, modulator idx, order index i

dir_high_low_theta = strcat(dir_main,sprintf('/%s/Resting_State/high_low_theta',monkey));
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_high_controls_SA.mat'))
coh_high_SA_mav = coh_all_c_sr_high; 
clear coh_all_c_sr_high
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_low_controls_SA.mat'))
coh_low_SA_mav = coh_all_c_sr_low; 
clear coh_all_c_sr_low

% --------------------------------

% %%%%%%%%%%%%%%%%%%%
% OTHER AREAS CONTROLS
% %%%%%%%%%%%%%%%%%%%

% ARCHIE
monkey = 'Archie';

% list of controls OA Archie
arc_ctrl_list_OA = importdata(strcat(dir_mod_ctrl_Arc,'/control_list_other_areas.txt')); % session, modulator idx, order index i

dir_high_low_theta = strcat(dir_main,sprintf('/%s/Resting_State/high_low_theta',monkey));
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_high_controls_OA.mat'))
coh_high_OA_arc = coh_all_c_sr_high;
clear coh_all_c_sr_high
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_low_controls_OA.mat'))
coh_low_OA_arc = coh_all_c_sr_low;
clear coh_all_c_sr_low

% MAVERICK
monkey = 'Maverick';
% list of controls OA Maverick
mav_ctrl_list_OA = importdata(strcat(dir_mod_ctrl_Mav,'/control_list_other_areas.txt')); % session, modulator idx, order index i

dir_high_low_theta = strcat(dir_main,sprintf('/%s/Resting_State/high_low_theta',monkey));
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_high_controls_OA.mat'))
coh_high_OA_mav = coh_all_c_sr_high;
clear coh_all_c_sr_high
load(strcat(dir_high_low_theta,'/coh_all_sess_sr_low_controls_OA.mat'))
coh_low_OA_mav = coh_all_c_sr_low;
clear coh_all_c_sr_low

% --------------------------------- 



diff_array = []; err_arr = [];
diff_array_SA = []; err_arr_SA = [];
diff_array_OA = []; err_arr_OA = [];

percentage = [0.1,0.2,0.3,0.4,0.5,1];
for percent = percentage
    
    % %%%%%%%%%%%%%%%%%%
    % SR coherence difference, high-low modulators' power
    % %%%%%%%%%%%%%%%%%
    
    mav_idx = round(size(mav_coh_H.coh_all_c_sr_high,1)*percent) % amount of modulators in Maverick
    arc_idx = round(size(arc_coh_H.coh_all_c_sr_high,1)*percent) % amount of modulators in Archie
    
    % coherences for high and low modulator power for Maverick and Archie
    coh_high = [mav_coh_H.coh_all_c_sr_high(1:mav_idx,:); arc_coh_H.coh_all_c_sr_high(1:arc_idx,:)];
    coh_low = [mav_coh_L.coh_all_c_sr_low(1:mav_idx,:); arc_coh_L.coh_all_c_sr_low(1:arc_idx,:)];
    
    % coherence difference
    diff = abs(coh_high) - abs(coh_low);
    diff_mean = mean(diff);
    diff_mean_theta = mean(diff_mean(9:19)); % average difference over theta band
    diff_array = [diff_array, diff_mean_theta]; % store differences in array
    
    % coherence difference error
    diff_var = var(diff); % variance
    diff_var_theta = diff_var(9:19);
    diff_var_mean = mean(diff_var_theta); % mean of the variance in theta
    std_diff = sqrt(diff_var_mean);
    err_diff = std_diff/sqrt(11*(mav_idx+arc_idx)) % 11 is the number of time points in the theta frequency
    err_arr = [err_arr, err_diff];
    
    
    % %%%%%%%%%%%%%%%%%%
    % SR coherence difference, high-low ctrl-SA power
    % %%%%%%%%%%%%%%%%%
    
    % ARCHIE 
    % -- select sessions corresponding to the first N modulators 
    sess_numb = unique(mod_arc(1:arc_idx,1)); % session label with top modulators
    
    ctrl_idx = [];
    for i=1:length(sess_numb)
        ctrl_idx = [ctrl_idx; find(arc_ctrl_list_SA(:,1)==sess_numb(i))];
    end
    
    n_ctrl_arc = length(ctrl_idx); % number of controls
    
    coh_H_SA_arc = coh_high_SA_arc(ctrl_idx,:);
    coh_L_SA_arc = coh_low_SA_arc(ctrl_idx,:);

    % MAVERICK 
    % -- select sessions corresponding to the first N modulators 
    sess_numb = unique(mod_mav(1:mav_idx,1)); % session label with top modulators
    
    ctrl_idx = [];
    for i=1:length(sess_numb)
        ctrl_idx = [ctrl_idx; find(mav_ctrl_list_SA(:,1)==sess_numb(i))];
    end
    
    n_ctrl_mav = length(ctrl_idx); % number of controls 

    
    coh_H_SA_mav = coh_high_SA_mav(ctrl_idx,:);
    coh_L_SA_mav = coh_low_SA_mav(ctrl_idx,:);

    % concatenate coherences 
    coh_high_SA = [coh_H_SA_arc; coh_H_SA_mav];
    coh_low_SA = [coh_L_SA_arc; coh_L_SA_mav];
    
    % coherence difference
    diff = abs(coh_high_SA) - abs(coh_low_SA);
    diff_mean = mean(diff);
    diff_mean_theta = mean(diff_mean(9:19)); % average difference over theta band
    diff_array_SA = [diff_array_SA, diff_mean_theta]; % store differences in array
    
    % coherence difference error
    diff_var = var(diff); % variance
    diff_var_theta = diff_var(9:19);
    diff_var_mean = mean(diff_var_theta); % mean of the variance in theta
    std_diff = sqrt(diff_var_mean);
    err_diff = std_diff/sqrt(11*(n_ctrl_mav + n_ctrl_arc)) % 11 is the number of time points in the theta frequency
    err_arr_SA = [err_arr_SA, err_diff];

    % %%%%%%%%%%%%%%%%%%
    % SR coherence difference, high-low ctrl-OA power
    % %%%%%%%%%%%%%%%%%
    
    % ARCHIE 
    % -- select sessions corresponding to the first N modulators 
    sess_numb = unique(mod_arc(1:arc_idx,1)); % session label with top modulators
    
    ctrl_idx = [];
    for i=1:length(sess_numb)
        ctrl_idx = [ctrl_idx; find(arc_ctrl_list_OA(:,1)==sess_numb(i))];
    end
    
    n_ctrl_arc = length(ctrl_idx); % number of controls
    
    coh_H_OA_arc = coh_high_OA_arc(ctrl_idx,:);
    coh_L_OA_arc = coh_low_OA_arc(ctrl_idx,:);

    % MAVERICK 
    % -- select sessions corresponding to the first N modulators 
    sess_numb = unique(mod_mav(1:mav_idx,1)); % session label with top modulators
    
    ctrl_idx = [];
    for i=1:length(sess_numb)
        ctrl_idx = [ctrl_idx; find(mav_ctrl_list_OA(:,1)==sess_numb(i))];
    end
    
    n_ctrl_mav = length(ctrl_idx); % number of controls 

    
    coh_H_OA_mav = coh_high_OA_mav(ctrl_idx,:);
    coh_L_OA_mav = coh_low_OA_mav(ctrl_idx,:);

    % concatenate coherences 
    coh_high_OA = [coh_H_OA_arc; coh_H_OA_mav];
    coh_low_OA = [coh_L_OA_arc; coh_L_OA_mav];
    
    % coherence difference
    diff = abs(coh_high_OA) - abs(coh_low_OA);
    diff_mean = mean(diff);
    diff_mean_theta = mean(diff_mean(9:19)); % average difference over theta band
    diff_array_OA = [diff_array_OA, diff_mean_theta]; % store differences in array
    
    % coherence difference error
    diff_var = var(diff); % variance
    diff_var_theta = diff_var(9:19);
    diff_var_mean = mean(diff_var_theta); % mean of the variance in theta
    std_diff = sqrt(diff_var_mean);
    err_diff = std_diff/sqrt(11*(n_ctrl_mav + n_ctrl_arc)) % 11 is the number of time points in the theta frequency
    err_arr_OA = [err_arr_OA, err_diff];
    

    
end

 

fig = figure;
errorbar(percentage,diff_array,err_arr); hold on
errorbar(percentage,diff_array_SA,err_arr_SA); hold on
errorbar(percentage, diff_array_OA,err_arr_OA); 

xlim([0.05 1.05])
ylim([0 0.18])
xlabel('Modulator ranking','FontName','Arial','FontSize',10);
ylabel('SR coherence difference (high-low)','FontName','Arial','FontSize',10);
title('SR coherence difference (high-low)','FontName','Arial','FontSize',10);
legend('modulator','ctrl-SA','ctrl-OA')

fname = strcat(dir_out,'/SR_coherence_diff_high_low_vs_modulator_ranking.pdf');
saveas(fig,fname);




