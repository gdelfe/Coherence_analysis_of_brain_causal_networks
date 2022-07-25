%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code compute plots the Sender-Receiver coherence for high/low theta
% power trials sorted according to:
% 1. Modulators' theta power
% 2. Controls (in the same region as the modulator, ctrl_SA)'s power
% 3. Controls (in other region but the modulator's, ctrl_OA)'s power
%
%    @ Gino Del Ferraro, June 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_HL_theta_Mav = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';
dir_HL_theta_Arc = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Archie/Resting_State/high_low_theta';
dir_both_monkeys = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys/high_low_theta';

% %%%%%%%%% Maverick High power %%%%%%
c_sr_mod_mav_H = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_high.mat')); % SR coherence - high theta - trial sorted by modulator power  
c_sr_mod_SA_mav_H = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_high_controls_SA.mat')); % SR coherence - high theta -  trials sorted by controls SA
c_sr_mod_OA_mav_H = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_high_controls_OA.mat')); % SR coherence - high theta - trials sorted by controls OA 

% %%%%%%%%% Maverick Low power %%%%%%
c_sr_mod_mav_L = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_low.mat')); % SR coherence - high theta - trial sorted by modulator power  
c_sr_mod_SA_mav_L = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_low_controls_SA.mat')); % SR coherence - high theta -  trials sorted by controls SA
c_sr_mod_OA_mav_L = load(strcat(dir_HL_theta_Mav,'/coh_all_sess_sr_low_controls_OA.mat')); % SR coherence - high theta - trials sorted by controls OA 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%  Archie High power %%%%%%
c_sr_mod_arc_H = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_high.mat')); % SR coherence - high theta - trial sorted by modulator power  
c_sr_mod_SA_arc_H = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_high_controls_SA.mat')); % SR coherence - high theta - trials sorted by controls SA
c_sr_mod_OA_arc_H = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_high_controls_OA.mat')); % SR coherence - high theta - trials sorted by controls OA

% %%%%%%%%% Archie Low power %%%%%%
c_sr_mod_arc_L = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_low.mat')); % SR coherence - high theta - trial sorted by modulator power  
c_sr_mod_SA_arc_L = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_low_controls_SA.mat')); % SR coherence - high theta -  trials sorted by controls SA
c_sr_mod_OA_arc_L = load(strcat(dir_HL_theta_Arc,'/coh_all_sess_sr_low_controls_OA.mat')); % SR coherence - high theta - trials sorted by controls OA 



% %%%% Concatenate results for Maverick and Archie 
% high power
c_sr_H = [c_sr_mod_mav_H.coh_all_c_sr_high; c_sr_mod_arc_H.coh_all_c_sr_high]; % concatenate data for the two monkeys - trial sorted by modulator power 
c_sr_SA_H = [c_sr_mod_SA_mav_H.coh_all_c_sr_high; c_sr_mod_SA_arc_H.coh_all_c_sr_high]; % concatenate data for the two monkeys - trial sorted by modulator power 
c_sr_OA_H = [c_sr_mod_OA_mav_H.coh_all_c_sr_high; c_sr_mod_OA_arc_H.coh_all_c_sr_high]; % concatenate data for the two monkeys - trial sorted by modulator power 
% low power
c_sr_L = [c_sr_mod_mav_L.coh_all_c_sr_low; c_sr_mod_arc_L.coh_all_c_sr_low]; % concatenate data for the two monkeys - trial sorted by modulator power 
c_sr_SA_L = [c_sr_mod_SA_mav_L.coh_all_c_sr_low; c_sr_mod_SA_arc_L.coh_all_c_sr_low]; % concatenate data for the two monkeys - trial sorted by modulator power 
c_sr_OA_L = [c_sr_mod_OA_mav_L.coh_all_c_sr_low; c_sr_mod_OA_arc_L.coh_all_c_sr_low]; % concatenate data for the two monkeys - trial sorted by modulator power 


% mean and SEM  ------------------------------
% high power 
mean_sr_H = mean(abs(c_sr_H));
mean_sr_SA_H = mean(abs(c_sr_SA_H));
mean_sr_OA_H = mean(abs(c_sr_OA_H));

err_sr_H = std(abs(c_sr_H))/sqrt(size(c_sr_H,1));
err_sr_SA_H = std(abs(c_sr_SA_H))/sqrt(size(c_sr_SA_H,1));
err_sr_OA_H = std(abs(c_sr_OA_H))/sqrt(size(c_sr_OA_H,1));

% low power
mean_sr_L = mean(abs(c_sr_L));
mean_sr_SA_L = mean(abs(c_sr_SA_L));
mean_sr_OA_L = mean(abs(c_sr_OA_L));

err_sr_L = std(abs(c_sr_L))/sqrt(size(c_sr_L,1));
err_sr_SA_L = std(abs(c_sr_SA_L))/sqrt(size(c_sr_SA_L,1));
err_sr_OA_L = std(abs(c_sr_OA_L))/sqrt(size(c_sr_OA_L,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fk = 200;
f = linspace(1,fk,409); % frequency values (range)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES MODULATORS vs CONTROLS %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all

% % trial sorted by modulators
shadedErrorBar(f,mean_sr_H,err_sr_H,'lineprops',{'color',[0, 102, 204]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_sr_L,err_sr_L,'lineprops',{'color',[128, 191, 255]/255},'patchSaturation',0.4); hold on

% % trial sorted by controls SA 
shadedErrorBar(f,mean_sr_SA_H,err_sr_SA_H,'lineprops',{'color',[102, 0, 204]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_sr_SA_L,err_sr_SA_L,'lineprops',{'color',[255, 102, 255]/255},'patchSaturation',0.4); hold on
% % % 
% % % % trial sorted by controls OA
shadedErrorBar(f,mean_sr_OA_H,err_sr_OA_H,'lineprops',{'color',[61, 92, 92]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,mean_sr_OA_L,err_sr_OA_L,'lineprops',{'color',[148, 184, 184]/255},'patchSaturation',0.4); hold on


grid on
title('Both animals: SR coherence for high/low power ctlr OA','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('high pow trial','low pow trial','FontSize',10)
set(gcf, 'Position',  [100, 600, 800, 500])
grid on
yticks([0 0.04 0.08 0.12 0.16 0.2 0.24])
xlim([0 95])
% ylim([0 0.24])

fname = strcat(dir_both_monkeys,'/coherency_SR_ctrl_OA_high_low_power_full.jpg');
saveas(fig,fname)
fname = strcat(dir_both_monkeys,'/coherency_SR_ctrl_OA_high_low_power_full.pdf');
saveas(fig,fname)


% coherence in theta range - high pow
theta_csr_H = mean_sr_H(:,9:20);
theta_csr_SA_H = mean_sr_SA_H(:,9:20);
theta_csr_OA_H = mean_sr_OA_H(:,9:20);

% theta coherence mean - high pow
theta_mean_csr_H = mean(theta_csr_H);
theta_mean_csr_SA_H = mean(theta_csr_SA_H);
theta_mean_csr_OA_H = mean(theta_csr_OA_H);

% SEM high power
err_csr_H = std(abs(c_sr_H(:,9:20)))/sqrt(size(c_sr_H,1));
err_mean_csr_H = sqrt(sum(err_csr_H.^2)/length(err_csr_H));

err_csr_SA_H = std(abs(c_sr_SA_H(:,9:20)))/sqrt(size(c_sr_SA_H,1));
err_mean_csr_SA_H = sqrt(sum(err_csr_SA_H.^2)/length(err_csr_SA_H));

err_csr_OA_H = std(abs(c_sr_OA_H(:,9:20)))/sqrt(size(c_sr_OA_H,1));
err_mean_csr_OA_H = sqrt(sum(err_csr_OA_H.^2)/length(err_csr_OA_H));

% ----------------------------------

% coherence in theta range - low pow
theta_csr_L = mean_sr_L(:,9:20);
theta_csr_SA_L = mean_sr_SA_L(:,9:20);
theta_csr_OA_L = mean_sr_OA_L(:,9:20);

% theta coherence mean - low pow
theta_mean_csr_L = mean(theta_csr_L);
theta_mean_csr_SA_L = mean(theta_csr_SA_L);
theta_mean_csr_OA_L = mean(theta_csr_OA_L);

% SEM low power
err_csr_L = std(abs(c_sr_L(:,9:20)))/sqrt(size(c_sr_L,1));
err_mean_csr_L = sqrt(sum(err_csr_L.^2)/length(err_csr_L));

err_csr_SA_L = std(abs(c_sr_SA_L(:,9:20)))/sqrt(size(c_sr_SA_L,1));
err_mean_csr_SA_L = sqrt(sum(err_csr_SA_L.^2)/length(err_csr_SA_L));

err_csr_OA_L = std(abs(c_sr_OA_L(:,9:20)))/sqrt(size(c_sr_OA_L,1));
err_mean_csr_OA_L = sqrt(sum(err_csr_OA_L.^2)/length(err_csr_OA_L));


diff_sr = theta_mean_csr_H - theta_mean_csr_L;
diff_sr_SA = theta_mean_csr_SA_H - theta_mean_csr_SA_L;
diff_sr_OA = theta_mean_csr_OA_H - theta_mean_csr_OA_L;

err_diff = sqrt(err_mean_csr_H.^2+err_mean_csr_L.^2);
err_diff_SA = sqrt(err_mean_csr_SA_H.^2+err_mean_csr_SA_L.^2);
err_diff_OA = sqrt(err_mean_csr_OA_H.^2+err_mean_csr_OA_L.^2);

data = [diff_sr,diff_sr_SA,diff_sr_OA];
err = [err_diff,err_diff_SA,err_diff_OA];

fig = figure;
bar(data)                
hold on
er = errorbar(data,err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
grid on 

hold off

fname = strcat(dir_both_monkeys,'/theta_choherence_diff_high_low.jpg');
saveas(fig,fname)
fname = strcat(dir_both_monkeys,'/theta_choherence_diff_high_low.pdf');
saveas(fig,fname)



