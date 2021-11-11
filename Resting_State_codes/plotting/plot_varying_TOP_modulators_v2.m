
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_fig = strcat(dir_main,'Figures_paper');


freq_band = 'theta_band';
recording_both = 'last_rec-rec001_002';
decoding = 'AUC';

fk = 200; W = 5;
f = linspace(1,fk,409); % frequency values (range)

% -- beta
% N_list = [10,20,30,40,50,100];
% N_mav_list = [4,8,12,16,20,41];
% N_arc_list = [5,10,15,20,25,51];

% -- theta
% N_list = [10,20,30,40,50,100];
N_list = 100;

% % ranges at which the peak occours for the modulators, ctrl SA, ctrl OA for MR and CR coherence 
% flist_mod_mr = {11:17,10:16,12:14,13:15,12:14,13:15};
% flist_SA_mr = {19:21,16:18,14:16,15:17,14:16,14:16};
% flist_OA_mr = {17:19,14:16,14:16,14:16,13:15,13:15};

flist_mod_mr = {7:15};
flist_SA_mr = {7:15};
flist_OA_mr = {7:15};

% % ranges at which the peak occours for the modulators, ctrl SA, ctrl OA for MS and CS coherence 
% flist_mod_ms = {16:18,15:17,16:18,16:18,15:17,15:17};
% flist_SA_ms = {12:14,12:14,15:17,15:17,14:16,14:16};
% flist_OA_ms = {12:14,15:17,15:17,16:18,15:17,15:17};

flist_mod_ms = {7:15};
flist_SA_ms = {7:15};
flist_OA_ms = {7:15};

% -- arrays to store the average values around the theta peak for MR, CR
mod_avg_mr = zeros(1,length(N_list));
err_mod_avg_mr = zeros(1,length(N_list));
ctrl_SA_avg_mr = zeros(1,length(N_list));
err_ctrl_SA_avg_mr = zeros(1,length(N_list));
ctrl_OA_avg_mr = zeros(1,length(N_list));
err_ctrl_OA_avg_mr = zeros(1,length(N_list));

% -- arrays to store the average values around the theta peak for MS, CS
mod_avg_ms = zeros(1,length(N_list));
err_mod_avg_ms = zeros(1,length(N_list));
ctrl_SA_avg_ms = zeros(1,length(N_list));
err_ctrl_SA_avg_ms = zeros(1,length(N_list));
ctrl_OA_avg_ms = zeros(1,length(N_list));
err_ctrl_OA_avg_ms = zeros(1,length(N_list));


for i=1:length(N_list)
    
    close all
    N = N_list(i); % --- max number of modulators    
    
    dir_both_monkeys = strcat(dir_main,sprintf('both_monkeys/%s/modulators_vs_controls/%s/%s',freq_band,recording_both,decoding));
    
    
    % %%%%%%%%% MODULATORS and CONTROLS for a given N %%%%%%
    load(strcat(dir_both_monkeys,sprintf('/modulators_N_%d.mat',N))); % structure mod
    load(strcat(dir_both_monkeys,sprintf('/controls_same_area_N_%d.mat',N))); % structure mod
    load(strcat(dir_both_monkeys,sprintf('/controls_other_areas_N_%d.mat',N))); % structure mod

    %%%%%%%%%%%%%%%%%%%
    % >>>>> MR / CR Coherence
    
    % average of MR at the pick and error bars 
    mod_avg_mr(i) = mean(modulators.mean_coh_mr(flist_mod_mr{1}));
    err_mod_avg_mr(i) = sqrt(sum(modulators.err_mr(flist_mod_mr{1}).^2))/length(flist_mod_mr{1});
    
    % average of CR same area at the pick and error bars 
    ctrl_SA_avg_mr(i) = mean(ctrl_SA.mean_coh_mr(flist_SA_mr{1}));
    err_ctrl_SA_avg_mr(i) = sqrt(sum(ctrl_SA.err_mr(flist_SA_mr{1}).^2))/length(flist_mod_mr{1});
    
    % average of CR same area at the pick and error bars
    ctrl_OA_avg_mr(i) = mean(ctrl_OA.mean_coh_mr(flist_OA_mr{1}));
    err_ctrl_OA_avg_mr(i) = sqrt(sum(ctrl_OA.err_mr(flist_OA_mr{1}).^2))/length(flist_mod_mr{1});
    
    %%%%%%%%%%%%%%%%%%%
    % >>>> MS / CS Coherence
    
    % average of MS at the pick and error bars
    mod_avg_ms(i) = mean(modulators.mean_coh_ms(flist_mod_ms{1}));
    err_mod_avg_ms(i) = sqrt(sum(modulators.err_ms(flist_mod_ms{1}).^2))/length(flist_mod_mr{1});
    
    % average of CS same area at the pick and error bars
    ctrl_SA_avg_ms(i) = mean(ctrl_SA.mean_coh_ms(flist_SA_ms{1}));
    err_ctrl_SA_avg_ms(i) = sqrt(sum(ctrl_SA.err_ms(flist_SA_ms{1}).^2))/length(flist_mod_mr{1});
    
    % average of CS same area at the pick and error bars
    ctrl_OA_avg_ms(i) = mean(ctrl_OA.mean_coh_ms(flist_OA_ms{1}));
    err_ctrl_OA_avg_ms(i) = sqrt(sum(ctrl_OA.err_ms(flist_OA_ms{1}).^2))/length(flist_mod_mr{1});
    
 
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- Figure: Coherence MR vs TOP modulators
fig = figure;
% errorbar(N_list,mod_avg_mr,err_mod_avg_mr,'color',[0, 51, 0]/255,'LineWidth',1); hold on 
% errorbar(N_list,ctrl_SA_avg_mr,err_ctrl_SA_avg_mr,'color',[26 198 1]/255,'LineWidth',1); hold on
% errorbar(N_list,ctrl_OA_avg_mr,err_ctrl_OA_avg_mr,'color',[102, 255, 217]/255,'LineWidth',1);

errorbar(N_list,mod_avg_mr,err_mod_avg_mr,'color',[28 199 139]/255,'LineWidth',1); hold on 
errorbar(N_list,ctrl_SA_avg_mr,err_ctrl_SA_avg_mr,'color',[50 250 93]/255,'LineWidth',1); hold on
errorbar(N_list,ctrl_OA_avg_mr,err_ctrl_OA_avg_mr,'color',[19 148 92]/255,'LineWidth',1);

% shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[28 199 139]/255 },'patchSaturation',0.5); hold on
% shadedErrorBar(f,ctrl_SA.mean_coh_mr,ctrl_SA.err_mr,'lineprops',{'color',[50 250 93]/255 },'patchSaturation',0.5); hold on
% shadedErrorBar(f,ctrl_OA.mean_coh_mr,ctrl_OA.err_mr,'lineprops',{'color',[19 148 92]/255 },'patchSaturation',0.5); 

% xlim([6,105])
xlim([85,110])
ylim([0.1 0.57])
% xticks([0 10 20 30 40 50 60 80 100])
grid on
set(gca,'FontSize',14)
% title('Coherence at the Theta peak vs TOP modulators','FontSize',11')
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Theta-Coherence','FontName','Arial','FontSize',15);
legend('Modulator - Receiver','Controls-SA - Receiver','Controls-OA - Receiver','FontSize',10)
set(gcf, 'Position',  [100, 600, 650, 450])

% fname = strcat(dir_both_monkeys,'/Coherence_MR_vs_TOP_modulators_both_monkeys.png');


fname = strcat(dir_fig,'/Coherence_MR_vs_TOP_modulators_both_monkeys_only100.pdf');
saveas(fig,fname);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- Figure: Coherence MS vs TOP modulators
fig = figure;
errorbar(N_list,mod_avg_ms,err_mod_avg_ms,'color',[0.4940, 0.1840, 0.5560],'LineWidth',1); hold on 
errorbar(N_list,ctrl_SA_avg_ms,err_ctrl_SA_avg_ms,'color',[255, 51, 153]/255,'LineWidth',1); hold on
errorbar(N_list,ctrl_OA_avg_ms,err_ctrl_OA_avg_ms,'color',[255, 128, 128]/255,'LineWidth',1);


xlim([85,110])
% xlim([5,105])
ylim([0.05 0.35])
% xticks([ 10 20 30 40 50 60 80 100])
grid on
% title('Coherence at the Theta peak vs TOP modulators','FontSize',11')
set(gca,'FontSize',14)
% title('Coherence at the Theta peak vs TOP modulators','FontSize',11')
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Theta-Coherence','FontName','Arial','FontSize',15);
legend('Modulator - Receiver','Controls-SA - Receiver','Controls-OA - Receiver','FontSize',10)

set(gcf, 'Position',  [100, 600, 650, 450])


% fname = strcat(dir_both_monkeys,'/Coherence_MS_vs_TOP_modulators_both_monkeys.png');
fname = strcat(dir_fig,'/Coherence_MS_vs_TOP_modulators_both_monkeys_only100.pdf');
saveas(fig,fname);













