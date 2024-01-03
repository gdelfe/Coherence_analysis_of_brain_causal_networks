


clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/toolbox/SVG_export/src')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_fig = strcat(dir_main,'Figures_paper');

freq_band = 'theta_band';
recording_both = 'last_rec-rec001_002';
decoding = 'AUC';
N = 100;
fk = 200;

dir_both_monkeys = strcat(dir_main,sprintf('both_monkeys/%s/modulators_vs_controls/%s/%s',freq_band,recording_both,decoding));

load(strcat(dir_both_monkeys,sprintf('/modulators_N_%d.mat',N))); % structure mod
load(strcat(dir_both_monkeys,sprintf('/controls_same_area_N_%d.mat',N)));
load(strcat(dir_both_monkeys,sprintf('/controls_other_areas_N_%d.mat',N)));
    
f = linspace(1,fk,size(modulators.mean_coh_ms,2)); % frequency values (range)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES MODULATORS vs CONTROLS %%%%%%%%%%%%


set(0,'DefaultFigureVisible','on')


% --- ELECTRODE-RECEIVER coherence   -------%


% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all


% shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,ctrl_SA.mean_coh_mr,ctrl_SA.err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,ctrl_OA.mean_coh_mr,ctrl_OA.err_mr,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); 


% shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[79, 250, 94]/255},'patchSaturation',0.6); hold on
shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[28 199 139]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_mr,ctrl_SA.err_mr,'lineprops',{'color',[50 250 93]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_mr,ctrl_OA.err_mr,'lineprops',{'color',[19 148 92]/255 },'patchSaturation',0.5); 
% shadedErrorBar(f,ctrl_OA.mean_coh_mr,ctrl_OA.err_mr,'lineprops',{'color',[64 122 97]/255 },'patchSaturation',0.5); 


% title(sprintf('Both animals: Abs MR coherence, %s - RS',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
legend({'Modulator - Receiver','Control-SA - Receiver','Control-OA - Receiver'},'FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
ylim([0 0.36])
grid on

% fname = strcat(dir_fig,sprintf('/coherency_MR_mod_vs_ctrl_both_monkeys_N_%d%%_%s',N,namefig));
% saveas(fig,fname)
fname = strcat(dir_fig,sprintf('/coherency_MR_mod_vs_ctrl_both_monkeys_N_%d.pdf',N));
saveas(fig,fname)


% --- ELECTRODE-SENDER coherence   -------%


fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_ms,ctrl_SA.err_ms,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_ms,ctrl_OA.err_ms,'lineprops',{'color',[255, 128, 128]/255},'patchSaturation',0.4); hold on

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
legend('Modulator - Receiver','Control-SA - Receiver','Control-OA - Receiver','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
ylim([0 0.25])
grid on

% fname = strcat(dir_both_monkeys,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namefig));
% saveas(fig,fname)
fname = strcat(dir_fig,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d.pdf',N));
saveas(fig,fname)



% --- ELECTRODE-RECEIVER coherence   -------%
% --- ELECTRODE-SENDER coherence   -------%


fig = figure;
hold all

shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[28 199 139]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
legend('Modulator - Receiver','Modulator - Sender','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
ylim([0 0.36])
grid on

% fname = strcat(dir_both_monkeys,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namefig));
% saveas(fig,fname)
fname = strcat(dir_fig,sprintf('/coherency_MR_vs_MS_mod_both_monkeys_N_%d.pdf',N));
saveas(fig,fname)


