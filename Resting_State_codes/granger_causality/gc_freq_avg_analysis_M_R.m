%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  @ Gino Del Ferraro, Oct 2024, NYU, Pesaran Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%


addpath('T:/People/Gino/Coherence_modulator_analysis/Gino_codes');
% Main directory path
dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
% dir_main = 'T:\People\Gino\Coherence_modulator_analysis\Shaoyu_data';


% ---- MAVERICK
monkey = 'Maverick';
dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_ctrl_SA_struct.mat"));
gc_mav = gc_data.gc;


% ---- ARCHIE
monkey = 'Archie';
dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_ctrl_SA_struct.mat"));
gc_arc = gc_data.gc;

XY_rate = 'RS_rate';
YX_rate = 'SR_rate';
cnt_tot = 'cnt_r_tot';

tot_rate_XY = (gc_mav.(XY_rate)*gc_mav.(cnt_tot) + gc_arc.(XY_rate)*gc_arc.(cnt_tot))/(gc_mav.(cnt_tot) + gc_arc.(cnt_tot))
tot_rate_YX = (gc_mav.(YX_rate)*gc_mav.(cnt_tot) + gc_arc.(YX_rate)*gc_arc.(cnt_tot))/(gc_mav.(cnt_tot) + gc_arc.(cnt_tot))


[gc_avg_SR_ctrl_SA] = average_gc_XY(gc_mav,gc_arc,'RS_all', 'SR_all');


% --- BOTH MONKEYS

figure;
shadedErrorBar(gc_avg.freq,gc_avg.RM_mean,gc_avg.RM_sem,'lineprops',{'color','blue'},'patchSaturation',0.4);
shadedErrorBar(gc_avg.freq,gc_avg.MR_mean,gc_avg.MR_sem,'lineprops',{'color','green'},'patchSaturation',0.4);

grid on
title('Granger - Both Monkeys','FontSize',11);
xlabel('freq (Hz)');
ylabel('Granger');
legend('RM','MR','FontSize',10)
grid on
xlim([0,80])
