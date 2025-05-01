

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
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_struct.mat"));
gc_mav = gc_data.gc;


% ---- ARCHIE
monkey = 'Archie';
dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_struct.mat"));
gc_arc = gc_data.gc;




% ---- MAVERICK
monkey = 'Maverick';
dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_ctrl_SA_struct.mat"));
gc_mav_SA = gc_data.gc;


% ---- ARCHIE
monkey = 'Archie';
dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_ctrl_SA_struct.mat"));
gc_arc_SA = gc_data.gc;


% ---- MAVERICK
monkey = 'Maverick';
dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_ctrl_OA_struct.mat"));
gc_mav_OA = gc_data.gc;


% ---- ARCHIE
monkey = 'Archie';
dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_ctrl_OA_struct.mat"));
gc_arc_OA = gc_data.gc;

XY_rate = 'RS_rate';
YX_rate = 'SR_rate';
cnt_tot = 'cnt_r_tot';
% 
% tot_rate_XY = (gc_mav.(XY_rate)*gc_mav.(cnt_tot) + gc_arc.(XY_rate)*gc_arc.(cnt_tot))/(gc_mav.(cnt_tot) + gc_arc.(cnt_tot))
% tot_rate_YX = (gc_mav.(YX_rate)*gc_mav.(cnt_tot) + gc_arc.(YX_rate)*gc_arc.(cnt_tot))/(gc_mav.(cnt_tot) + gc_arc.(cnt_tot))
% 

% [gc_avg] = average_gc_XY(gc_mav,gc_arc,'RS_all', 'SR_all');
gc_avg = average_gc_across_monkeys(gc_mav,gc_arc,'RS','SR');
[gc_avg_SA] = average_gc_XY(gc_mav_SA,gc_arc_SA,'RS_all', 'SR_all');
[gc_avg_OA] = average_gc_XY(gc_mav_OA,gc_arc_OA,'RS_all', 'SR_all');


% --- BOTH MONKEYS

figure;
% shadedErrorBar(gc_avg.freq,gc_avg.XY_mean,gc_avg.XY_sem,'lineprops',{'color','blue'},'patchSaturation',0.4);
shadedErrorBar(gc_avg_SA.freq,gc_avg.RS.mean,gc_avg.RS.sem,'lineprops',{'color','blue'},'patchSaturation',0.4)
shadedErrorBar(gc_avg_SA.freq,gc_avg_SA.XY_mean,gc_avg_SA.XY_sem,'lineprops',{'color','green'},'patchSaturation',0.4);
shadedErrorBar(gc_avg_OA.freq,gc_avg_OA.XY_mean,gc_avg_OA.XY_sem,'lineprops',{'color','magenta'},'patchSaturation',0.4);


grid on
title('Granger - Both Monkeys - SR ctrl','FontSize',11);
xlabel('freq (Hz)');
ylabel('Granger');
legend('RS','R ctrl SA-S','R ctrl OA-S','FontSize',10)
grid on
xlim([0,80])





figure;
% shadedErrorBar(gc_avg.freq,gc_avg.XY_mean,gc_avg.XY_sem,'lineprops',{'color','blue'},'patchSaturation',0.4);
shadedErrorBar(gc_avg_SA.freq,gc_avg.SR.mean,gc_avg.SR.sem,'lineprops',{'color','blue'},'patchSaturation',0.4)
shadedErrorBar(gc_avg_SA.freq,gc_avg_SA.YX_mean,gc_avg_SA.YX_sem,'lineprops',{'color','green'},'patchSaturation',0.4);
shadedErrorBar(gc_avg_OA.freq,gc_avg_OA.YX_mean,gc_avg_OA.YX_sem,'lineprops',{'color','magenta'},'patchSaturation',0.4);


grid on
title('Granger - Both Monkeys - SR ctrl','FontSize',11);
xlabel('freq (Hz)');
ylabel('Granger');
legend('SR','SR ctrl SA','SR ctrl OA','FontSize',10)
grid on
xlim([0,80])
