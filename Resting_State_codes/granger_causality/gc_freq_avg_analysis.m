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


[gc_SR_stat] = average_gc_XY_YX(gc_mav,'RS','SR');

figure;
shadedErrorBar(gc_SR_stat.freq, gc_SR_stat.RS.mean, gc_SR_stat.RS.sem,'lineprops',{'color','blue'},'patchSaturation',0.4);
shadedErrorBar(gc_SR_stat.freq, gc_SR_stat.SR.mean, gc_SR_stat.SR.sem,'lineprops',{'color','green'},'patchSaturation',0.4);

grid on
title('Granger - Maverick','FontSize',11);
xlabel('freq (Hz)');
ylabel('Granger');
legend('RS','SR','FontSize',10)
% set(gcf, 'Position',  [100, 600, 1000, 600])
grid on
xlim([0,300])


% ---- ARCHIE

monkey = 'Archie';
dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
gc_data = load(fullfile(dir_main,dir_name, "gc_SR_struct.mat"));
gc_arc = gc_data.gc;


[gc_SR_stat] = average_gc_XY_YX(gc_arc,'RS','SR')

figure;
shadedErrorBar(gc_SR_stat.freq,gc_SR_stat.RS.mean,gc_SR_stat.RS.sem,'lineprops',{'color','blue'},'patchSaturation',0.4);
shadedErrorBar(gc_SR_stat.freq, gc_SR_stat.SR.mean, gc_SR_stat.SR.sem,'lineprops',{'color','green'},'patchSaturation',0.4);

grid on
title('Granger - Archie','FontSize',11);
xlabel('freq (Hz)');
ylabel('Granger');
legend('RS','SR','FontSize',10)
grid on
xlim([0,300])



% --- BOTH MONKEYS

gc_avg = average_gc_across_monkeys(gc_mav,gc_arc,'RS','SR')

figure;
shadedErrorBar(gc_avg.freq,gc_avg.RS.mean,gc_avg.RS.sem,'lineprops',{'color','blue'},'patchSaturation',0.4);
shadedErrorBar(gc_avg.freq,gc_avg.SR.mean,gc_avg.SR.sem,'lineprops',{'color','green'},'patchSaturation',0.4);

grid on
title('Granger - Both Monkeys','FontSize',11);
xlabel('freq (Hz)');
ylabel('Granger');
legend('RS','SR','FontSize',10)
grid on
xlim([0,300])


