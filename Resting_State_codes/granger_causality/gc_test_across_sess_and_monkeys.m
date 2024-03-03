
clear all; close all;


dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
Lag = 10; % maxLag for the computation of GC test 

% Monkey Maverick -----------------
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_GC_results = strcat(dir_RS_Theta,sprintf('/GC_results/'));

load(strcat(dir_GC_results,sprintf('GC_rate_%d.mat',Lag)))
gc_maverick = gc_struct;

% Monkey Archie   -----------------
monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_GC_results = strcat(dir_RS_Theta,sprintf('/GC_results/'));

load(strcat(dir_GC_results,sprintf('GC_rate_%d.mat',Lag)))
gc_archie = gc_struct;

SR = 0; RS = 0;
SM = 0; MS = 0; RM = 0; MR = 0;
cnt_m = 0; % total modulator count 
% Maverick
for i=1:size(gc_maverick,2)

    SR = SR + gc_maverick(i).SR; % SR granger rate
    RS = RS + gc_maverick(i).RS; % RS granger rate 

    if size(gc_maverick(i).mod,2) > 0 % if there is at least one modulator for this session
        for ch = size(gc_maverick(i).mod,2)

            SM = SM + gc_maverick(i).mod(ch).SM; % SM granger rate
            MS = MS + gc_maverick(i).mod(ch).MS;
            RM = RM + gc_maverick(i).mod(ch).RM; % RM granger rate
            MR = MR + gc_maverick(i).mod(ch).MR;

            cnt_m = cnt_m + 1;
        end
    end
end

% Archie
for i=1:size(gc_archie,2)

    SR = SR + gc_archie(i).SR;
    RS = RS + gc_archie(i).RS;

    if size(gc_archie(i).mod,2) > 0 % if there is at least one modulator for
        for ch = size(gc_archie(i).mod,2)

            SM = SM + gc_archie(i).mod(ch).SM;
            MS = MS + gc_archie(i).mod(ch).MS;
            RM = RM + gc_archie(i).mod(ch).RM;
            MR = MR + gc_archie(i).mod(ch).MR;

            cnt_m = cnt_m + 1;
        end
    end
end

% number of total sessions, i.e. total number of SR pairs
cnt_sr = size(gc_maverick,2) + size(gc_archie,2);

% Store results across monkeys and sessions into a structure. Divide to
% normalize the results

gc_avg.SR = SR/cnt_sr;
gc_avg.RS = RS/cnt_sr;

gc_avg.SM = SM/cnt_m;
gc_avg.MS = MS/cnt_m;
gc_avg.RM = RM/cnt_m;
gc_avg.MR = MR/cnt_m;

dir_out = strcat(dir_main,sprintf('both_monkeys/theta_band/gc_results/'));
if ~exist(dir_out, 'dir')
    mkdir(dir_out)
end
save(strcat(dir_out,sprintf('GC_rate_average_%d.mat',Lag)),'gc_avg');










