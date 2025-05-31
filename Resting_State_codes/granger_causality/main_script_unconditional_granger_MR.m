% Main script: Granger Causality Analysis (Modulator <-> Receiver)

clear all; close all;
set(0,'DefaultFigureVisible','off')
set(0,'DefaultLineLineWidth',2)

addpath('./utils');
addpath('/vol/brains/bd5/People/Gino/Coherence_modulator_analysis/Gino_codes_v2')

dir_main = '/vol/brains/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
directory_model = 'Resting_state/granger_gino/model_order/MR';
directory_graphs = 'Resting_state/granger/granger_graphs/MR';
dir_output = 'Resting_state/granger/granger_results/MR';

name_struct_input = '/session_data_lfp.mat';
recording = 'last_recording';
freq_band = 'theta_band';
monkey = 'Maverick';
output_filename = "gc_MR_struct.mat";

Lag = 50; fs = 250; fres = 1000;

% parameters for VAR
regmode = 'LWR'; icregmode = 'LWR'; morder = 'AIC'; momax = 50;
acmaxlags = 2000; tstat = 'F'; alpha = 0.05; mhtc = 'FDR';

% File paths and session info
session_info_file = fullfile(dir_main, monkey, 'Resting_state', freq_band, 'Sessions_with_modulator_info_movie.txt');
fid = fopen(session_info_file);
sess_info = textscan(fid, '%d%s%s'); fclose(fid);

gc = {}; sess_cnt = 0; cnt_m_tot = 0;

for i = 1:size(sess_info{1},1)
    close all
    Sess = sess_info{1}(i);
    dir_RS_Theta = fullfile(dir_main, monkey, 'Resting_state', freq_band);
    dir_Modulators = fullfile(dir_RS_Theta, sprintf('Sess_%d/Modulators', Sess));
    load(strcat(dir_Modulators,name_struct_input));

    mod_Ch = sess_data_lfp.mod_idx;

    cnt_m = 0; % local modulator count per session
    for ch = mod_Ch
        if ismember(sess_data_lfp.receiver_idx, ch)
            continue;  % Skip if receiver is also a modulator
        end

        cnt_m = cnt_m + 1;
        cnt_m_tot = cnt_m_tot + 1;
        lfp_M = squeeze(sess_data_lfp.lfp_E(ch,:,:));
        lfp_R = sess_data_lfp.lfp_R;

        outliers_MR = unique([sess_data_lfp.outliers_E(cnt_m).idx, sess_data_lfp.outliers_R]);
        lfp_M(outliers_MR,:) = [];
        lfp_R(outliers_MR,:) = [];

        % % Remove DC offset
        % lfp_M = lfp_M - mean(lfp_M, 2);
        % lfp_R = lfp_R - mean(lfp_R, 2);

        lfp_M = detrend(lfp_M')';  % Remove linear trends
        lfp_R = detrend(lfp_R')';

        [b, a] = butter(2, 1/(fs/2), 'high');  % 1 Hz cutoff
        lfp_M = filtfilt(b, a, lfp_M')';
        lfp_R = filtfilt(b, a, lfp_R')';

        X = permute(cat(3, lfp_M, lfp_R), [3,2,1]);
        ntrials = size(X,3);
        nobs = size(X,2);

        dir_model = fullfile(dir_main, monkey, directory_model);
        if ~exist(dir_model, 'dir'), mkdir(dir_model); end

        [morder_val, AIC, BIC] = estimate_model_order(X, momax, icregmode, fs, i, dir_model);
        if strcmpi(morder,'AIC')
            morder = morder_val;
        end

        [A, SIG] = estimate_var_model(X, morder, regmode);
        [G, info] = compute_autocovariance(A, SIG, acmaxlags);

        [F, pval, sig] = compute_time_domain_gc(G, X, morder, nobs, ntrials, tstat, alpha, mhtc);
        sig(isnan(sig)) = 0;

        dir_granger = fullfile(dir_main, monkey, directory_graphs);
        if ~exist(dir_granger, 'dir'), mkdir(dir_granger); end

        plot_time_domain_gc(F, 'M', 'R', pval, sig, i, alpha, dir_granger);

        [GC, freq_axis] = autocov_to_spwcgc(G, fres, [], fs);
        assert(~isbad(GC, false),'spectral GC calculation failed');

        plot_frequency_domain_gc(GC, fs, i, dir_granger, 'M', 'R');

        %  Save GC if significant
        if sig(1,2)==1
            idx = length(gc) + 1;
            gc{idx}.MR = GC;
            gc{idx}.freq = freq_axis;
            gc{idx}.pval = pval(1,2);
            gc{idx}.sess = Sess;
            gc{idx}.sess_i = i;
            gc{idx}.mod_ch = ch;
            gc{idx}.morder = morder;
        end
        if sig(2,1)==1
            idx = length(gc) + 1;
            gc{idx}.RM = GC;
            gc{idx}.freq = freq_axis;
            gc{idx}.pval = pval(2,1);
            gc{idx}.sess = Sess;
            gc{idx}.sess_i = i;
            gc{idx}.mod_ch = ch;
            gc{idx}.morder = morder;
        end

        % Check numeric accuracy
        Fint = smvgc_to_mvgc(GC);
        amax = maxabs(F+Fint)/2; if amax < 1e-5, amax = 1; end
        mre = maxabs(F-Fint)/amax;
        if mre > 1e-5
            fprintf(2,'WARNING: high max relative error ~ %.0e\n',mre);
        end
    end

    sess_cnt = sess_cnt + 1;
end

gc{1}.sess_tot = sess_cnt;
gc{1}.cnt_modulators = cnt_m_tot;

MR_count = sum(cellfun(@(g) isfield(g, 'MR'), gc(2:end)));
RM_count = sum(cellfun(@(g) isfield(g, 'RM'), gc(2:end)));
total_modulators = cnt_m_tot;

gc{1}.MR_rate = MR_count / total_modulators;
gc{1}.RM_rate = RM_count / total_modulators;
gc{1}.MR_count = MR_count;
gc{1}.RM_count = RM_count;

dir_name = fullfile(monkey, dir_output);
file_path = make_dir_get_file_path(dir_main, dir_name, output_filename);
save(file_path,'gc');