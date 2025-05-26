% Main script: Granger Causality Analysis

clear all; close all;
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

addpath('T:/People/Gino/Coherence_modulator_analysis/Gino_codes');
addpath('./utils');

dir_main = '/vol/brains/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
directory_model = 'Resting_state/granger_gino/model_order/SR';
directory_graphs = 'Resting_state/granger/granger_graphs/SR';
dir_output = 'Resting_state/granger/granger_results/SR';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat';
recording = 'last_recording';
freq_band = 'theta_band';
monkey = 'Maverick';
output_filename = "gc_SR_struct.mat";
U =  'S';
V = 'R';
Lag = 50;
fs = 1000;
fres = 500;

% parameters for VAR
regmode = 'OLS'; icregmode = 'LWR'; morder = 'AIC'; momax = 50;
acmaxlags = 2000; tstat = 'F'; alpha = 0.05; mhtc = 'FDR';

% File paths and session info

session_info_file = fullfile(dir_main, monkey, 'Resting_state', freq_band, 'Sessions_with_modulator_info_movie.txt');
session_info_file = strrep(session_info_file, '\\', '/');
fid = fopen(session_info_file);
if fid == -1
    error('Failed to open file: %s', session_info_file);
end
sess_info = textscan(fid, '%d%s%s');
fclose(fid);

gc = {}; sess_cnt = 0; sr_cnt = 1; rs_cnt = 1;

for i = 1:size(sess_info{1},1)

    close all
    Sess = sess_info{1}(i);
    dir_RS_Theta = fullfile(dir_main, monkey, 'Resting_state', freq_band);
    dir_Modulators = fullfile(dir_RS_Theta, sprintf('Sess_%d/Modulators', Sess));
    load(strcat(dir_Modulators,name_struct_input));

    dir_Mod_recording = fullfile(dir_Modulators, recording);
    if ~exist(dir_Mod_recording, 'dir'), mkdir(dir_Mod_recording); end

    sess_cnt = sess_cnt + 1;
    lfp_S = sess_data_lfp.lfp_S;
    lfp_R = sess_data_lfp.lfp_R;
    outliers_SR = unique([sess_data_lfp.outliers_S, sess_data_lfp.outliers_R]);
    lfp_S(outliers_SR,:) = [];
    lfp_R(outliers_SR,:) = [];

    lfp_S = detrend(lfp_S')';  % Remove linear trends
    lfp_R = detrend(lfp_R')';

    lfp_S = resample(lfp_S, 1, 4);  % Downsample
    lfp_R = resample(lfp_R, 1, 4);
    fs = 250;

    X = permute(cat(3, lfp_S, lfp_R), [3,2,1]);
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
    var_acinfo(info,true);

    [F, pval, sig] = compute_time_domain_gc(G, X, morder, nobs, ntrials, tstat, alpha, mhtc);
    sig(isnan(sig)) = 0;

    dir_granger = fullfile(dir_main, monkey, directory_graphs);
    if ~exist(dir_granger, 'dir'), mkdir(dir_granger); end

    plot_time_domain_gc(F, U, V, pval, sig, i, alpha, dir_granger);

    cd = mean(F(~isnan(F)));
    fprintf('\ncausal density = %f\n',cd);

    [GC, freq_axis] = autocov_to_spwcgc(G,fres, [], fs);
    if sig(1,2)==1
        gc(rs_cnt).RS = GC;
        gc(rs_cnt).freq = freq_axis;
        gc(rs_cnt).pval = pval(1,2);
        gc(rs_cnt).sess = Sess;
        gc(rs_cnt).sess_i = i;
        rs_cnt = rs_cnt + 1;
    end
    if sig(2,1)==1
        gc(sr_cnt).SR = GC;
        gc(sr_cnt).freq = freq_axis;
        gc(sr_cnt).pval = pval(2,1);
        gc(sr_cnt).sess = Sess;
        gc(sr_cnt).sess_i = i;
        sr_cnt = sr_cnt + 1;
    end

    % Save data needed for permutation test
    session_data(i).GC = GC;
    session_data(i).X = X;
    session_data(i).freq = freq_axis;
    session_data(i).sess_id = Sess;
    session_data(i).morder = morder;
    session_data(i).ntrials = ntrials;


    assert(~isbad(GC, false),'spectral GC calculation failed');
    plot_frequency_domain_gc(GC, fs, i, dir_granger, U, V);

    Fint = smvgc_to_mvgc(GC);
    amax = maxabs(F+Fint)/2; if amax < 1e-5, amax = 1; end
    mre = maxabs(F-Fint)/amax;
    if mre < 1e-5
        fprintf('OK (max relative error ~ %.0e)\n',mre);
    else
        fprintf(2,'WARNING: high max relative error ~ %.0e\n',mre);
    end
end

RS_rate = (rs_cnt-1)/sess_cnt;
SR_rate = (sr_cnt-1)/sess_cnt;
gc(1).RS_rate = RS_rate;
gc(1).SR_rate = SR_rate;
gc(1).rs_cnt_tot = rs_cnt-1;
gc(1).sr_cnt_tot = sr_cnt-1;
gc(1).sess_tot = sess_cnt;

dir_name = fullfile(monkey, dir_output);
file_path = make_dir_get_file_path(dir_main, dir_name, output_filename);
save(file_path,'gc');

file_path = make_dir_get_file_path(dir_main, dir_name, 'gc_session_data.mat');
save(file_path, 'session_data');

