% Granger Permutation Test Script
% --------------------------------------------------
% Performs:
% 1. Per-session permutation-based GC significance test (trial shuffling)
% 2. Group-level null GC distribution and threshold (averaged across sessions)

% Load session GC data saved from main_script_unconditional_granger_SR.m
dir_main = '/vol/brains/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_output = 'Resting_state/granger/granger_results/SR';
monkey = 'Archie';

% Load session GC data saved from main_script_unconditional_granger_SR.m
load(fullfile(dir_main, monkey, dir_output,'gc_session_data.mat'));


nPerm = 200;                 
alpha_perm = 0.05;           % Significance level
nSessions = length(session_data);

% Preallocate cell arrays to handle variable frequency lengths
all_GC_SR = cell(nSessions, 1);
all_GC_RS = cell(nSessions, 1);
perm_GC_SR = cell(nSessions, 1);
perm_GC_RS = cell(nSessions, 1);
GC_thresh_SR = cell(nSessions, 1);
GC_thresh_RS = cell(nSessions, 1);

for i = 1:nSessions
    fprintf('Processing session %d/%d (ID: %d)\n', i, nSessions, session_data(i).sess_id);

    % Frequency axis and resolution for this session
    freq_axis_real = session_data(i).freq;
    fres_real = length(freq_axis_real);

    % Real GC
    GC = session_data(i).GC;
    fres_gc = size(GC, 3);
    min_real_len = min(fres_real, fres_gc);
    all_GC_SR{i} = squeeze(GC(1,2,1:min_real_len))';
    all_GC_RS{i} = squeeze(GC(2,1,1:min_real_len))';

    % Trial shuffling
    X = session_data(i).X;
    morder = session_data(i).morder;
    ntrials = session_data(i).ntrials;

    perm_GC_SR{i} = zeros(nPerm, min_real_len);
    perm_GC_RS{i} = zeros(nPerm, min_real_len);

    for p = 1:nPerm
        if mod(p, 50) == 0
            fprintf('  Permutation %d/%d\n', p, nPerm);
        end

        perm_idx = randperm(ntrials);
        X_perm = X;
        X_perm(2,:,:) = X(2,:,perm_idx);  % shuffle R trials only

        % Suppress internal output
        evalc("[A_p, SIG_p] = estimate_var_model(X_perm, morder, 'OLS');");
        evalc("[G_p, ~] = compute_autocovariance(A_p, SIG_p, 2000);");
        GC_perm = autocov_to_spwcgc(G_p, [], [], 250);  % fs = 250 Hz

        fres_perm = size(GC_perm, 3);
        min_len = min(min_real_len, fres_perm);

        gc_sr = squeeze(GC_perm(1,2,1:min_len));
        gc_rs = squeeze(GC_perm(2,1,1:min_len));

        % Ensure correct shape (1 x min_len)
        if size(gc_sr, 1) > size(gc_sr, 2)
            gc_sr = gc_sr';
        end
        if size(gc_rs, 1) > size(gc_rs, 2)
            gc_rs = gc_rs';
        end

        perm_GC_SR{i}(p,1:min_len) = gc_sr;
        perm_GC_RS{i}(p,1:min_len) = gc_rs;
    end

    % Compute per-session threshold (95th percentile)
    GC_thresh_SR{i} = prctile(perm_GC_SR{i}, 100*(1 - alpha_perm), 1);
    GC_thresh_RS{i} = prctile(perm_GC_RS{i}, 100*(1 - alpha_perm), 1);
end

% Align for group-level analysis
common_len = min(cellfun(@length, all_GC_SR));
GC_mat_SR = cell2mat(cellfun(@(x) x(1:common_len), all_GC_SR, 'UniformOutput', false));
GC_mat_RS = cell2mat(cellfun(@(x) x(1:common_len), all_GC_RS, 'UniformOutput', false));

GC_mean_SR = mean(GC_mat_SR, 1);
GC_mean_RS = mean(GC_mat_RS, 1);

% Compute group-level permutation thresholds
perm_mat_SR = cell2mat(cellfun(@(x) x(:,1:common_len), perm_GC_SR, 'UniformOutput', false));
perm_mat_RS = cell2mat(cellfun(@(x) x(:,1:common_len), perm_GC_RS, 'UniformOutput', false));
GC_mean_thresh_SR = prctile(perm_mat_SR, 100*(1 - alpha_perm), 1);
GC_mean_thresh_RS = prctile(perm_mat_RS, 100*(1 - alpha_perm), 1);

dir_name = fullfile(monkey, dir_output);
file_path = make_dir_get_file_path(dir_main, dir_name, 'gc_permutation_results.mat');

save(file_path, 'all_GC_SR', 'all_GC_RS', 'GC_thresh_SR', 'GC_thresh_RS', ...
     'GC_mean_SR', 'GC_mean_RS', 'GC_mean_thresh_SR', 'GC_mean_thresh_RS', 'freq_axis');

fprintf('Permutation test completed and results saved.\n');
