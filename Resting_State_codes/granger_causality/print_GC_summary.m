% Print GC Summary Statistics from Saved Results
% --------------------------------------------------
% Loads the results from the permutation-based GC analysis and prints:
% - Number of sessions
% - Frequency resolution
% - Number and percentage of significant frequencies at the group level
% - Peak GC values and corresponding frequencies
% - Full lists of significant frequency bands

set(0,'DefaultFigureVisible','on')

% Define paths
monkey = 'Maverick';
dir_main = '/vol/brains/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
dir_output = 'Resting_state/granger/granger_results/SR';
dir_name = fullfile(monkey, dir_output);
file_path = fullfile(dir_main, dir_name, 'gc_permutation_results.mat');

% Load output
load(file_path);

% Print basic info
disp('--- Summary of Granger Causality Results ---');

nSessions = length(all_GC_SR);
fprintf('Number of sessions: %d\n', nSessions);

% Frequency resolution
fprintf('Common frequency resolution (used for group-level): %d points\n', length(GC_mean_SR));

% Frequency range (assumes linear spacing)
dfreq = 250 / (2 * length(GC_mean_SR));
freq_axis = linspace(0, 250/2, length(GC_mean_SR));

% Significant frequencies at group level
sig_SR = GC_mean_SR > GC_mean_thresh_SR;
sig_RS = GC_mean_RS > GC_mean_thresh_RS;

fprintf('\nSignificant frequencies (S->R): %d/%d (%.1f%%)\n', sum(sig_SR), length(sig_SR), 100*mean(sig_SR));
fprintf('Significant frequencies (R->S): %d/%d (%.1f%%)\n', sum(sig_RS), length(sig_RS), 100*mean(sig_RS));

% Print frequency bands where GC is significant
fprintf('\nFrequencies where GC (S->R) is significant:\n');
disp(freq_axis(sig_SR));

fprintf('Frequencies where GC (R->S) is significant:\n');
disp(freq_axis(sig_RS));

% Print peak GC values and corresponding frequencies
[max_val_SR, idx_SR] = max(GC_mean_SR);
[max_val_RS, idx_RS] = max(GC_mean_RS);

fprintf('\nPeak GC (S->R): %.4f at %.2f Hz\n', max_val_SR, freq_axis(idx_SR));
fprintf('Peak GC (R->S): %.4f at %.2f Hz\n', max_val_RS, freq_axis(idx_RS));

% Optional: plot group-level GC curves and thresholds
figure;
plot(freq_axis, GC_mean_SR, 'b', 'LineWidth', 2); hold on;
plot(freq_axis, GC_mean_thresh_SR, 'r--', 'LineWidth', 1.5);
xlim([0, 100])
title('Group Mean GC S\rightarrowR'); xlabel('Frequency (Hz)'); ylabel('GC'); legend('Mean GC', 'Threshold'); grid on;

figure;
plot(freq_axis, GC_mean_RS, 'g', 'LineWidth', 2); hold on;
plot(freq_axis, GC_mean_thresh_RS, 'r--', 'LineWidth', 1.5);
xlim([0, 100])
title('Group Mean GC R\rightarrowS'); xlabel('Frequency (Hz)'); ylabel('GC'); legend('Mean GC', 'Threshold'); grid on;
