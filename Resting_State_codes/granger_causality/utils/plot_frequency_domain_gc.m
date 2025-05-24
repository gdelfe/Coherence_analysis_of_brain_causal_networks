
% % plot_frequency_domain_gc.m
% function plot_frequency_domain_gc(GC, fs, sess_id, dir_granger, U, V)
%     fig = figure(3); clf;
%     set(gcf, 'Position', [300, 300, 1000, 1000]);
%     sgtitlex(sprintf('Pairwise Granger causality - frequency domain\nSess = %d', sess_id));
%     plot_spw_Gino(GC, fs, [0,50], U, V);
%     saveas(fig, fullfile(dir_granger, sprintf('Sess_%d_granger_frequency.jpg', sess_id)));
% end

function plot_frequency_domain_gc(GC, fs, sess_id, dir_granger, U, V, threshold)
    % plot_frequency_domain_gc
    % Inputs:
    % - GC: spectral Granger causality (2x2xF)
    % - fs: sampling frequency
    % - sess_id: session index
    % - dir_granger: output directory
    % - U, V: channel names ('S', 'R')
    % - threshold: optional threshold line (same size as frequency vector or scalar)

    if nargin < 7
        threshold = [];
    end

    % Frequency axis
    nFreqs = size(GC, 3);
    freq_axis = linspace(0, fs/2, nFreqs);

    % Convert to dB
    GC_dB = 10 * log10(GC);
    GC_dB(isinf(GC_dB)) = NaN;  % handle log(0) cases

    % Plot
    fig = figure(3); clf;
    set(gcf, 'Position', [300, 300, 1200, 600]);

    subplot(1,2,1)
    plot(freq_axis, squeeze(GC_dB(1,2,:)), 'b', 'LineWidth', 2); hold on;
    if ~isempty(threshold)
        if numel(threshold) == 1
            yline(10*log10(threshold), 'r--', 'Threshold');
        else
            plot(freq_axis, 10*log10(threshold), 'r--');
        end
    end
    xlabel('Frequency (Hz)'); ylabel('GC S → R (dB)');
    title(sprintf('%s → %s (Sess %d)', U, V, sess_id)); xlim([0 100]); grid on;

    subplot(1,2,2)
    plot(freq_axis, squeeze(GC_dB(2,1,:)), 'g', 'LineWidth', 2); hold on;
    if ~isempty(threshold)
        if numel(threshold) == 1
            yline(10*log10(threshold), 'r--', 'Threshold');
        else
            plot(freq_axis, 10*log10(threshold), 'r--');
        end
    end
    xlabel('Frequency (Hz)'); ylabel(sprintf('GC %s → %s (dB)',V,U));
    title(sprintf('%s → %s (Sess %d),',V, U, sess_id)); xlim([0 100]); grid on;

    sgtitle(sprintf('Frequency-domain Granger Causality (Session %d)', sess_id));
    saveas(fig, fullfile(dir_granger, sprintf('Sess_%d_granger_frequency.jpg', sess_id)));
end
