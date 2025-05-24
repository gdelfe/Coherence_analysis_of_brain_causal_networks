
% plot_frequency_domain_gc.m
function plot_frequency_domain_gc(GC, fs, sess_id, dir_granger, U, V)
    fig = figure(3); clf;
    set(gcf, 'Position', [300, 300, 1000, 1000]);
    sgtitlex(sprintf('Pairwise Granger causality - frequency domain\nSess = %d', sess_id));
    plot_spw_Gino(GC, fs, [0,50], U, V);
    saveas(fig, fullfile(dir_granger, sprintf('Sess_%d_granger_frequency.jpg', sess_id)));
end