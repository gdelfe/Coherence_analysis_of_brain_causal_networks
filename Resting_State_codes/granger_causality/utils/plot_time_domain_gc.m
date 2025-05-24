% plot_time_domain_gc.m
function plot_time_domain_gc(F, U, V, pval, sig, sess_id, alpha, dir_granger)
    fig = figure(2); clf;
    set(gcf, 'Position', [100, 100, 1400, 400]);
    sgtitlex(sprintf('Pairwise Granger causality - time domain \nSession = %d', sess_id));
    subplot(1,4,1); plot_pw_Gino(F, U, V); title('Pairwise-conditional GC');
    subplot(1,4,2); plot_pw_Gino(pval,U, V); title('p-values');
    subplot(1,4,3); plot_pw_Gino(sig,U, V); title(['Significant at p = ' num2str(alpha)]);
    subplot(1,4,4); Granger_graph(sig, {U,V},'Granger interaction');
    saveas(fig, fullfile(dir_granger, sprintf('Sess_%d_granger_graph.jpg', sess_id)));
end
