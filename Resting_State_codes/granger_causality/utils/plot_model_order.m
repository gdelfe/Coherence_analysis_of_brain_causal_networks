
% plot_model_order.m
function plot_model_order(AIC, BIC, fs, sess_id, dir_model)
    fig = figure(1); clf;
    plot_tsdata([AIC BIC]', {'AIC','BIC'}, 1/fs);
    title(sprintf('Model order estimation\nSess = %d', sess_id));
    saveas(fig, fullfile(dir_model, sprintf('Model_order_sess_%d.jpg', sess_id)));
end