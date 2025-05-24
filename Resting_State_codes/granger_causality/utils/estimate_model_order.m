
% estimate_model_order.m
function [morder, AIC, BIC] = estimate_model_order(X, momax, icregmode, fs, sess_id, dir_model)
    [AIC, BIC, moAIC, moBIC] = tsdata_to_infocrit(X, momax, icregmode);
    fig = figure(1); clf;
    plot_tsdata([AIC BIC]', {'AIC', 'BIC'}, 1/fs);
    title(sprintf('Model order estimation\nSess = %d', sess_id));
    saveas(fig, fullfile(dir_model, sprintf('Model_order_sess_%d.jpg', sess_id)));
    fprintf('\nbest model order (AIC) = %d\n', moAIC);
    fprintf('best model order (BIC) = %d\n', moBIC);
    morder = moAIC;
end