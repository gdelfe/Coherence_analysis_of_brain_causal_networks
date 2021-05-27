function [data] = mean_coh_across_monkeys(mav,arc)

    % --- mean coherences
data.mean_coh_ms = (arc.modulators.mean_coh_ms + mav.modulators.mean_coh_ms)/2;  % modulator - sender
data.mean_coh_mr = (arc.modulators.mean_coh_mr + mav.modulators.mean_coh_mr)/2;  % modulator - receiver
data.mean_coh_sr = (arc.modulators.mean_coh_sr + mav.modulators.mean_coh_sr)/2;  % sender - receiver 

% --- Error bars
data.err_ms = (arc.modulators.err_ms + mav.modulators.err_ms)/sqrt(2);
data.err_mr = (arc.modulators.err_mr + mav.modulators.err_mr)/sqrt(2);
data.err_sr = (arc.modulators.err_sr + mav.modulators.err_sr)/sqrt(2);


end