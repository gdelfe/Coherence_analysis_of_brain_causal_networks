function [data] = mean_coh_across_monkeys_SA(mav,arc)

    % --- mean coherences
data.mean_coh_ms = (arc.ctrl_SA.mean_coh_ms + mav.ctrl_SA.mean_coh_ms)/2;  % modulator - sender
data.mean_coh_mr = (arc.ctrl_SA.mean_coh_mr + mav.ctrl_SA.mean_coh_mr)/2;  % modulator - receiver
data.mean_coh_sr = (arc.ctrl_SA.mean_coh_sr + mav.ctrl_SA.mean_coh_sr)/2;  % sender - receiver 

% --- Error bars
data.err_ms = (arc.ctrl_SA.err_ms + mav.ctrl_SA.err_ms)/sqrt(2);
data.err_mr = (arc.ctrl_SA.err_mr + mav.ctrl_SA.err_mr)/sqrt(2);
data.err_sr = (arc.ctrl_SA.err_sr + mav.ctrl_SA.err_sr)/sqrt(2);


end