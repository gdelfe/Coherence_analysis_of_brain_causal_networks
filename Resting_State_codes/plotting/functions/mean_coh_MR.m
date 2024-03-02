function [data] = mean_coh_MR(mod)


% -- structures to matrices
mod_mat = cell2mat(struct2cell(mod)); % transform struct to mat for modulators

% -- assign fields to matrices
coh_ms = sq(mod_mat(1,:,:))'; % 1st field, c_ms
coh_mr = sq(mod_mat(2,:,:))'; %  2nd field, c_mr
spec_m = sq(mod_mat(3,:,:))'; %  3rd field, spec_m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS COHERENCES and SPECTRUM vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = size(coh_mr,1); % -- number of electrodes
data.num_el = M;

% ----------- COHERENCE

% --- mean coherences
data.mean_coh_ms = mean(abs(coh_ms));  % modulator - sender
data.mean_coh_mr = mean(abs(coh_mr));  % modulator - receiver

% --- std coherences
data.std_coh_ms = std(abs(coh_ms));  % modulator - sender
data.std_coh_mr = std(abs(coh_mr)); % modulator - receiver

% --- Error bars
data.err_ms = data.std_coh_ms/sqrt(M);
data.err_mr = data.std_coh_mr/sqrt(M);


% ----------- SPECTRUM

% --- mean spectrum
data.mean_spec_m = mean(spec_m);  % mean spectrum modulator

% --- std spectrum
data.std_spec_m = std(spec_m);    

% --- Error bars spectrum
data.err_S_m = data.std_spec_m/sqrt(M); % -- all


end