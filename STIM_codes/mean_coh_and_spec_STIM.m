function [data] = mean_coh_and_spec_STIM(stim)


stim_mat = cell2mat(struct2cell(stim)); % transform struct to mat for sender-receiver


coh_mr = sq(stim_mat(1,:,:))'; % 1st field, c_mr
spec_m = sq(stim_mat(2,:,:))'; %  2nd field, spec_m
spec_r = sq(stim_mat(3,:,:))'; %  3rd field, spec_r

coh_mr_H = sq(stim_mat(4,:,:))'; % 1st field, c_mr Hits
spec_m_H = sq(stim_mat(5,:,:))'; %  2nd field, spec_m Hits
spec_r_H = sq(stim_mat(6,:,:))'; %  3rd field, spec_r Hits

coh_mr_M = sq(stim_mat(7,:,:))'; % 1st field, c_mr Hits
spec_m_M = sq(stim_mat(8,:,:))'; %  2nd field, spec_m Misses
spec_r_M = sq(stim_mat(9,:,:))'; %  3rd field, spec_r Misses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS COHERENCES and SPECTRUM vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = size(stim,2); % -- number of electrodes
data.num_el = M;
S = size(coh_sr,1); % -- number of senders (receivers)
data.num_send = S;

% ----------- COHERENCE

% --- mean coherences
data.mean_coh_mr = mean(abs(coh_mr));  % modulator - sender
data.mean_coh_mr_H = mean(abs(coh_mr_H));  % modulator - receiver
data.mean_coh_mr_M = mean(abs(coh_mr_M));  % sender - receiver

% --- std coherences
data.std_mr = std(abs(coh_mr));  % modulator - sender
data.std_mr_H = std(abs(coh_mr_H)); % modulator - receiver
data.std_mr_M = std(abs(coh_mr_M));  % modulator - receiver

% --- Error bars choerencies
data.err_coh_mr = data.std_mr/sqrt(M);
data.err_coh_mr_H = data.std_mr_H/sqrt(M);
data.err_coh_mr_M = data.std_mr_M/sqrt(M);

% ----------- SPECTRUM

% --- mean spectrum
data.mean_spec_m = mean(abs(spec_m));  % mean spectrum modulator
data.mean_spec_r = mean(abs(spec_r));  % mean spectrum receiver
data.mean_spec_m_H = mean(abs(spec_m_H));  % mean spectrum modulator Hits
data.mean_spec_r_H = mean(abs(spec_r_H));  % mean spectrum receiver Hits
data.mean_spec_m_M = mean(abs(spec_m_M));  % mean spectrum modulator Misses
data.mean_spec_r_M = mean(abs(spec_r_M));  % mean spectrum receiver Misses

% --- std spectrum
data.std_spec_m = std(abs(spec_m));  % mean spectrum modulator
data.std_spec_r = std(abs(spec_r));  % mean spectrum receiver
data.std_spec_m_H = std(abs(spec_m_H));  % mean spectrum modulator Hits
data.std_spec_r_H = std(abs(spec_r_H));  % mean spectrum receiver Hits
data.std_spec_m_M = std(abs(spec_m_M));  % mean spectrum modulator Misses
data.std_spec_r_M = std(abs(spec_r_M));  % mean spectrum receiver Misses

% --- Error bars spectrum
data.err_S_m = data.std_spec_m/sqrt(M); % -- all
data.err_S_r = data.std_spec_r/sqrt(M);
data.err_S_m_H = data.std_spec_m_H/sqrt(M); % -- Hits
data.err_S_r_H = data.std_spec_r_H/sqrt(M);
data.err_S_m_M = data.std_spec_m_M/sqrt(M); % -- misses
data.err_S_r_M = data.std_spec_r_M/sqrt(M);



end