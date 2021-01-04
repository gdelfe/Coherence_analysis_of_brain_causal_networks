function [data] = mean_coh_and_spec_RS(mod,stim)


% -- structures to matrices
mod_mat = cell2mat(struct2cell(mod)); % transform struct to mat for modulators
stim_mat = cell2mat(struct2cell(stim)); % transform struct to mat for sender-receiver 

% -- assign fields to matrices
coh_ms = sq(mod_mat(1,:,:))'; % 1st field, c_ms
coh_mr = sq(mod_mat(2,:,:))'; %  2nd field, c_mr
spec_m = sq(mod_mat(3,:,:))'; %  3rd field, spec_m

coh_sr = sq(stim_mat(1,:,:))'; % 1st field, c_sr
spec_s = sq(stim_mat(2,:,:))'; %  2nd field, spec_s
spec_r = sq(stim_mat(3,:,:))'; %  3rd field, spec_r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS COHERENCES and SPECTRUM vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = size(coh_mr,1); % -- number of electrodes
data.num_el = M;
S = size(coh_sr,1); % -- number of senders (receivers)
data.num_send = S;

% ----------- COHERENCE

% --- mean coherences
data.mean_coh_ms = mean(abs(coh_ms));  % modulator - sender
data.mean_coh_mr = mean(abs(coh_mr));  % modulator - receiver
data.mean_coh_sr = mean(abs(coh_sr));  % sender - receiver 

% --- std coherences
data.std_coh_ms = std(abs(coh_ms));  % modulator - sender
data.std_coh_mr = std(abs(coh_mr)); % modulator - receiver
data.std_coh_sr = std(abs(coh_sr));  % modulator - receiver

% --- Error bars
data.err_ms = data.std_coh_ms/sqrt(M);
data.err_mr = data.std_coh_mr/sqrt(M);
data.err_sr = data.std_coh_sr/sqrt(S);

% ----------- SPECTRUM

% --- mean spectrum
data.mean_spec_m = mean(abs(spec_m));  % mean spectrum modulator
data.mean_spec_s = mean(abs(spec_s));  % mean spectrum sender
data.mean_spec_r = mean(abs(spec_r));  % mean spectrum receiver


% --- std spectrum
data.std_spec_m = std(abs(spec_m));  
data.std_spec_s = std(abs(spec_s));  
data.std_spec_r = std(abs(spec_r));  


% --- Error bars spectrum
data.err_S_m = data.std_spec_m/sqrt(M); % -- all
data.err_S_s = data.std_spec_s/sqrt(S);
data.err_S_r = data.std_spec_r/sqrt(S);



end