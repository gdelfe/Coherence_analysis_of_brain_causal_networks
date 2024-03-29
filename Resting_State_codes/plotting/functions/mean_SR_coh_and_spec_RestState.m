function [data] = mean_SR_coh_and_spec_RestState(stim)


% -- structures to matrices
stim_mat = cell2mat(struct2cell(stim)); % transform struct to mat for sender-receiver 

% -- assign fields to matrices
coh_sr = sq(stim_mat(1,:,:))'; % 1st field, c_sr
spec_s = sq(stim_mat(2,:,:))'; %  2nd field, spec_s
spec_r = sq(stim_mat(3,:,:))'; %  3rd field, spec_r


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS COHERENCES and SPECTRUM vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


S = size(coh_sr,1); % -- number of senders (receivers)
data.num_send = S;

% ----------- COHERENCE  sender - receiver 

% --- mean coherences
data.mean_coh_sr = mean(abs(coh_sr));  
% --- std coherences
data.std_coh_sr = std(abs(coh_sr)); 
% --- Error bars
data.err_sr = data.std_coh_sr/sqrt(S);

% ----------- SPECTRUM

% --- mean spectrum
data.mean_spec_s = mean(spec_s);  % mean spectrum sender
data.mean_spec_r = mean(spec_r);  % mean spectrum receiver

% --- std spectrum
data.std_spec_s = std(spec_s);  
data.std_spec_r = std(spec_r);  

% --- Error bars spectrum
data.err_S_s = data.std_spec_s/sqrt(S);
data.err_S_r = data.std_spec_r/sqrt(S);



end