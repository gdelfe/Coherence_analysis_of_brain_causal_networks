function [data] = mean_coh_and_spec_across_monkeys(struct_mod_mav,struct_stim_mav,struct_mod_arc,struct_stim_arc)

% MAVERICK 
% -- structures to matrices
mod_mav = cell2mat(struct2cell(struct_mod_mav)); % transform struct to mat for modulators
stim_mav = cell2mat(struct2cell(struct_stim_mav)); % transform struct to mat for sender-receiver 

% -- assign fields to matrices
coh_ms_M = sq(mod_mav(1,:,:))'; % 1st field, c_ms
coh_mr_M = sq(mod_mav(2,:,:))'; %  2nd field, c_mr
spec_m_M = sq(mod_mav(3,:,:))'; %  3rd field, spec_m

coh_sr_M = sq(stim_mav(1,:,:))'; % 1st field, c_sr
spec_s_M = sq(stim_mav(2,:,:))'; %  2nd field, spec_s
spec_r_M = sq(stim_mav(3,:,:))'; %  3rd field, spec_r


% ARCHIE 
% -- structures to matrices
mod_arc = cell2mat(struct2cell(struct_mod_arc)); % transform struct to mat for modulators
stim_arc = cell2mat(struct2cell(struct_stim_arc)); % transform struct to mat for sender-receiver 

% -- assign fields to matrices
coh_ms_A = sq(mod_arc(1,:,:))'; % 1st field, c_ms
coh_mr_A = sq(mod_arc(2,:,:))'; %  2nd field, c_mr
spec_m_A = sq(mod_arc(3,:,:))'; %  3rd field, spec_m

coh_sr_A = sq(stim_arc(1,:,:))'; % 1st field, c_sr
spec_s_A = sq(stim_arc(2,:,:))'; %  2nd field, spec_s
spec_r_A = sq(stim_arc(3,:,:))'; %  3rd field, spec_r


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS COHERENCES and SPECTRUM vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = size(coh_mr_M,1) + size(coh_mr_A,1); % -- number of electrodes
data.num_el = M;
S = size(coh_sr_M,1) + size(coh_sr_A,1); % -- number of senders (receivers)
data.num_send = S;

% concatenate coherences data
coh_ms = [coh_ms_A; coh_ms_M];
coh_mr = [coh_mr_A; coh_mr_M];
coh_sr = [coh_sr_A; coh_sr_M];

% concatenate spectrum data
spec_s = [spec_s_A; spec_s_M];
spec_r = [spec_r_A; spec_r_M];
spec_m = [spec_m_A; spec_m_M];


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
data.mean_spec_m = mean(spec_m);  % mean spectrum modulator
data.mean_spec_s = mean(spec_s);  % mean spectrum sender
data.mean_spec_r = mean(spec_r);  % mean spectrum receiver


% --- std spectrum
data.std_spec_m = std(spec_m);  
data.std_spec_s = std(spec_s);  
data.std_spec_r = std(spec_r);  


% --- Error bars spectrum
data.err_S_m = data.std_spec_m/sqrt(M); % -- all
data.err_S_s = data.std_spec_s/sqrt(S);
data.err_S_r = data.std_spec_r/sqrt(S);



end