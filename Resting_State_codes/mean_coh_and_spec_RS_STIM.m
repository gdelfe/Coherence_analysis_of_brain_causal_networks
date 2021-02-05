function [data] = mean_coh_and_spec_RS_STIM(resting,stim)


% -- structures to matrices
RS_mat = cell2mat(struct2cell(resting)); % transform struct to mat for modulators
stim_new = rmfield(stim,'n_trials_RS_STIM');
stim_mat = cell2mat(struct2cell(stim_new)); % transform struct to mat for sender-receiver 

% -- assign fields to matrices
RS_coh_mr = sq(RS_mat(1,:,:))'; % 1st field, c_mr
RS_spec_m = sq(RS_mat(2,:,:))'; %  2nd field, s_m
RS_spec_r = sq(RS_mat(3,:,:))'; %  3rd field, s_r

STIM_coh_mr = sq(stim_mat(1,:,:))'; % 1st field, c_sr
STIM_spec_m = sq(stim_mat(4,:,:))'; %  4th field, dmt spec_s
STIM_spec_r = sq(stim_mat(5,:,:))'; %  5th field, dmt spec_r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN ABS COHERENCES and SPECTRUM vs FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = size(RS_coh_mr,1); % -- number of electrodes
data.num_el = M;

% ----------- COHERENCE RS

% --- mean coherences
data.RS_mean_coh_mr = mean(abs(RS_coh_mr));  % modulator - receiver
% --- std coherences
data.RS_std_coh_mr = std(abs(RS_coh_mr)); % modulator - receiver
% --- Error bars
data.RS_err_mr = data.RS_std_coh_mr/sqrt(M);

% ----------- COHERENCE STIM

% --- mean coherences
data.STIM_mean_coh_mr = mean(abs(STIM_coh_mr));  % modulator - receiver
% --- std coherences
data.STIM_std_coh_mr = std(abs(STIM_coh_mr)); % modulator - receiver
% --- Error bars
data.STIM_err_mr = data.STIM_std_coh_mr/sqrt(M);

% ----------- SPECTRUM RS

% --- mean spectrum
data.RS_mean_spec_m = mean(RS_spec_m);  % mean spectrum modulator
data.RS_mean_spec_r = mean(RS_spec_r);  % mean spectrum receiver

% --- std spectrum
data.RS_std_spec_m = std(RS_spec_m);  
data.RS_std_spec_r = std(RS_spec_r);  

% --- Error bars spectrum
data.RS_err_S_m = data.RS_std_spec_m/sqrt(M); % -- all
data.RS_err_S_r = data.RS_std_spec_r/sqrt(M); % Note: the receiver should be less than M, but for the way they are saved they are M = 41. It's because for each modulator I save the receiver as well. Therefore some receivers repeats


% ----------- SPECTRUM STIM

% --- mean spectrum
data.STIM_mean_spec_m = mean(STIM_spec_m);  % mean spectrum modulator
data.STIM_mean_spec_r = mean(STIM_spec_r);  % mean spectrum receiver

% --- std spectrum
data.STIM_std_spec_m = std(STIM_spec_m);  
data.STIM_std_spec_r = std(STIM_spec_r);  

% --- Error bars spectrum
data.STIM_err_S_m = data.STIM_std_spec_m/sqrt(M); % -- all
data.STIM_err_S_r = data.STIM_std_spec_r/sqrt(M); % Note: the receiver should be less than M, but for the way they are saved they are M = 41. It's because for each modulator I save the receiver as well. Therefore some receivers repeats


end