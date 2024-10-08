
% Function to compute the granger test, across all trials, for two time
% series X and Y, which must have the same number of trials  (always the
% case in our experiments).
% OUTPUT: h_rate: rate of GC test passed across trials

function h_rate = granger_causality_X_Y(lfp_X, lfp_Y,lag)

% Assum% Add the MVGC toolbox to the MATLAB path
addpath('/vol/bd5/People/Gino/Coherence_modulator_analysis/Gino_codes');

numTrials = size(lfp_X, 1); % Number of trials

h = 0;
for trial = 1:numTrials

    % Extract data for the current trial. MVGC expects data in the form [variables x time points x trials].
    % Since we're processing one trial at a time in this loop, trials dimension will be 1.
    X = lfp_X(trial, :)'; % X data
    Y = lfp_Y(trial, :)'; % Y data

    gc_h = gctest(X,Y,'NumLags',lag); % gc_h is 1 if the test is passed, 0 otherwise
    h = h + gc_h;

end

h_rate = h/numTrials;

end 