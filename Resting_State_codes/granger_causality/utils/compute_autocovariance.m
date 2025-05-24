% compute_autocovariance.m
function [G, info] = compute_autocovariance(A, SIG, acmaxlags)
    ptic('*** var_to_autocov... ');
    [G, info] = var_to_autocov(A, SIG, acmaxlags);
    ptoc;
end