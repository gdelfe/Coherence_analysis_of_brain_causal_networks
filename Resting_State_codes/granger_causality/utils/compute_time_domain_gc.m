
% compute_time_domain_gc.m
function [F, pval, sig] = compute_time_domain_gc(G, X, morder, nobs, ntrials, tstat, alpha, mhtc)
    ptic('*** autocov_to_pwcgc... ');
    F = autocov_to_pwcgc(G);
    ptoc;
    assert(~isbad(F,false), 'GC calculation failed');
    nvars = size(X,1);
    pval = mvgc_pval(F, morder, nobs, ntrials, 1, 1, nvars-2, tstat);
    sig = significance(pval, alpha, mhtc);
end