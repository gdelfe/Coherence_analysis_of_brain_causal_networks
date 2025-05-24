% estimate_var_model.m
function [A, SIG] = estimate_var_model(X, morder, regmode)
    ptic('\n*** tsdata_to_var... ');
    [A, SIG] = tsdata_to_var(X, morder, regmode);
    ptoc;
    assert(~isbad(A), 'VAR estimation failed');
end