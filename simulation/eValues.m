function evalues = eValues(pvals, alpha)
    % pvals: for all edges in a split of test set (I_test^i0(k0))
    % alpha: significance level for calculating e-values, default to 0.1
    % return: e_values from a single local split test (I_test^i0(k0))

    if nargin < 2
        alpha = 0.1;
    end

    ntest = length(pvals);

    % Initialize output
    evalues = zeros(1, ntest);
    
     % BH
    [sortedPvals, sortedIdx] = sort(pvals);
    threshold = alpha * (1:ntest) / ntest;
    
    l_set = find(sortedPvals <= threshold);
    
    if isempty(l_set)
        R_i0j0 = [];
    else
        R_i0j0 = sortedIdx(1:max(l_set));
    end

    evalues(R_i0j0) = ntest / max(length(R_i0j0), 1) / alpha;

end

