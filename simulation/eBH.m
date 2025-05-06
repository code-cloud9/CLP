function reject_Total = eBH(evals, alpha , Mbatch_size)
    % evals: a vector of e-values with length equals length(missing)
    % alpha: significance level for eBH, default to 0.2
    % batch_size: maximum batch size of each local evalue rejection
    % returns: a vector of eBH(alpha) rejected tests in the whole test set
    
    if nargin <3
        Mbatch_size = length(evals);
    if nargin < 2
        alpha = 0.2;
    end
    end
    n = length(evals);
    nbatch =  ceil(n/Mbatch_size);
    Partition = cell(1, nbatch);
    for i =1:nbatch
     Partition{i}=i:nbatch:n;
    end
    
    reject_Total = [];
    
    for i =1:nbatch
    
    Idx_batch = Partition{i};
    evals_batch = evals(Idx_batch);

    ntest = length(evals_batch);
    [sortedEvals, sortedIdx] = sort(evals_batch, 'descend');
    thresholds = (ntest / alpha) ./ (1:ntest);

    k_set = find(sortedEvals >= thresholds);
    
    if isempty(k_set)
        rejectedTests = [];
    else
        max_k = max(k_set);
        rejectedTests = sortedIdx(1:max_k);
    end
    end

    reject_Total = [reject_Total, Idx_batch(rejectedTests) ];

end

