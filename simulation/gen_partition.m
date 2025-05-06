function [train_idx, calib_idx, test_idx] = gen_partition(n, i0, M, ptrain)
    % n: number of nodes
    % i0: index of the selected missing entry's row
    % M: missing pattern (1 indicates missing)
    % ptrain: proportion of the training set in observed entries in i0 row
    %
    % return: two vectors of training/calib indices in i0 row (1 for chosen)
    %         one vector of missing indices in i0 row (1 for missing)
    
    if nargin < 4  
        ptrain = 0.5;
    end

    % Generate the partitions
    nomissing = 1 - M(i0, :);
    nomissing(i0) = 0;  % no self-loops position in training set
    bi_idx = binornd(1, ptrain, [1, n]);  % binomial random variables
    train_idx = nomissing .* bi_idx;
    calib_idx = nomissing .* (1 - bi_idx);
    test_idx = M(i0, :);
end

