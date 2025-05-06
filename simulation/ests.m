function est_res = ests(n, i0, A, M, train_idx_i0_row, j0, calib_idx_set_j0)
    % i0: index of the selected missing entry's row
    % A: the observed adjacency matrix with missings = 0
    % M: missing pattern (1 indicates missing)
    % train_idx_i0_row: a vector of training indices in i0 row (1 for chosen)
    % j0: the choosen test index on i0 row
    % calib_idx_set_j0: the index set of calibration for j0
    %
    % return: est_res (R^n): estimation of i0 row's (calib_idx_set_j0) U (j0)
    % using Gaussian kernel, other position = nan

    M(logical(eye(size(M)))) = 1;

    % Initialize output
    est_res = nan(1, n);

    % Construct indicator Omega_j0
    Omega_j0 = train_idx_i0_row;
    Omega_j0(calib_idx_set_j0)=true;
    Omega_j0(j0)=true;
    Omega_j0(i0)=false;
    Omega_j0 = logical(Omega_j0);

    Calib_plus = [calib_idx_set_j0,j0];
    Calib_plus = Calib_plus(Calib_plus~=i0);

    % Construct matched support
    matchIndices = all(M(Omega_j0, Calib_plus) == 0, 2);
    support = find(Omega_j0==1);
    Idx_support = support(matchIndices);
    nIdx = length(support);
    
    % Compute distance and estimation:
    combined_cols = [calib_idx_set_j0,j0];
    hat_d=zeros(n,n);
    train_set = find(train_idx_i0_row==1);
    ntrain = length(train_set);

    for i = combined_cols
        hat_d = zeros(1,n);
        for j = train_set
            diff = A(Idx_support,i)-A(Idx_support,j);
            temptrain = train_set(train_set~=j);
            hat_d(j) = sum(abs(diff'*(A(Idx_support,temptrain ) )),"all")/((ntrain-1)*nIdx);
        end
         weights = exp(-hat_d(train_set).^2 / 2) / sqrt(2*pi); 
         est_res(i) = sum(A(i0, train_set).*weights) / ...
                     sum(weights);
    end
    
    
        
    
    

