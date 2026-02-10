function est_res = ests_new(n, i0, A, M, train_idx_i0_row, j0, calib_idx_set_j0)

    M(logical(eye(size(M)))) = 1; 
    est_res = nan(1, n);
    combined_cols = [calib_idx_set_j0,j0];
    train_set = find(train_idx_i0_row==1);

    for i = combined_cols
        hat_d = zeros(1,n);
        for j = train_set
            hat_d(j) =  1 ;
        end
         weights = exp(-hat_d(train_set).^2 / 2) / sqrt(2*pi);
         est_res(i) = sum(A(i0, train_set).*weights) / ...
                     sum(weights);
    end
    
    
        
    
    

