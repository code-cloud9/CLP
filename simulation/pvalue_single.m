function pvalue_i0j0 = pvalue_single(i0, j0, A, C, est_res_i0_row_j0split, calib_idx_j0)
    % i0, j0: indices of the selected missing entry
    % A: the observed adjacency matrix with missings = 0
    % C: the null values
    % est_res_i0_row_j0split (R^n): estimation of calib and teat set in i0 row, NA
    %                       if missing or not in (calib_idx_set_j0) U (j0)
    % calib_idx_j0 (R^n): a vector of calib indices in i0 row (1 for chosen) for j0
    %
    % return: p-value for (i0,j0)

    % Calculate calib_scores
    calib_scores = A(i0, :) - est_res_i0_row_j0split;
    calib_idx = calib_idx_j0 == 1;
    calib_scores_j0 = calib_scores(calib_idx);

    % Calculate test_score_at_i0j0
    test_score_j0 = C(i0, j0) - est_res_i0_row_j0split(j0); 

    % Compute p-value for (i0,j0)
    ncalib = length(calib_scores_j0);
    pvalue_i0j0 = (1 + sum(calib_scores_j0 < test_score_j0)) / (ncalib + 1);
end


