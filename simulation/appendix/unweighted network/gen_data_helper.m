function [matrix_A, matrix_M, matrix_C] = gen_data_helper(alt_spar, matrix_f, n, sig, mis_prop, rho, isBinary,isHetero, isDiag)
    
    matrix_q = mis_prop * ones(n, n);  

    diag_indices = 1:n+1:n^2; 

    if isHetero
        hetero = 2*rand(n,n);
        matrix_M = binornd(1, hetero.*matrix_q);
    else
        matrix_M = binornd(1, matrix_q);
    end
    matrix_M(diag_indices) = 0; 

    if ~isBinary
        matrix_A =  rho * (matrix_f + 0.2 * rand(n) - 0.1);
    else
	matrix_f = matrix_f + 0.2 * rand(n) - 0.1;
        matrix_A =  (2 * rho * rescale(matrix_f)-1).*(matrix_M) + ...
            (2 * binornd(1,rho * rescale(matrix_f))-1).*(1- matrix_M);
    end 
    matrix_A(diag_indices) = isDiag * matrix_A(diag_indices); 

    matrix_Alt = binornd(1,alt_spar*ones(n,n)); 
    matrix_C = (matrix_M).*matrix_A-sig*(matrix_M).*matrix_Alt;

