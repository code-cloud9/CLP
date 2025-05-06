function [matrix_A, matrix_M, matrix_C] = gen_data_unic(setting, n, sig, mis_prop, rho)
    % n: number of nodes
    % mis_prop: missing proportion, defaults to 0.2
    % rho: optional sparse parameter, defaults to 1
    % sig: signal-th quantile of the observed outcomes
    % return: matrix_A: the true adjacency matrix without missing entries, diagonal equals 0 as no self-loops
    %         matrix_M: missing pattern (1 indicates missing), diagonal equals 0 as no self-loops
    %         matrix_C: the null values
    
    if nargin < 4
        mis_prop = 0.2;
    end

    if nargin < 5
        rho = 1;
    end

    if setting == 4 % directed unweighted network

        % generate A
        xi = rand(1,n);
        matrix_temp = 1/2*(xi' + xi);
        matrix_f = matrix_temp + 0.5*(rand(n, n)-0.5);
        matrix_A = rho *  2*(matrix_f > sig) - 1; 

        diag_indices = 1:n+1:n^2; 
        matrix_A(diag_indices) = 0; 

        % generate M (1 indicates missing)
        matrix_q = mis_prop * ones(n, n);  
        hetero = 2*rand(n,n);
        matrix_M = binornd(1, hetero.*matrix_q);
        matrix_M(diag_indices) = 0; 

        % generate C
        matrix_Alt = 0.00*binornd(1,ones(n,n)); 
        matrix_C = (matrix_M).*matrix_Alt; 

    else 

        error('Invalid setting');
    end
end


