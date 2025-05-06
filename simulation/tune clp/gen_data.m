function [matrix_A, matrix_M, matrix_C] = gen_data(setting, n, sig, mis_prop, rho)
    % n: number of nodes
    % mis_prop: missing proportion, defaults to 0.2
    % rho: optional sparse parameter, defaults to 1
    % sig: signal strength
    % return: matrix_A: the true adjacency matrix without missing entries, diagonal equals 0 as no self-loops
    %         matrix_M: missing pattern (1 indicates missing), diagonal equals 0 as no self-loops
    %         matrix_C: the null values
    
    alt_spar = 0.3; % sparsity of alternative
    if nargin < 4
        mis_prop = 0.2;
    end

    if nargin < 5
        rho = 1;
    end

    if setting == 1 
        
        % generate A
        xi = rand(1, n);
        matrix_temp = (xi').^3 + 2 * xi.^3;
        matrix_f = matrix_temp + 0.2 * rand(n) - 0.1;
        matrix_A = rho * matrix_f;
        diag_indices = 1:n+1:n^2; 
        matrix_A(diag_indices) = 0; 

        % generate M (1 indicates missing)
        matrix_q = mis_prop * ones(n, n);  
        matrix_M = binornd(1, matrix_q);
        matrix_M(diag_indices) = 0;

        % generate C
        matrix_Alt = binornd(1,alt_spar*ones(n,n));
        matrix_C = (matrix_M).*matrix_A-sig*(matrix_M).*matrix_Alt;  

    elseif setting == 2
        
        % generate A
        xi = rand(1, n);
        denominator = (2*xi' - 0.5).^3 + (xi - 0.5).^3 + 0.01;
        matrix_temp = cos(0.1 ./ denominator) .* max(xi', xi).^(2/3);
        matrix_f = matrix_temp + 0.2 * rand(n) - 0.1;
        matrix_A = rho * matrix_f;
        diag_indices = 1:n+1:n^2; 
        matrix_A(diag_indices) = 0; 

        % generate M (1 indicates missing)
        matrix_q = mis_prop * ones(n, n);  
        matrix_M = binornd(1, matrix_q);
        matrix_M(diag_indices) = 0; 

        % generate C
        matrix_Alt = binornd(1,alt_spar*ones(n,n));
        matrix_C = (matrix_M).*matrix_A-sig*(matrix_M).*matrix_Alt; % 

    elseif setting == 3 

        % generate A
        xi = rand(1, n);
        matrix_temp = (3*xi'.^2 + xi.^2) .* cos(1 ./ (2*xi'.^4 + xi.^4));
        matrix_f = matrix_temp + 0.2 * rand(n) - 0.1;
        matrix_A = rho * matrix_f;
        diag_indices = 1:n+1:n^2; 
        matrix_A(diag_indices) = 0; 

        % generate M (1 indicates missing)
        matrix_q = mis_prop * ones(n, n);  
        matrix_M = binornd(1, matrix_q);
        matrix_M(diag_indices) = 0; 

        % generate C
        matrix_Alt = binornd(1,alt_spar*ones(n,n));
        matrix_C = (matrix_M).*matrix_A-sig*(matrix_M).*matrix_Alt; % 

    else 

        error('Invalid setting');
    end
end


