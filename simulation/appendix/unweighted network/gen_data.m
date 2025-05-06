function [matrix_A, matrix_M, matrix_C] = gen_data(setting, n, sig, mis_prop, rho, isBinary,isHetero, isDiag)
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

    if nargin<6
        isBinary = false;
    end

     if nargin<7
        isHetero = false;
    end

    if nargin< 8
        isDiag = 0;
    end

    if setting == 1 

        % graphon
        xi = rand(1, n);
        xj = rand(1, n);
        matrix_f = xi'^3 + 2 * xj^3;
        [matrix_A, matrix_M, matrix_C] = gen_data_helper(alt_spar, matrix_f, n, sig, mis_prop, rho, isBinary,isHetero, isDiag);


    elseif setting == 2 

        % graphon
        xi = rand(1, n);
        denominator = (2*xi' - 0.5).^3 + (xi - 0.5).^3 + 0.01;
        matrix_f = cos(0.1 ./ denominator) .* max(xi', xi).^(2/3);
        [matrix_A, matrix_M, matrix_C] = gen_data_helper(alt_spar, matrix_f, n, sig, mis_prop, rho, isBinary,isHetero, isDiag);


    elseif setting == 3 
        
        % graphon
        xi = rand(1, n);
        matrix_f = (3*xi'.^2 + xi.^2) .* cos(1 ./ (2*xi'.^4 + xi.^4));
        [matrix_A, matrix_M, matrix_C] = gen_data_helper(alt_spar, matrix_f, n, sig, mis_prop, rho, isBinary,isHetero, isDiag);

    else 

        error('Invalid setting');
    end
end


