function Est_M = lr_init(matrix_A, matrix_M, r)
% Apply GD to estimate a low rank-r matrix
% matrix_A total matrix
% matrix_M missing matrix

Omega = 1-matrix_M;
Obv = Omega.*matrix_A;
[n1,n2] = size(matrix_M);
p = (n1*n2-sum(matrix_M,"all"))/(n1*n2);
[X,~,Y] = svds(1/p*Obv,r);
Epoch = 6000;
lamb = 0.1;
eta = 0.0001;
res = zeros(1,Epoch);
resp = Inf;
for i = 1:Epoch

    Xt = X - eta/p*( Omega.*(X*Y'-Obv)*Y +lamb*X );
    Yt = - eta/p*( (Omega.*(X*Y'-Obv))'*X +lamb*Y );
    X = Xt;
    Y = Yt;

    res(i) = norm(matrix_A-X*Y',"fro");
    if res(i)>resp
        break;
    end
    resp = res(i);
end
Est_M = X*Y';
end
