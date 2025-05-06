clear

% Parameters
graphon_setting = 1; % {1,2,3} 
node_num = 200;
missing_prob = 0.2;
train_proportion = 0.4;
alpha_2 = 0.2; 
sig = 0.1:0.1:0.7; % quantile

% Initialize the result table
reps = 100;  
nsig=length(sig); 
res_example = table(nan(reps, nsig), nan(reps, nsig), nan(reps, nsig), nan(reps, nsig), 'VariableNames', {'TP', 'FP', 'FN', 'TN'});

r = 2;
% Loop for #reps repetitions
parfor rep_idx = 1:reps

    disp(rep_idx);
    TP = zeros(1,nsig);
    FP = zeros(1,nsig);
    FN = zeros(1,nsig);
    TN = zeros(1,nsig);

    [matrix_A_eg, matrix_M_eg, matrix_C_eg] = gen_data_unic(graphon_setting, node_num, 1, missing_prob);  

    matrix_q = train_proportion * ones(node_num, node_num);  
    matrix_Train = binornd(1, matrix_q).*(1-matrix_M_eg);

    
    matrix_Cal = (1-matrix_Train).*(1-matrix_M_eg);

    hat_A = lr_init(matrix_A_eg, 1-matrix_Train, r); 

    for idx_sig = 1:nsig
    
    matrix_C_eg = sig(idx_sig) * ones(node_num,node_num);

    ResdueR = (matrix_Cal.*(hat_A-matrix_A_eg)); 

    Calset =  [ResdueR(ResdueR>0)' ,Inf];
    Ntotal =length(Calset);

    Test_M =  matrix_M_eg.*(hat_A - matrix_C_eg);
    Testvar = Test_M(Test_M~=0);
    Ntest = length(Testvar);

    ptest = zeros(1,Ntest);
    for i = 1:Ntest
        ptest(i) = sum(Calset>Testvar(i) )/Ntotal;
    end

    [sortedPvals, sortedIdx] = sort(ptest);
    threshold = alpha_2 * (1:Ntest) / Ntest;
    
    l_set = find(sortedPvals <= threshold);
    
    if isempty(l_set)
        rej_set = [];
    else
        rej_set = sortedIdx(1:max(l_set));
    end

    [row, col] = find(matrix_M_eg == 1);
    temp_mis_loc = [row, col];
    missingLocation = sortrows(temp_mis_loc, [1 2]);
    df_missing_idx = table((1:size(missingLocation, 1))', missingLocation(:, 1), ...
        missingLocation(:, 2), 'VariableNames', {'id', 'row', 'colmn'});
    
    rej_location = zeros(node_num, node_num);
    

    for i = rej_set
        rej_location(df_missing_idx.row(i), df_missing_idx.colmn(i)) = 1;
    end
    
    conditionTrue = matrix_M_eg.*matrix_A_eg > matrix_M_eg.*matrix_C_eg;
    conditionFalse = matrix_M_eg.*(~conditionTrue);
    TP(idx_sig) = sum(sum(conditionTrue & rej_location));
    FP(idx_sig) = sum(sum(conditionFalse & rej_location));
    FN(idx_sig) = sum(sum(conditionTrue  & ~rej_location));
    TN(idx_sig) = sum(sum(conditionFalse & ~rej_location));

    end

    res_example(rep_idx,:) = {TP,FP,FN,TN};

end


% Save results
save(['cmc_unic_stg' num2str(graphon_setting) '.mat'], 'res_example');

% Calculate fdp
ratio1 = res_example.FP ./ (res_example.TP + res_example.FP);
ratio1(isnan(ratio1)) = 0;
fdp = mean(ratio1,1);
disp(fdp);

fdpsderr = std(ratio1, 1) / 10; 
disp(fdpsderr);

% Calculate power
ratio2 = res_example.TP ./ (res_example.TP + res_example.FN);
ratio2(isnan(ratio2)) = 0;
power = mean(ratio2,1);
disp(power);

powersderr = std(ratio2, 1) / 10; 
disp(powersderr);

