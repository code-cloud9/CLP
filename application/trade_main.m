%% clp
clear;  

% Parameters
matrix_A_eg = readmatrix('trade_matrix.csv');
node_num = size(matrix_A_eg, 2);
n = node_num;
diag_indices = 1:n+1:n^2;
matrix_A_eg(diag_indices) = 0;

missing_prob = 0.1;
train_proportion = 0.4;
r0 = 25;  
alpha_1 = 0.1;  % for BH procedure {0.05,0.1,0.15}
alpha_2 = 0.2;  % for eBH procedure {0.1,0.2,0.3}
m_deRand = 20;  % # replications for derandomization
c_null_value = 0.2;

% Initialize the result table
reps = 100;  
res_example = table(nan(reps, 1), nan(reps, 1), nan(reps, 1), nan(reps, 1), 'VariableNames', {'TP', 'FP', 'FN', 'TN'});

% Loop for #reps repetitions
for rep_idx = 1:reps

    matrix_q = missing_prob * ones(n, n);
    matrix_M_eg = binornd(1, matrix_q);
    matrix_M_eg(diag_indices) = 0; 

    matrix_C_eg = ones(n, n).* c_null_value; 
    row_indices = find(any(matrix_M_eg == 1, 2));
    
    % ----- derandomize: get #m_deRand e_list(s) and take the average ------vvv
    
    e_list_deRand = zeros(1, sum(matrix_M_eg(:)));
    e_total = sum(matrix_M_eg(:)) ;
    
    parfor m = 1:m_deRand
        
        e_list = nan(1,e_total);
        loc_track = 0;
        
        for i0 = row_indices'
              
            [vector_train_idx, vector_calib_idx, vector_test_idx] = gen_partition(node_num, ...
                i0, matrix_M_eg, train_proportion);
                 
             r1 = floor(sum(vector_calib_idx) / r0); 
             r2 = ceil(sum(vector_test_idx) / r1);

            splits_test = cell(1, r2);
            for i = 1:r2
                splits_test{i} = zeros(1, length(vector_test_idx));
            end

            ones_idx_test = find(vector_test_idx == 1);
            for i = 1:length(ones_idx_test)
                split_idx_test = mod(i-1, r2) + 1;
                splits_test{split_idx_test}(ones_idx_test(i)) = 1;
            end
            
            evalues_i0 = nan(1, node_num);
            for ii = 1:r2
                one_split_test = find(splits_test{ii} == 1); 
                pvalues_inaTestsplit = nan(1, sum(splits_test{ii}));

                num_parti = sum(splits_test{ii}); 
                splits_calib_j0 = cell(1, num_parti);
                for i = 1:num_parti
                    splits_calib_j0{i} = zeros(1, length(vector_calib_idx));
                end
                ones_idx_calib_j0 = find(vector_calib_idx == 1);
                ones_idx_shuffled = ones_idx_calib_j0(randperm(length(ones_idx_calib_j0)));
                for i = 1:length(ones_idx_shuffled)
                    split_idx_calib_j0 = mod(i-1, num_parti) + 1;
                    splits_calib_j0{split_idx_calib_j0}(ones_idx_shuffled(i)) = 1;
                end
                
                track_calib_split = 1;
                track_j0_count = 1;
                for j0 = one_split_test
                    calib_idx_j0 = splits_calib_j0{track_calib_split};
                    track_calib_split = track_calib_split + 1;
            
                    calib_idx_set_j0 = find(calib_idx_j0 == 1);
                    est_res_j0 = ests(node_num, i0, matrix_A_eg, matrix_M_eg, vector_train_idx, j0, calib_idx_set_j0);
                    pvalue_i0j0 = pvalue_single(i0, j0, matrix_A_eg, matrix_C_eg, est_res_j0, calib_idx_j0);
                    pvalues_inaTestsplit(track_j0_count) = pvalue_i0j0;
                    track_j0_count = track_j0_count + 1;
                end 

                evalues_inaTestsplit = eValues(pvalues_inaTestsplit, alpha_1);
                evalues_i0(one_split_test) = evalues_inaTestsplit;
            end
            
            non_na_count = sum(~isnan(evalues_i0));
            e_list(loc_track+1 : loc_track+non_na_count) = evalues_i0(~isnan(evalues_i0));
            loc_track = loc_track + non_na_count;
        end
        e_list_deRand = e_list_deRand + e_list;
    end
    
    e_list_deRand = e_list_deRand / m_deRand;
    
    % ----- derandomize get many e_list and take the average ------^^^

    
    batch_size = node_num^2; 
    eBH_rej_set = eBH(1.0/alpha_1 * e_list_deRand, alpha_2, batch_size);

    [row, col] = find(matrix_M_eg == 1);

    temp_mis_loc = [row, col];
    missingLocation = sortrows(temp_mis_loc, [1 2]);
    
    df_missing_idx = table((1:size(missingLocation, 1))', missingLocation(:, 1), ...
        missingLocation(:, 2), 'VariableNames', {'id', 'row', 'colmn'});
    
    rej_location = zeros(node_num, node_num);
    

    for i = eBH_rej_set
        rej_location(df_missing_idx.row(i), df_missing_idx.colmn(i)) = 1;
    end
    
    conditionTrue = matrix_M_eg.*matrix_A_eg > matrix_M_eg.*matrix_C_eg;
    conditionFalse = matrix_M_eg.*(~conditionTrue);
    TP = sum(sum(conditionTrue & rej_location));
    FP = sum(sum(conditionFalse & rej_location));
    FN = sum(sum(conditionTrue  & ~rej_location));
    TN = sum(sum(conditionFalse & ~rej_location));

    res_example.TP(rep_idx) = TP;
    res_example.FP(rep_idx) = FP;
    res_example.FN(rep_idx) = FN;
    res_example.TN(rep_idx) = TN;

end


% Calculate FDP (False Discovery Proportion)
fdp = mean(res_example.FP ./ (res_example.TP + res_example.FP),1);
fdp(isnan(fdp)) = 0;
disp(fdp);

% Calculate power
power = mean(res_example.TP ./ (res_example.TP + res_example.FN),1);
disp(power);



%% cmc

clear;

% Parameters
matrix_A_eg = readmatrix('trade_matrix.csv');
node_num = size(matrix_A_eg, 2);
n = node_num;
diag_indices = 1:n+1:n^2; 
matrix_A_eg(diag_indices) = 0;

missing_prob = 0.1;
train_proportion = 0.4;
alpha_2 = 0.1;  % target {0.1,0.2,0.3}
c_null_value = 0.2;

% Initialize the result table
reps = 100;  
res_example = table(nan(reps, 1), nan(reps, 1), nan(reps, 1), nan(reps, 1), 'VariableNames', {'TP', 'FP', 'FN', 'TN'});

r = 2;
for rep_idx = 1:reps

    disp(rep_idx);

    matrix_q = missing_prob * ones(n, n);  
    matrix_M_eg = binornd(1, matrix_q);
    matrix_M_eg(diag_indices) = 0;
    matrix_C_eg = ones(n, n).* c_null_value; 

    % cmc procedure
    matrix_q = train_proportion * ones(node_num, node_num);  
    matrix_Train = binornd(1, matrix_q).*(1-matrix_M_eg);
    
    matrix_Cal = (1-matrix_Train).*(1-matrix_M_eg); 

    hat_A = lr_init(matrix_A_eg, 1-matrix_Train, r);

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
    TP = sum(sum(conditionTrue & rej_location));
    FP = sum(sum(conditionFalse & rej_location));
    FN = sum(sum(conditionTrue  & ~rej_location));
    TN = sum(sum(conditionFalse & ~rej_location));

    res_example.TP(rep_idx) = TP;
    res_example.FP(rep_idx) = FP;
    res_example.FN(rep_idx) = FN;
    res_example.TN(rep_idx) = TN;

 end


% Calculate FDP (False Discovery Proportion)
fdp = mean(res_example.FP ./ (res_example.TP + res_example.FP),1);
fdp(isnan(fdp)) = 0;
disp(fdp);

% Calculate power
power = mean(res_example.TP ./ (res_example.TP + res_example.FN),1);
disp(power);




%% dss

clear;

% Parameters
matrix_A_eg = readmatrix('trade_matrix.csv');
c_null_value = 0.2;
matrix_A_eg(matrix_A_eg < 0.2) = -1;
matrix_A_eg(matrix_A_eg >= 0.2) = 1;
node_num = size(matrix_A_eg, 2);
n = node_num;
diag_indices = 1:n+1:n^2;
matrix_A_eg(diag_indices) = 0;

missing_prob = 0.1;
train_proportion = 0.4;
alpha_2 = 0.1;  % target {0.1,0.2,0.3}
c_null_value = 0.2;

% Initialize the result table
reps = 100;  
res_example = table(nan(reps, 1), nan(reps, 1), nan(reps, 1), nan(reps, 1), 'VariableNames', {'TP', 'FP', 'FN', 'TN'});

r = 2;
% Loop for #reps repetitions
for rep_idx = 1:reps

    disp(rep_idx);

    matrix_q = missing_prob * ones(n, n); 
    matrix_M_eg = binornd(1, matrix_q);
    matrix_M_eg(diag_indices) = 0; 

    % generate C
    matrix_C_eg = ones(n, n).* c_null_value; 

    % dss procedure
    obs_neg = ((1-matrix_M_eg).*matrix_A_eg)<0;
    matrix_q = 0.5 *obs_neg;
    matrix_Cal = binornd(1, matrix_q).*(1-matrix_M_eg);
    matrix_Train = (1-matrix_Cal).*(1-matrix_M_eg); 
    hat_A = lr_init(matrix_A_eg, 1-matrix_Train, r); 
    CalA = (matrix_Cal.*(hat_A)); 
    Calset =  [CalA(matrix_Cal>0)' ,Inf];
    Ntotal =length(Calset);

    Test_M =  matrix_M_eg.*(hat_A );
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
    TP = sum(sum(conditionTrue & rej_location));
    FP = sum(sum(conditionFalse & rej_location));
    FN = sum(sum(conditionTrue  & ~rej_location));
    TN = sum(sum(conditionFalse & ~rej_location));

    res_example.TP(rep_idx) = TP;
    res_example.FP(rep_idx) = FP;
    res_example.FN(rep_idx) = FN;
    res_example.TN(rep_idx) = TN;

 end

% Calculate FDP (False Discovery Proportion)
fdp = mean(res_example.FP ./ (res_example.TP + res_example.FP),1);
fdp(isnan(fdp)) = 0;
disp(fdp);

% Calculate power
power = mean(res_example.TP ./ (res_example.TP + res_example.FN),1);
disp(power);














%%
clear;

reps = 100; 
trade_result = nan(reps, 12);

loaded_data = load('clp_01.mat');
res_example = loaded_data.res_example;
trade_result(:, 1) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 7) = res_example.TP ./ (res_example.TP + res_example.FN);

loaded_data = load('clp_02.mat');
res_example = loaded_data.res_example;
trade_result(:, 3) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 9) = res_example.TP ./ (res_example.TP + res_example.FN);

loaded_data = load('clp_03.mat');
res_example = loaded_data.res_example;
trade_result(:, 5) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 11) = res_example.TP ./ (res_example.TP + res_example.FN);

loaded_data = load('cmc_01_r2.mat');
res_example = loaded_data.res_example;
trade_result(:, 2) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 8) = res_example.TP ./ (res_example.TP + res_example.FN);

loaded_data = load('cmc_02_r2.mat');
res_example = loaded_data.res_example;
trade_result(:, 4) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 10) = res_example.TP ./ (res_example.TP + res_example.FN);

loaded_data = load('cmc_03_r2.mat');
res_example = loaded_data.res_example;
trade_result(:, 6) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 12) = res_example.TP ./ (res_example.TP + res_example.FN);

loaded_data = load('dss_01_r2.mat');
res_example = loaded_data.res_example;
trade_result(:, 13) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 16) = res_example.TP ./ (res_example.TP + res_example.FN);

loaded_data = load('dss_02_r2.mat');
res_example = loaded_data.res_example;
trade_result(:, 14) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 17) = res_example.TP ./ (res_example.TP + res_example.FN);

loaded_data = load('dss_03_r2.mat');
res_example = loaded_data.res_example;
trade_result(:, 15) = res_example.FP ./ (res_example.TP + res_example.FP);
trade_result(:, 18) = res_example.TP ./ (res_example.TP + res_example.FN);


writematrix(trade_result, 'trade_result_with_dss.csv');
