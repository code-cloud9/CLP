clear;

% Parameters
graphon_setting = 4;
node_num = 200;
missing_prob = 0.1;
train_proportion = 0.4;
r0 = 50;  
alpha_1 = 0.1;  % for BH procedure
alpha_2 = 0.2;  % for eBH procedure
m_deRand = 20;  % # replications for derandomization
sig = 0.0:0.05:0.2; % threshold

% Initialize the result table
reps = 100;  
nsig=length(sig); 
res_example = table(nan(reps, nsig), nan(reps, nsig), nan(reps, nsig), nan(reps, nsig), 'VariableNames', {'TP', 'FP', 'FN', 'TN'});

% Loop for #reps repetitions
for rep_idx = 1:reps

    disp(rep_idx);
    TP = zeros(1,nsig);
    FP = zeros(1,nsig);
    FN = zeros(1,nsig);
    TN = zeros(1,nsig);

    for idx_sig = 1:nsig

    [matrix_A_eg, matrix_M_eg, matrix_C_eg] = gen_data_unic(graphon_setting, node_num, sig(idx_sig), missing_prob);  

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

                ones_idx_calib_i0 = find(vector_calib_idx == 1);
                testset = find(vector_test_idx == 1);
                ntesti0 = length(testset);
                est_res = ests_new(node_num, i0, matrix_A_eg, matrix_M_eg, vector_train_idx, testset, ones_idx_calib_i0);
                ptestset = zeros(1,ntesti0);
                track_j0_count =1 ;
                for j0 =testset
                    ptestset(track_j0_count) = pvalue_single(i0, j0, matrix_A_eg, matrix_C_eg, est_res, vector_calib_idx) ;
                    track_j0_count = track_j0_count + 1;
                end
                
             evalues_i0 = eValues(ptestset, alpha_1);
            
            
            non_na_count = sum(~isnan(evalues_i0));
            e_list(loc_track+1 : loc_track+non_na_count) = evalues_i0(~isnan(evalues_i0)) ;
            loc_track = loc_track + non_na_count;
        end
        
            e_list_deRand = e_list_deRand + e_list;
    end
    
    e_list_deRand = e_list_deRand / m_deRand;

    % ----- derandomize get many e_list and take the average ------^^^

    batch_size = node_num^2;
    eBH_rej_set = eBH(1/alpha_1 * e_list_deRand, alpha_2, batch_size);

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
    TP(idx_sig) = sum(sum(conditionTrue & rej_location));
    FP(idx_sig) = sum(sum(conditionFalse & rej_location));
    FN(idx_sig) = sum(sum(conditionTrue  & ~rej_location));
    TN(idx_sig) = sum(sum(conditionFalse & ~rej_location));

    end

    res_example(rep_idx,:) = {TP,FP,FN,TN};

 end


% Save results (does not need here)
% save(['clp_new_unic_stg' num2str(graphon_setting) '.mat'], 'res_example');


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
