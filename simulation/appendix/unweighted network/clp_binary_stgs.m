clear;

% Parameters
graphon_setting = 1; % {1,2,3}
node_num = 200;
missing_prob = 0.2;
train_proportion = 0.4;
r0 = 25;  
alpha_1 = 0.1;  % for BH procedure
alpha_2 = 0.2;  % for eBH procedure
m_deRand = 20;  % # replications for derandomization
sig = 0.1:0.1:0.7; % signal strength

isBinary = true;
isHetero = true;

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
        
        [matrix_A_eg, matrix_M_eg, matrix_C_eg] = gen_data(graphon_setting, node_num, sig(idx_sig), missing_prob,1,isBinary,isHetero);  
  
        row_indices = find(any(matrix_M_eg == 1, 2));
        
        % ----- derandomize: get #m_deRand e_list(s) and take the average ------vvv
        
        e_list_deRand = zeros(1, sum(matrix_M_eg(:)));
        e_total = sum(matrix_M_eg(:)) ;
        
        parfor m = 1:m_deRand
            
            % disp(m);
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


% Save results
save(['clp_binary_stg' num2str(graphon_setting) '.mat'], 'res_example');

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



