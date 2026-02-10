clear;

% Parameters
graphon_setting = 1;
node_num = 200;
missing_prob = 0.2;
train_proportion = 0.4;
r0 = 25;
reps = 100;

mean_pvals = nan(reps, node_num);

parfor rep_idx = 1:reps

    disp(rep_idx);

    mean_pvals_rep = nan(1, node_num);

    [matrix_A_eg, matrix_M_eg, ~] = ...
        gen_data_unic(graphon_setting, node_num, 1, missing_prob);
    
    % Nulls
    matrix_C_eg = matrix_A_eg;
    row_indices = find(any(matrix_M_eg == 1, 2));

    for i0 = row_indices'

        [vector_train_idx, vector_calib_idx, vector_test_idx] = ...
            gen_partition(node_num, i0, matrix_M_eg, train_proportion);

        r1 = floor(sum(vector_calib_idx) / r0);
        if r1 < 1
            continue;
        end

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

        pvals_row = [];
        for ii = 1:r2

            one_split_test = find(splits_test{ii} == 1);
            if isempty(one_split_test)
                continue;
            end

             num_parti = sum(splits_test{ii}); 
            splits_calib_j0 = cell(1, num_parti);

            for k = 1:num_parti
                splits_calib_j0{k} = zeros(1, length(vector_calib_idx));
            end

            ones_idx_calib = find(vector_calib_idx == 1);
            ones_idx_shuffled = ones_idx_calib(randperm(length(ones_idx_calib)));

            for k = 1:length(ones_idx_shuffled)
                split_idx = mod(k-1, num_parti) + 1;
                splits_calib_j0{split_idx}(ones_idx_shuffled(k)) = 1;
            end

            for k = 1:num_parti

                j0 = one_split_test(k);
                calib_idx_j0 = splits_calib_j0{k};
                calib_idx_set = find(calib_idx_j0 == 1);

                est_res_j0 = ests(node_num, i0, matrix_A_eg, matrix_M_eg, ...
                                  vector_train_idx, j0, calib_idx_set);

                pval = pvalue_single(i0, j0, matrix_A_eg, matrix_C_eg, ...
                                     est_res_j0, calib_idx_j0);

                pvals_row(end+1) = pval; %#ok<SAGROW>
            end
        end

        if ~isempty(pvals_row)
            mean_pvals_rep(i0) = mean(pvals_row);
        end
    end
    mean_pvals(rep_idx, :) = mean_pvals_rep;
end


% Save the mean p-values to a file for further analysis
save('mean_pvals.mat', 'mean_pvals');

% Load the mean p-values from the saved file for further analysis
loaded_data = load('mean_pvals.mat');
mean_pvals = loaded_data.mean_pvals;


% Heatmap
valid_rows = all(~isnan(mean_pvals), 1);
mean_pvals_clean = mean_pvals(:, valid_rows);
corr_mat = corr(mean_pvals_clean); 
corr_plot = corr_mat;

figure;
imagesc(corr_plot);
colorbar;
caxis([-0.32 0.32]);
xlabel('Row index');
ylabel('Row index');
title('Correlation matrix');

n = 256;
darkred   = [0.55, 0.00, 0.00];
white     = [1.00, 1.00, 1.00];
steelblue = [0.27, 0.51, 0.71];
cmap1 = [ ...
    linspace(darkred(1), white(1), n/2)', ...
    linspace(darkred(2), white(2), n/2)', ...
    linspace(darkred(3), white(3), n/2)'  ];
cmap2 = [ ...
    linspace(white(1), steelblue(1), n/2)', ...
    linspace(white(2), steelblue(2), n/2)', ...
    linspace(white(3), steelblue(3), n/2)' ];
cmap = [cmap1; cmap2];
colormap(cmap);


% Subset of heatmap
valid_rows = all(~isnan(mean_pvals), 1);
mean_pvals_clean = mean_pvals(:, valid_rows);

corr_mat = corr(mean_pvals_clean);
num_show = 10;
row_idx = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200];
corr_sub = corr_mat(row_idx, row_idx);
figure;
imagesc(corr_sub);
axis square;
colorbar;

max_abs = 0.32;
caxis([-max_abs max_abs]);

set(gca, 'Color', [0.9 0.9 0.9]);
n = 256;
darkred   = [0.55, 0.00, 0.00];
white     = [1.00, 1.00, 1.00];
steelblue = [0.27, 0.51, 0.71];

cmap1 = [ ...
    linspace(darkred(1), white(1), n/2)', ...
    linspace(darkred(2), white(2), n/2)', ...
    linspace(darkred(3), white(3), n/2)' ];

cmap2 = [ ...
    linspace(white(1), steelblue(1), n/2)', ...
    linspace(white(2), steelblue(2), n/2)', ...
    linspace(white(3), steelblue(3), n/2)' ];

colormap([cmap1; cmap2]);

xticks(1:num_show);
yticks(1:num_show);
xticklabels(row_idx);
yticklabels(row_idx);

xlabel('Selected row index');
ylabel('Selected row index');
title('Correlation matrix of 10 rows');  

for i = 1:num_show
    for j = 1:num_show
        if ~isnan(corr_sub(i,j))
            if abs(corr_sub(i,j)) > 0.5 * max_abs
                txt_color = 'w';   
            else
                txt_color = 'k';  
            end

            text(j, i, sprintf('%.2f', corr_sub(i,j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'Color', txt_color);
        end
    end
end

