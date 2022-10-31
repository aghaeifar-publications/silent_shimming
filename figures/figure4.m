
clear; clc;
result_dir   = '/DATA/aaghaeifar/rawdata/silent_shimming/results';
load(fullfile(result_dir, 'vars_zscale=1.mat'));

C = coef(:,2:end,:); % exclude Freq.
D = pdist2(C(:,:,1), C(:,:,1), 'euclidean');

% sort based on euclidean distance, 
D_ordered       = inf(size(D,1), 1);
ind_slc_ordered = zeros(size(D,1)+1, 1);
for k=1:size(D,1)
    ind_ordered     = zeros(size(D,1)+1, 1);
    ind_ordered(1)  = k;
    ind_ordered(end)= ind_ordered(1);

    D_temp          = D;
    D_temp(1:1+size(D_temp,1):end) = inf; % set diagonal elements to inf
    D_temp(:, ind_ordered(1)) = inf; % the starting index is alreaded used and can't be chosen again

    for i = 1:numel(ind_ordered)-2
        [~, I] = sort(D_temp(ind_ordered(i), :));
        ind_ordered(i+1) = I(1);
        D_temp(:, I(1)) = inf;
%         imagesc(D_temp); drawnow; pause(0.3);
    end
    
    D_ordered_temp = zeros(size(D,1), 1);
    for i=1:numel(D_ordered_temp)
        D_ordered_temp(i) = D(ind_ordered(i), ind_ordered(i+1));
    end

    if sum(D_ordered_temp) < sum(D_ordered)
        D_ordered = D_ordered_temp;
        ind_slc_ordered = ind_ordered;
    end
end

% assignment problem
% https://www.mathworks.com/matlabcentral/answers/302378-how-to-select-one-value-from-each-column-and-one-value-from-each-row-and-get-minimal-sum
% https://stackoverflow.com/questions/14034594/choose-one-element-from-each-row-and-column-in-matrix-and-the-sum-is-minimized
D_temp = D;
D_temp(1:1+size(D_temp,1):end) = 2 * max(D(:));
[assignment, cost] = munkres(D_temp);

D_ordered_Hungarian = zeros(size(D,1), 1);
for i=1:numel(D_ordered_Hungarian)
    D_ordered_Hungarian(i) = D(i, assignment(i));
end




% plots
D_not_ordered = [diag(D(1:end-1,2:end)); D(end, 1)];

figure(4); clf
plot(D_not_ordered, 'LineWidth', 2); hold on;
plot(D_ordered, '--', 'LineWidth', 2); 
plot(D_ordered_Hungarian, '--*', 'LineWidth', 2);
set(gca, 'Color', 'None', 'Box', 'off', 'LineWidth', 2)

legend('Sequential', 'Basic Sorting', 'Hungarian Algorithm', 'Color', 'None', 'Box', 'off');


[sum(D_not_ordered), sum(D_ordered), sum(D_ordered_Hungarian)]

index_slice_all_methods = [1:size(D,1); ...
                            ind_slc_ordered(1:end-1)'; ...
                           assignment ];
disp(index_slice_all_methods)

