
clear; clc;
result_dir   = '/DATA/aaghaeifar/rawdata/silent_shimming/results';
load(fullfile(result_dir, 'vars_zscale=1.mat'));

C = coef(:,2:end,:); % exclude Freq.
D = zeros(size(C, 1), size(C, 1), size(C, 3));
for i=1:size(C, 3) % size(C, 3) = number of regularizations
    D(:,:,i)  = pdist2(C(:,:,i), C(:,:,i), 'euclidean');
end

figure(3); clf;
line_style = {'', '--', '--o', '--*'};
ind_reg = [1 5 7 9];
for i=1:numel(ind_reg)
    subplot(2,4, 1:4);
    D1 = D(:,:,ind_reg(i));
    ascending_current = [diag(D1(1:end-1, 2:end)); D1(end, 1)];
    plot(ascending_current, line_style{i}); hold on;

    subplot(2,4, 4+i);
    imagesc(D1); colormap(inferno); axis off; axis image;  clim([0 18])

end
subplot(2,4, 1:4);
legend('0', '1e2', '1e4', '1e6')