
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
iRSS = zeros(1,size(D, 3));

for i=1:size(D, 3)
    subplot(2,4, 1:4);
    D1 = D(:,:,i);
    ascending_current = [diag(D1(1:end-1, 2:end)); D1(end, 1)];
    plot(ascending_current, line_style{i}, 'LineWidth', 2); hold on;
    iRSS(i) = sum(ascending_current);

    subplot(2,4, 4+i);
    imagesc(D1); colormap(inferno); axis off; axis image;  clim([0 18])

end

subplot(2,4, 1:4);
xlabel('Slice Acquisition Order');
ylabel('RSS of Inter-Slice Shims Current Change');
ylim([0 18]);
xlim([0 43])
legend('0', '1e2', '1e4', '1e6', 'Color', 'None', 'Box', 'off')
set(gca, 'Color', 'None', 'Box', 'off', 'LineWidth', 2)

sd_wholebrain = zeros(1, 4);
file_names{1} = fullfile(result_dir, '..', 'data', 'resliced_nii', 'r_b0.nii');
file_names{2} = fullfile(result_dir, 'shimmed_lambda=0_zscale=1.nii');
file_names{3} = fullfile(result_dir, 'shimmed_lambda=100_zscale=1.nii');
file_names{4} = fullfile(result_dir, 'shimmed_lambda=10000_zscale=1.nii');
file_names{5} = fullfile(result_dir, 'shimmed_lambda=1000000_zscale=1.nii');

mask = niftiread(fullfile(result_dir, '..', 'data', 'resliced_nii', 'r_mask.nii'));

for i=1:numel(file_names)
    b0 = niftiread(file_names{i});
    temp = b0 .* mask;
    temp(temp == 0) = nan;
    sd_wholebrain(i) = std(temp, [], "all", "omitnan");
end

disp(sd_wholebrain)