
clear; clc;

figure(2); clf;

result_dir   = '/DATA/aaghaeifar/rawdata/silent_shimming/results';
z_factor     = [1 5 9];
slc_tickness = 3; % mm
max_tickness = max(z_factor) * slc_tickness;
plt_samples  = linspace(-max_tickness/2, max_tickness/2, 100*max_tickness);
slc_samples  = numel(plt_samples)/max_tickness * slc_tickness;

slc = [106, 125];
temp = niftiread(fullfile(result_dir, '..', 'data', 'resliced_nii', 'r_mask.nii'));
mask = temp(:,:,slc);
temp = niftiread(fullfile(result_dir, '..', 'data', 'resliced_nii', 'r_b0.nii'));
vol = temp(:,:,slc) .* mask;
iRSS = zeros(1,numel(z_factor));

for i=1:numel(z_factor)
    filename = ['vars_zscale=' num2str(z_factor(i)) '.mat'];
    load(fullfile(result_dir, filename));

    C = coef(:, 2:end, 1); % Exclude Freq. and pick no regularization case
    D = pdist2(C, C, 'euclidean');
    current_change_ascending_order = [diag(D(1:end-1, 2:end)); D(end, 1)];

    subplot(2,3,4:6); 
    plot(current_change_ascending_order, 'LineWidth', 2); ylim([0 18]); hold on;
    iRSS(i) = sum(current_change_ascending_order);

    z = z_factor(i) * slc_tickness; % 3 is number of samples per slice with no z_factor. see prepare_SODA(...)
    gaussian_std     = z/2;
    smoothing_kernel = gaussmf(plt_samples, [gaussian_std, 0]);

    cutoff_kernel    = ones(1, (slc_samples-1) * z_factor(i) + 1);
    cutoff_kernel    = padarray(cutoff_kernel, [0 (numel(smoothing_kernel)-numel(cutoff_kernel))/2], 0, 'both');

    subplot(2,3,1);
    plot(plt_samples, smoothing_kernel .* cutoff_kernel, 'LineWidth', 2); hold on; ylim([0 1])

    %%%%%%%%%%
    filename = ['shimmed_lambda=0_zscale=' num2str(z_factor(i)) '.nii'];
    temp = niftiread(fullfile(result_dir, filename));
    vol = cat(3, vol, temp(:,:,slc));
end

subplot(2,3,2:3);
vol = flip(permute(vol, [2 1 3]), 1);
vol = vol(:,:, [1:2:end, 2:2:end]); 
vol(vol == 0) = nan;
win = centerCropWindow3d(size(vol),[floor(0.7*size(vol,1)) floor(0.7*size(vol,2)) size(vol,3)]);
vol = imcrop3(vol,win);
v = montage(vol, "Size",[2 4]); v = v.CData;
imagesc(v, 'AlphaData', ~isnan(v)); axis off; axis image; clim([-200 200]); colormap(jet); colorbar('box', 'off');
std(reshape(vol, [], size(vol, 3)), [], 1, "omitnan")


subplot(2,3, 4:6); 
legend('1', '5', '9', 'Color', 'None', 'Box', 'off')
set(gca, 'Color', 'None', 'Box', 'off', 'LineWidth', 2)
xlabel('Slice Acquisition Order')


subplot(2,3, 1);
title('Gaussian profile for adjacent slices weightening')
set(gca, 'Color', 'None', 'Box', 'off', 'LineWidth', 2);
mx = 16.5;
xlim([-mx mx]);
for i=-6:6
    x = [-slc_tickness/2, -slc_tickness/2, slc_tickness/2, slc_tickness/2] + i*slc_tickness;
    y = [0 1 1 0];
    patch(x, y, rem(abs(i),2), 'LineWidth', 1.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end
%%
