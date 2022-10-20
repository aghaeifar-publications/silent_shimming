%
% Silent slicewise shimming
% Author: Ali aghaeifar (ali.aghaeifar@tuebingen.mpg.de) 
%

clear;
clc;

% parameters
param.z_scalefactor = 3;
param.reslc_pre     = 'r_';

% dual-echo B0 mapping file in ismrmrd format  
fpath.b0_mrd        = '/DATA/aaghaeifar/rawdata/silent_shimming/data/meas_MID00504_FID00909_gre_S_T_3D.mrd';
fpath.b0_mrd        = '/DATA/aaghaeifar/rawdata/b0phantom/meas_MID00575_FID29449_aa_B0Phantom.mrd';

% slice orientation data (SODA), exported from protocol
fpath.soda_txt      = '/DATA/aaghaeifar/rawdata/silent_shimming/data/gre_s_t.txt';   
fpath.soda_txt      = '/DATA/aaghaeifar/rawdata/silent_shimming/data/manual.txt';   

% home
fpath.home_dir      = '/DATA/aaghaeifar/rawdata/silent_shimming/';
% where to keep SODA output
fpath.soda_dir      = fullfile(fpath.home_dir, 'soda_nii');
% where to keep resliced and transformed data to standard space
fpath.resliced_dir  = fullfile(fpath.home_dir, 'resliced_nii');
% b0 map
fpath.b0.phase_nii  = fullfile(fpath.home_dir, 'b0.nii');   % unwrapped and in Hz
fpath.b0.mag_nii    = fullfile(fpath.home_dir, 'mag.nii');  %
fpath.b0.mask_nii   = fullfile(fpath.home_dir, 'mask.nii'); % binary mask
% standard space; all files will be transformed to the standard space
fpath.std_space_nii = fullfile(fpath.home_dir, 'std_space.nii');
% basis maps of shims
fpath.shimsmap_nii  = fullfile(fpath.home_dir, 'shimsmap.nii'); 

%% reconstruct & prepare b0
prepare_B0(fpath.b0_mrd, fpath.b0)

%% load & prepare adjustment
prepare_SODA(fpath.soda_txt, fpath.soda_dir, param.z_scalefactor);

%% reslice to standard space -> shimpsmap
[filepath, name, ext] = fileparts(fpath.shimsmap_nii);
flist = cell(2,1);
flist{1} = fpath.std_space_nii;
flist{2} = spm_select('ExtFPList', filepath, [name ext], inf);
spm_reslice(char(flist), struct('mean', false, 'which', 1, 'prefix', param.reslc_pre , 'mask', false, 'interp', 2));
movefile(fullfile(filepath, [param.reslc_pre name ext]), fpath.resliced_dir);

clear flist filepath name ext

%% reslice to standard space -> b0 map and SODA
delete(fullfile(fpath.resliced_dir, [param.reslc_pre  'scaled_slc*.nii']));
delete(fullfile(fpath.resliced_dir, [param.reslc_pre  'slc*.nii']));

flist = cell(6,1);
flist{1} = fpath.std_space_nii;
flist{2} = fpath.b0.phase_nii;
flist{3} = fpath.b0.mag_nii;
flist{4} = fpath.b0.mask_nii;
flist{5} = spm_select('FPList', fpath.soda_dir, '^scaled_slc.*.nii$'); % scaled SODA
flist{6} = spm_select('FPList', fpath.soda_dir, '^slc.*.nii$');  % original SODA
spm_reslice(char(flist), struct('mean', false, 'which', 1, 'prefix', param.reslc_pre , 'mask', false, 'interp', 2));

% move resliced files to a separate folder
for i=2:4
    [filepath, name, ext] = fileparts(flist{i});
    movefile(fullfile(filepath, [param.reslc_pre name ext]), fpath.resliced_dir);
end

movefile(fullfile(fpath.soda_dir, [param.reslc_pre '*.nii']), fpath.resliced_dir);

% combine individual SODA for demonestration
if size(flist{5}, 1) > 1
    Vi = spm_vol(spm_select('FPList', fpath.resliced_dir, ['^' param.reslc_pre 'scaled_slc.*.nii$']));
    spm_imcalc(Vi, fullfile(fpath.home_dir, 'soda_scaled_all.nii'), 'sum(X)', {1, -1, -7,  16}); % {'dmtx', 'mask', 'interp', 'dtype'} = {read to matrix, nan should be zeroed, 7th degree sinc interpolation, float output}
    
    Vi = spm_vol(spm_select('FPList', fpath.resliced_dir, ['^' param.reslc_pre 'slc.*.nii$']));
    spm_imcalc(Vi, fullfile(fpath.home_dir, 'soda_all.nii'), 'sum(X)', {1, -1, -7,  16}); % {'dmtx', 'mask', 'interp', 'dtype'} = {read to matrix, nan should be zeroed, 7th degree sinc interpolation, float output}
else
    copyfile(flist{5}, fullfile(fpath.home_dir, 'soda_scaled_all.nii'));
    copyfile(flist{6}, fullfile(fpath.home_dir, 'soda_all.nii'));
end

clear flist Vi i name filepath ext ans

%% It is the time to do it! 
% load maps
[~, f, e] = fileparts(fpath.b0.phase_nii);
hdr       = niftiinfo(fullfile(fpath.resliced_dir, [param.reslc_pre, f, e]));
b0        = double(niftiread(hdr));

[~, f, e] = fileparts(fpath.b0.mask_nii);
mask      = niftiread(fullfile(fpath.resliced_dir, [param.reslc_pre, f, e]));
mask      = imbinarize(mask, 0.5); % due to reslice, there are some non binary values. 

[~, f, e] = fileparts(fpath.shimsmap_nii);
shimsmap  = double(niftiread(fullfile(fpath.resliced_dir, [param.reslc_pre, f, e])));

soda_lsit        = spm_select('FPList', fpath.resliced_dir, ['^' param.reslc_pre 'slc.*.nii$']);
scaled_soda_list = spm_select('FPList', fpath.resliced_dir, ['^' param.reslc_pre 'scaled_slc.*.nii$']);

% current boundries
current_ul = 3 * ones(size(shimsmap, 4), 1);
current_ll = -current_ul;

% regularization weightings
lambda_all   = [0, 1e1 1e2 1e3 1e4 1e5 1e6];

n_slice      = size(scaled_soda_list, 1);
coef         = zeros(n_slice, size(shimsmap, 4), numel(lambda_all)); % shim currents with regularized coefficient 
sd           = zeros(n_slice, numel(lambda_all) + 1);
shimsmap_row = reshape(shimsmap, [], size(shimsmap, 4));

% delete(findall(0)); % close if there is any figure open
f = waitbar(0,'Please wait...');
% loop over soda and calculate shims
for cslc=1:n_slice
% for cslc=1:5:n_slice    
    soda      = niftiread(strtrim(scaled_soda_list(cslc,:)));
    mask_comb = mask .* soda;
    mask_comb(mask_comb <= 0) = nan;
    
    for i = 1:numel(lambda_all)
        lambda = lambda_all(i);
        if cslc == 1 % there is no previous slice for cslc=1, thus, reset lambda
            lambda = 0;
            old_current = zeros(size(current_ul));
        else            
            old_current = coef(cslc-1,:,i)';
        end

        coef(cslc,:,i) = rclsq(b0, shimsmap, mask_comb, current_ll, current_ul, lambda, old_current);
    end
    waitbar(cslc/n_slice, f, 'Optimization is running...');
end
close(f)

clear f e current_ul current_ll scaled_soda_list lambda cslc i soda mask_comb old_current shimsmap f

%% Calculate shimmed field and std of residual inhomogeneities
% delete(findall(0)); % close if there is any figure open
f = waitbar(0,'Please wait...');
for i = 1:numel(lambda_all)
    shimmed = zeros(size(b0));
    for cslc=1:n_slice
        % calculated shimmed volume and std 
        soda        = niftiread(strtrim(soda_lsit(cslc,:)));
        mask_comb   = mask .* soda;
        sd(cslc,1)  = std(b0 .* mask_comb, [], "all", "omitnan");
    
        temp        = mask_comb(:) .* (b0(:) + shimsmap_row * coef(cslc,:,i)');
        temp        = reshape(temp, size(b0));
        temp(isnan(temp)) = 0;
        shimmed     = shimmed + temp;
        sd(cslc, i+1)  = std(temp .* mask_comb, [], "all", "omitnan");

        waitbar(cslc/n_slice, f, 'Running...');
    end

    niftiwrite(shimmed, fullfile(fpath.home_dir, ['shimmed_lambda=' num2str(lambda_all(i)) '.nii']), hdr)
end
close(f)

clear temp mask b0 hdr shimmed soda_lsit i cslc n_slice mask_comb shimsmap_row soda lambda f

%% Method 1, sort based on euclidean distance
C = coef(:,2:end,:); % exclude Freq.
D = zeros(size(C, 1), size(C, 1), size(C, 3));
for i=1:size(C, 3)
    D(:,:,i)  = pdist2(C(:,:,i), C(:,:,i));
end

montage(D)
clim([1.2*min(D(:)), 0.8*max(D(:))]);
colormap jet

%% a) find 
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

D_not_ordered = [diag(D(1:end-1,2:end)); D(end, 1)];

plot(D_ordered); hold on;
plot(D_not_ordered); hold off;
legend('Ordered', 'Not Ordered')

%% b) assignment problem
% https://www.mathworks.com/matlabcentral/answers/302378-how-to-select-one-value-from-each-column-and-one-value-from-each-row-and-get-minimal-sum
% https://stackoverflow.com/questions/14034594/choose-one-element-from-each-row-and-column-in-matrix-and-the-sum-is-minimized
D_temp = D;
D_temp(1:1+size(D_temp,1):end) = 2 * max(D(:));

[assignment,cost] = munkres(D_temp);

D_ordered_Hungarian = zeros(size(D,1), 1);
for i=1:numel(D_ordered_Hungarian)
    D_ordered_Hungarian(i) = D(i, assignment(i));
end

[sum(D_ordered), sum(D_not_ordered), sum(D_ordered_Hungarian)]

plot(D_ordered); hold on;
% plot(D_not_ordered); hold on;
plot(D_ordered_Hungarian); hold off;
legend('Ordered', 'Hungarian');

%% 












