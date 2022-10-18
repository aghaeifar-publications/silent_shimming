%
% Silent slicewise shimming
% Author: Ali aghaeifar (ali.aghaeifar@tuebingen.mpg.de) 
%

clear;
clc;

% parameters
param.mask_erode_size      = 5;
param.tickness_scalefactor = 3.0;

% paths 
fpath.b0_mrd        = '/DATA/aaghaeifar/rawdata/silent_shimming/data/meas_MID00504_FID00909_gre_S_T_3D.mrd';
fpath.b0_mrd        = '/DATA/aaghaeifar/rawdata/b0phantom/meas_MID00575_FID29449_aa_B0Phantom.mrd';

fpath.soda_txt      = '/DATA/aaghaeifar/rawdata/silent_shimming/data/gre_s_t.txt';   
fpath.soda_txt     = '/DATA/aaghaeifar/rawdata/silent_shimming/data/manual.txt';   

fpath.home_dir      = '/DATA/aaghaeifar/rawdata/silent_shimming/';
fpath.soda_dir      = fullfile(fpath.home_dir, 'soda_nii');
fpath.resliced_dir  = fullfile(fpath.home_dir, 'resliced_nii');
fpath.b0_nii        = fullfile(fpath.home_dir, 'b0.nii');
fpath.mag_nii       = fullfile(fpath.home_dir, 'mag.nii');
fpath.mask_nii      = fullfile(fpath.home_dir, 'mask.nii');
fpath.std_space_nii = fullfile(fpath.home_dir, 'std_space.nii');
fpath.shimsmap_nii  = fullfile(fpath.home_dir, 'shimsmap.nii'); 

%% reconstruct & prepare b0
myMRD = recoMRD_B0(fpath.b0_mrd);

b0      = myMRD.img_b0;
% unwrap phase
b0_uw   = UnWrap_mex_interface(b0);
% convert radian to Hz
b0_hz   = myMRD.get_b0hz(b0_uw);
% mask the brain 
mask    = bet2_mex_interface(myMRD.img_mag);
% erode mask to remove surface noisy voxels
erd_krn = ones((size(mask) > 1) * (param.mask_erode_size-1) + 1);
mask_e  = imerode(mask, erd_krn);

% store results as nifti
myMRD.make_nifti(fpath.b0_nii, b0_hz);
myMRD.make_nifti(fpath.mag_nii, myMRD.img_mag);
myMRD.make_nifti(fpath.mask_nii, mask_e);
% check niftis
% mrd_files = {fpath.mag_nii, fpath.b0_nii, fpath.mask_nii};
% spm_check_registration(char(mrd_files));
clear mask_e erd_krn mask b0_hz b0_uw b0

%% load & prepare adjustment
myADJ = adjustments_obj(fpath.soda_txt);

% cleaning output directories
delete(fullfile(fpath.soda_dir, '*.nii'));
delete(fullfile(fpath.resliced_dir, 'r_sslc*.nii'));
delete(fullfile(fpath.resliced_dir, 'r_slc*.nii'));

for i=1:numel(myADJ.slc_soda.T)
    FoV          = myADJ.slc_soda.sorted{i}.FoV;
    FoV_scaled   = FoV;
    FoV_scaled.z = FoV.z * param.tickness_scalefactor;

    matrixSize          = myADJ.slc_soda.sorted{i}.matrixSize;
    matrixSize_scaled   = matrixSize;
    matrixSize_scaled.z = 2 * floor(matrixSize.z * param.tickness_scalefactor/2) + 1; % rounded to the nearest odd number

    % applying gussian kernel to profile in z direction to weight included adjacent areas 
    gaussian_std     = floor(matrixSize_scaled.z/2);
    smoothing_kernel = gaussmf(-floor(matrixSize_scaled.z/2):floor(matrixSize_scaled.z/2), [gaussian_std, 0]);
    volume           = ones(matrixSize.x, matrixSize.y, matrixSize.z);
    volume_scaled    = ones(matrixSize_scaled.x, matrixSize_scaled.y, matrixSize_scaled.z);
    volume_scaled    = permute(permute(volume_scaled, [3, 1, 2]) .* smoothing_kernel', [2, 3, 1]);

    ImageOrientation      = myADJ.slc_soda.T{i}(1:3, 1:2);
    ImageOrientation(:,2) = -ImageOrientation(:, 2); % experimentally discovered
    ImagePosition         = myADJ.slc_soda.T{i}(1:3, 4);


    filename         = fullfile(fpath.soda_dir, ['slc' num2str(i, '%03d') '.nii']);
    affine_transform = nii_tools.make_affine(ImageOrientation, ImagePosition, FoV, matrixSize);        
    nii_tools.create(volume, filename, affine_transform)

    filename         = fullfile(fpath.soda_dir, ['sslc' num2str(i, '%03d') '.nii']);
    affine_transform = nii_tools.make_affine(ImageOrientation, ImagePosition, FoV_scaled, matrixSize_scaled);        
    nii_tools.create(volume_scaled, filename, affine_transform)
end

% flist = spm_select('FPList', fpath.soda_dir, '^sslc.*.nii$');
% spm_check_registration(flist));
clear ImageOrientation ImagePosition smoothing_kernel gaussian_std affine_transform volume FoV matrixSize i filename volume_scaled FoV_scaled matrixSize_scaled

%% reslice to shimsmap space, approach 1, using spm
delete(fullfile(fpath.resliced_dir, 'r_sslc*.nii'));
prefix = 'r_';

clear flist;
flist{1} = [fpath.shimsmap_nii, ',1'];
flist{1} = fpath.std_space_nii;
flist{2} = fpath.b0_nii;
flist{3} = fpath.mag_nii;
flist{4} = fpath.mask_nii;
flist{5} = spm_select('ExtFPList', fpath.home_dir, '^shimsmap.*.nii$', inf);
flist{6} = spm_select('FPList', fpath.soda_dir, '^sslc.*.nii$'); % scaled SODA
flist{7} = spm_select('FPList', fpath.soda_dir, '^slc.*.nii$');  % original SODA
spm_reslice(char(flist), struct('mean', false, 'which', 1, 'prefix', prefix, 'mask', false, 'interp', 2));

% move resliced files to a separate folder
movefile(fullfile(fpath.home_dir, 'r_*.nii'), fpath.resliced_dir);
movefile(fullfile(fpath.soda_dir, 'r_*.nii'), fpath.resliced_dir);

% combine individual SODA for demonestration
if size(flist{6}, 1) > 1
    Vi = spm_vol(spm_select('FPList', fpath.resliced_dir, ['^' prefix 'sslc.*.nii$']));
    spm_imcalc(Vi, fullfile(fpath.home_dir, 'soda_scaled_all.nii'), 'sum(X)', {1, -1, -7,  16}); % {'dmtx', 'mask', 'interp', 'dtype'} = {read to matrix, nan should be zeroed, 7th degree sinc interpolation, float output}
    
    Vi = spm_vol(spm_select('FPList', fpath.resliced_dir, ['^' prefix 'slc.*.nii$']));
    spm_imcalc(Vi, fullfile(fpath.home_dir, 'soda_all.nii'), 'sum(X)', {1, -1, -7,  16}); % {'dmtx', 'mask', 'interp', 'dtype'} = {read to matrix, nan should be zeroed, 7th degree sinc interpolation, float output}
else
    copyfile(flist{6}, fullfile(fpath.home_dir, 'soda_scaled_all.nii'));
    copyfile(flist{7}, fullfile(fpath.home_dir, 'soda_all.nii'));
end

clear flist Vi

%% reslice to shimsmap space, approach 2
% delete(fullfile(fpath.resliced_dir, 'rr_*.nii'));
% prefix = 'rr_';
% 
% nii_info = getfield(niftiinfo(fpath.shimsmap_nii), 'raw');
% affine_ref = [nii_info.srow_x; nii_info.srow_y; nii_info.srow_z; [0, 0, 0, 1]];
% size_ref   = nii_info.dim(2:4);
% 
% clear flist;
% flist{1} = fpath.b0_nii;
% flist{2} = fpath.mag_nii;
% flist{3} = fpath.mask_nii;
% flist{4} = spm_select('FPList', fpath.soda_dir, '^slc.*.nii$');
% 
% flist = char(flist);
% for i=1:size(flist,1)
%     output_vol = transform_to_x_space(strtrim(flist(i,:)), affine_ref, size_ref, [], 0);
%     [filepath, name, ext] = fileparts(flist(i,:));
%     nii_tools.create(output_vol, fullfile(filepath, [prefix name ext]), affine_ref);
%     progressmsg(i/size(flist,1), 'Reslice to shimsmap: ') ;
% end
% 
% % move resliced files to a separate folder
% movefile(fullfile(fpath.home_dir, 'rr_*.nii'), fpath.resliced_dir);
% movefile(fullfile(fpath.soda_dir, 'rr_*.nii'), fpath.resliced_dir);
% 
% % go to world coordinate for demonestration
% res = 2;
% fov = 300;
% [mag_std_space, ~, affine_transform]  = transform_to_std_space(fpath.mag_nii, fov, res, 'linear', 0);
% 
% flist = spm_select('FPList', fpath.resliced_dir, ['^' prefix 'slc.*.nii$']);
% soda_std_space = zeros(fov/res, fov/res, fov/res);
% for i=1:size(flist,1)
%     soda_std_space = soda_std_space + transform_to_std_space(flist(i,:), fov, res, 'linear', 0);
%     progressmsg(i/size(flist,1), 'Transform to standard space: ') ;
% end
% 
% nii_tools.create(soda_std_space, fullfile(fpath.home_dir, 'soda_all_std.nii'), affine_transform);
% nii_tools.create(mag_std_space, fullfile(fpath.home_dir, 'mag_std.nii'), affine_transform);
% clear res fov soda_std_space mag_std_space

%% It is the time to do it! 
% load maps
[~, f, e] = fileparts(fpath.b0_nii);
hdr       = niftiinfo(fullfile(fpath.resliced_dir, [prefix, f, e]));
b0        = double(niftiread(hdr));

[~, f, e] = fileparts(fpath.mask_nii);
mask      = niftiread(fullfile(fpath.resliced_dir, [prefix, f, e]));
mask      = imbinarize(mask, 0.5); % due to reslice, there are some non binary values. 

[~, f, e] = fileparts(fpath.shimsmap_nii);
shimsmap  = double(niftiread(fullfile(fpath.resliced_dir, [prefix, f, e])));


% current boundries
up = 2 * ones(size(shimsmap, 4), 1);
ll = -up;

% loop over soda and calculate shims
coef         = zeros(numel(myADJ.slc_soda.T), size(shimsmap, 4));
shimmed      = zeros(size(b0));
sd           = zeros(numel(myADJ.slc_soda.T), 2);
shimsmap_row = reshape(shimsmap, [], size(shimsmap, 4));

for i=1:numel(myADJ.slc_soda.T)
    soda      = niftiread(fullfile(fpath.resliced_dir, [prefix 'sslc' num2str(i, '%03d') '.nii']));
    mask_comb = mask .* soda;
    mask_comb(mask_comb <= 0) = nan;
    coef(i,:) = shim_it(b0, shimsmap, mask_comb, ll, up);

    % calculated shimmed volume and std 
    soda      = niftiread(fullfile(fpath.resliced_dir, [prefix 'slc' num2str(i, '%03d') '.nii']));
    mask_comb = mask .* soda;

    temp      = mask_comb(:) .* (b0(:) + shimsmap_row * coef(i,:)');
    temp      = reshape(temp, size(b0));
    temp(isnan(temp)) = 0;
    shimmed   = shimmed + temp;
    sd(i,:)   = [std(temp .* mask_comb, [], "all", "omitnan"), ...
                 std(b0   .* mask_comb, [], "all", "omitnan")];
end
niftiwrite(shimmed, fullfile(fpath.home_dir, 'shimmed.nii'), hdr)

clear temp shimsmap mask b0 f e hdr




