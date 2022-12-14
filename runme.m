%
% Silent slicewise shimming
% Author: Ali aghaeifar (ali.aghaeifar@tuebingen.mpg.de) 
%

clear;
clc;

% parameters
param.z_scalefactor = 3;
param.reslc_pre     = 'r_';
param.max_current   = 3;
param.lambda        = [0 1e-1 1 1e1 1e2 1e3 1e4 1e5 1e6]; % regularization weightings

% home
fpath.dir.home      = '/DATA/aaghaeifar/rawdata/silent_shimming';
fpath.dir.data      = fullfile(fpath.dir.home, 'data');
fpath.dir.results   = fullfile(fpath.dir.home, 'results');
% where to keep SODA output
fpath.dir.soda      = fullfile(fpath.dir.data, 'soda_nii');
% where to keep resliced and transformed data to standard space
fpath.dir.resliced  = fullfile(fpath.dir.data, 'resliced_nii');

% dual-echo B0 mapping file in ismrmrd format  
fpath.file.b0_mrd   = fullfile(fpath.dir.data, 'mrd', 'meas_MID00575_FID29449_aa_B0Phantom.mrd');
% slice orientation data (SODA), exported from protocol
fpath.file.soda_txt = fullfile(fpath.dir.data, 'soda_txt', 'manual.txt');   
% b0 map
fpath.file.phase_nii  = fullfile(fpath.dir.data, 'b0.nii');   % unwrapped and in Hz
fpath.file.mag_nii    = fullfile(fpath.dir.data, 'mag.nii');  %
fpath.file.mask_nii   = fullfile(fpath.dir.data, 'mask.nii'); % binary mask
% standard space; all files will be transformed to the standard space
fpath.file.std_space_nii = fullfile(fpath.dir.data, 'std_space.nii');
% basis maps of shims
fpath.file.shimsmap_nii  = fullfile(fpath.dir.data, 'shimsmap.nii'); 

addpath(fullfile(fpath.dir.home, 'src'));

if isfolder(fpath.dir.results) == false
    mkdir(fpath.dir.results);
end

%% reconstruct & prepare b0
% prepare_B0(fpath.file.b0_mrd, fpath.file)

%% load & prepare adjustment
prepare_SODA(fpath.file.soda_txt, fpath.dir.soda, param.z_scalefactor);

%% reslice shimpsmap to standard space 
[filepath, name, ext] = fileparts(fpath.file.shimsmap_nii);
clear flist;
flist{1} = fpath.file.std_space_nii;
flist{2} = spm_select('ExtFPList', filepath, [name ext], inf);
spm_reslice(char(flist), struct('mean', false, 'which', 1, 'prefix', param.reslc_pre , 'mask', false, 'interp', 2));
movefile(fullfile(filepath, [param.reslc_pre name ext]), fpath.dir.resliced);

clear flist filepath name ext

%% reslice B0 map to standard space 

clear flist;
flist{1} = fpath.file.std_space_nii; % first file is destination space
flist{2} = fpath.file.phase_nii;
flist{3} = fpath.file.mag_nii;
flist{4} = fpath.file.mask_nii;
spm_reslice(char(flist), struct('mean', false, 'which', 1, 'prefix', param.reslc_pre , 'mask', false, 'interp', 0)); % nearest neighbour 

% move resliced files to a separate folder
for i=2:4
    [d, f, e] = fileparts(flist{i});
    movefile(fullfile(d, [param.reslc_pre f e]), fpath.dir.resliced);
end

%% reslice b0 map and SODA to standard space
delete(fullfile(fpath.dir.resliced, [param.reslc_pre  'scaled_slc*.nii']));
delete(fullfile(fpath.dir.resliced, [param.reslc_pre  'slc*.nii']));

clear flist;
flist{1} = fpath.file.std_space_nii; % first file is destination space
flist{2} = spm_select('FPList', fpath.dir.soda, '^scaled_slc.*.nii$');  % scaled SODA
flist{3} = spm_select('FPList', fpath.dir.soda, '^slc.*.nii$');         % original SODA
spm_reslice(char(flist), struct('mean', false, 'which', 1, 'prefix', param.reslc_pre , 'mask', false, 'interp', 0));

% move resliced files to a separate folder
movefile(fullfile(fpath.dir.soda, [param.reslc_pre '*.nii']), fpath.dir.resliced);

% combine individual SODA for demonestration
if size(flist{2}, 1) > 1
    Vi = spm_vol(spm_select('FPList', fpath.dir.resliced, ['^' param.reslc_pre 'scaled_slc.*.nii$']));
    spm_imcalc(Vi, fullfile(fpath.dir.data, ['soda_all_zscale=' num2str(param.z_scalefactor) '.nii']), 'sum(X)', {1, -1, -7,  16}); % {'dmtx', 'mask', 'interp', 'dtype'} = {read to matrix, nan should be zeroed, 7th degree sinc interpolation, float output}
    
    Vi = spm_vol(spm_select('FPList', fpath.dir.resliced, ['^' param.reslc_pre 'slc.*.nii$']));
    spm_imcalc(Vi, fullfile(fpath.dir.data, 'soda_all.nii'), 'sum(X)', {1, -1, -7,  16}); % {'dmtx', 'mask', 'interp', 'dtype'} = {read to matrix, nan should be zeroed, 7th degree sinc interpolation, float output}
else
    copyfile(flist{2}, fullfile(fpath.dir.data, ['soda_all_zscale=' num2str(param.z_scalefactor) '.nii']));
    copyfile(flist{3}, fullfile(fpath.dir.data, 'soda_all.nii'));
end

clear flist Vi i name filepath ext ans

%% It is the time to do it! 
% load maps
[~, f, e] = fileparts(fpath.file.phase_nii);
b0        = niftiread(fullfile(fpath.dir.resliced, [param.reslc_pre f e])); b0 = double(b0);
[~, f, e] = fileparts(fpath.file.shimsmap_nii);
shimsmap  = niftiread(fullfile(fpath.dir.resliced, [param.reslc_pre f e])); shimsmap = double(shimsmap);
[~, f, e] = fileparts(fpath.file.mask_nii);
mask      = niftiread(fullfile(fpath.dir.resliced, [param.reslc_pre f e]));
mask      = imbinarize(mask, 0.5); % due to reslice, there are some non binary values. 

soda_lsit        = spm_select('FPList', fpath.dir.resliced, ['^' param.reslc_pre 'slc.*.nii$']);
scaled_soda_list = spm_select('FPList', fpath.dir.resliced, ['^' param.reslc_pre 'scaled_slc.*.nii$']);

% current boundries
current_ul = param.max_current * ones(size(shimsmap, 4), 1);
current_ll = -current_ul;


n_slice      = size(scaled_soda_list, 1);
coef         = zeros(n_slice, size(shimsmap, 4), numel(param.lambda)); % shim currents with regularized coefficient 
sd           = zeros(n_slice, numel(param.lambda) + 1);
shimsmap_row = reshape(shimsmap, [], size(shimsmap, 4));

% delete(findall(0)); % close if there is any figure open
f = waitbar(0,'Please wait...');
% loop over soda and calculate shims
for cslc=1:n_slice
% for cslc=1:5:n_slice    
    soda      = niftiread(strtrim(scaled_soda_list(cslc,:)));
    mask_comb = mask .* soda;
    mask_comb(mask_comb <= 0) = nan;
    
    for i = 1:numel(param.lambda)
        lambda = param.lambda(i);
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

%% Calculate shimmed field and std of residual inhomogeneities
% delete(findall(0)); % close if there is any figure open

[~, f, e] = fileparts(fpath.file.phase_nii);
hdr = niftiinfo(fullfile(fpath.dir.resliced, [param.reslc_pre f e]));
f   = waitbar(0,'Please wait...');
for i = 1:numel(param.lambda)
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

    niftiwrite(shimmed, fullfile(fpath.dir.results, ['shimmed_lambda=' num2str(param.lambda(i)) '_zscale=' num2str(param.z_scalefactor) '.nii']), hdr)
end
close(f)

save(fullfile(fpath.dir.results, ['vars_zscale=' num2str(param.z_scalefactor)]), 'coef', 'param', 'sd');
clearvars -except coef sd













