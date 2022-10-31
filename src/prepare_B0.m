function prepare_B0(filename, output_names)

erode_size    = 5;

myMRD   = recoMRD_B0(filename);

% mask the brain 
mask    = bet2_mex_interface(myMRD.img_mag);
% erode mask to remove surface noisy voxels
erd_krn = ones((size(mask) > 1) * (erode_size-1) + 1);
mask_e  = imerode(mask, erd_krn);

% b0 map
b0      = myMRD.img_b0 ;
% unwrap phase
b0_uw   = UnWrap_mex_interface(b0) .* mask_e; % unwrapping fails if I mask b0 in advance
% convert radian to Hz
b0_hz   = myMRD.get_b0hz(b0_uw);


% store results as nifti
myMRD.make_nifti(output_names.phase_nii, b0_hz);
myMRD.make_nifti(output_names.mag_nii, myMRD.img_mag);
myMRD.make_nifti(output_names.mask_nii, mask_e);

disp('Done.');


% check niftis
% mrd_files = {output_names.phase_nii, output_names.b0_nii, output_names.mask_nii};
% spm_check_registration(char(mrd_files));
