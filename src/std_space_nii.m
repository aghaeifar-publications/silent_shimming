% create standard space nifti file

res = 1.5;
FoV = 255;

% define a nifti header for standard space 
hdr_std.sform_code= 1;
hdr_std.qform_code= 1;
hdr_std.quatern_b = 0;
hdr_std.quatern_c = 0;
hdr_std.quatern_d = 0;
hdr_std.qoffset_x = -FoV/2;
hdr_std.qoffset_y = -FoV/2;
hdr_std.qoffset_z = -FoV/2;
hdr_std.srow_x    = [res, 0, 0, hdr_std.qoffset_x];
hdr_std.srow_y    = [0, res, 0, hdr_std.qoffset_y];
hdr_std.srow_z    = [0, 0, res, hdr_std.qoffset_z];
hdr_std.dim       = [3, FoV/res, FoV/res, FoV/res, 1, 1, 1, 1]; 
hdr_std.pixdim    = [1.0, res, res, res, 1.0, 1.0, 1.0, 1.0];

size_std   = single(hdr_std.dim(2:4));
affine_std = [hdr_std.srow_x; hdr_std.srow_y; hdr_std.srow_z; 0 0 0 1];

nii_tools.create(ones(hdr_std.dim(2:4)), fullfile('/DATA/aaghaeifar/rawdata/silent_shimming/', 'std_space.nii'), affine_std);