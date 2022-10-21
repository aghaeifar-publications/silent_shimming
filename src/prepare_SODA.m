function prepare_SODA(filename, output_dir, z_factor)

myADJ = adjustments_obj(filename);

% cleaning output directories
delete(fullfile(output_dir, '*.nii'));

for i=1:numel(myADJ.slc_soda.T)
    % original SODA
    FoV          = myADJ.slc_soda.sorted{i}.FoV;
    FoV_scaled   = FoV;
    matrixSize   = myADJ.slc_soda.sorted{i}.matrixSize;
    volume       = ones(matrixSize.x, matrixSize.y, matrixSize.z);

    % scaled SODA
    FoV_scaled.z        = FoV.z * z_factor;
    matrixSize_scaled   = matrixSize;
    matrixSize_scaled.z = z_factor * (matrixSize.z - 1) + 1; % rounded to the nearest odd number

    % applying gussian kernel to profile in z direction to weight included adjacent areas 
    gaussian_std     = floor(matrixSize_scaled.z/2);
    smoothing_kernel = gaussmf(-floor(matrixSize_scaled.z/2):floor(matrixSize_scaled.z/2), [gaussian_std, 0]);    
    volume_scaled    = ones(matrixSize_scaled.x, matrixSize_scaled.y, matrixSize_scaled.z);
    volume_scaled    = permute(permute(volume_scaled, [3, 1, 2]) .* smoothing_kernel', [2, 3, 1]);

    ImageOrientation      = myADJ.slc_soda.T{i}(1:3, 1:2);
    ImageOrientation(:,2) = -ImageOrientation(:, 2); % experimentally discovered
    ImagePosition         = myADJ.slc_soda.T{i}(1:3, 4);

    % save original SODA
    filename         = fullfile(output_dir, ['slc' num2str(i, '%03d') '.nii']);
    affine_transform = nii_tools.make_affine(ImageOrientation, ImagePosition, FoV, matrixSize);        
    nii_tools.create(volume, filename, affine_transform)

    % save scaled and smoothed SODA in Z direction
    filename         = fullfile(output_dir, ['scaled_slc' num2str(i, '%03d') '.nii']);
    affine_transform = nii_tools.make_affine(ImageOrientation, ImagePosition, FoV_scaled, matrixSize_scaled);        
    nii_tools.create(volume_scaled, filename, affine_transform)
end

disp('Done.');