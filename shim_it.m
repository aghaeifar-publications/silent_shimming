function coef = shim_it(b0map, basismaps, mask, current_ll, current_ul)

basismaps_size = size(basismaps);         
if numel(b0map) ~= prod(basismaps_size(1:end-1)) || numel(mask) ~= numel(b0map)
    error('size mismatch! aborting...');
end

% b0map_bu     = b0map;
% mask_bu      = mask;
% basismaps_bu = basismaps;

basismaps = reshape(basismaps, prod(basismaps_size(1:end-1)), basismaps_size(end));
b0map     = b0map(:);
mask      = mask(:);

% remove nan
ind             = isnan(b0map) | isnan(mask);
basismaps(ind,:)= [];
mask(ind)       = [];
b0map(ind)      = [];

if any(isnan(basismaps(:,1)), 'all')
    warning('basismapss volume does not cover whole b0map (%.2f%%).\n', 100*sum(isnan(basismaps(1,:))) / numel(b0map));
end

ind             = isnan(basismaps(:,1));
basismaps(ind,:)= [];
mask(ind)       = [];
b0map(ind)      = [];

if isempty(b0map)
    warning('Empty slice to shim.');
    coef = zeros(size(current_ll));
else
    % form L2 weighted norm
    positive_ind = mask>0; % should not be negative
    mask_w       = sqrt(mask(positive_ind));
    b0map        = double(b0map(positive_ind) .* mask_w);
    basismaps    = double(basismaps(positive_ind, :) .* mask_w);
    
    % Run optimization
    options = optimoptions('lsqlin', 'Display', 'iter-detailed', 'Algorithm', 'interior-point', 'MaxIter', 100000); % sometimes the method gives not results as well as fmincon without this options
    coef = lsqlin(basismaps, -b0map, [], [], [], [], current_ll, current_ul, [], options);
end

%% some additional outputs 
% basismaps_bu = reshape(basismaps_bu, prod(basismaps_size(1:end-1)), basismaps_size(end));
% shimmed = b0map_bu(:) + basismaps_bu * coef;
% shimmed = reshape(shimmed, size(b0map_bu));
% sd      = [std(b0map_bu .* mask_bu, [], "all", "omitnan"), ...
%            std(shimmed  .* mask_bu, [], "all", "omitnan")];