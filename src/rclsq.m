function coef = rclsq(b0map, basismaps, mask, current_ll, current_ul, lambda, old_current)
%
% regularized constrained least square optimization for problems of the form:
%   min || A*x + b ||^2 + λ * ||x - c||^2
% there is already a solver to solve this problem, but I couldn't manage to compile it, maybe you can: https://www.cs.ubc.ca/~mpf/bcls/index.html
%

basismaps_size = size(basismaps);         
if numel(b0map) ~= prod(basismaps_size(1:end-1)) || numel(mask) ~= numel(b0map)
    error('size mismatch! aborting...');
end

current_ll= current_ll(:);
current_ul= current_ul(:);

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
    options = optimoptions('fmincon', 'UseParallel', true, 'SpecifyObjectiveGradient', true, 'Display', 'iter-detailed', 'Algorithm', 'interior-point', 'MaxFunEvals', 10^6, 'MaxIter', 100000);
    coef = fmincon(@(x) cost_fun(x, b0map, basismaps, lambda, old_current), zeros(size(current_ll)), [], [], [], [], current_ll, current_ul, [], options);
    coef = transpose(coef); 
end

end


function [err, grad] = cost_fun(x, b0map, basismaps, lambda, old_current)
    err = norm(basismaps * x + b0map).^2 + lambda * norm(x - old_current).^2;
    if nargout > 1 % gradient required
        grad = 2 * transpose(basismaps) * (b0map + basismaps*x) + 2 * lambda * (x - old_current);
    end
end