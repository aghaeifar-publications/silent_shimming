clear;

% import result of step1
load('./results/vars_zscale=3.mat');

%% 
C = coef(:,2:end,:); % exclude Freq.
D = zeros(size(C, 1), size(C, 1), size(C, 3));
for i=1:size(C, 3) % size(C, 3) = number of regularizations
    D(:,:,i)  = pdist2(C(:,:,i), C(:,:,i), 'euclidean');
end

% figure 1 in abstract
D1 = D(:,:,1); % index 1 points to no regularization (lambda = 0)
% Sequential ordering : [1 2 3 4 5 6 ... end 1 ... ]
ascending_current = [diag(D1(1:end-1, 2:end)); D1(end, 1)];
% Interleaved ordering [2 4 6 8 ... ind_end_e 1 3 5 7 ... ind_end_o 2 ...]
ind_end_e = max(2:2:size(D1,1));
ind_end_o = max(1:2:size(D1,1));
interleaved_current = diag(D1(1:end-2, 3:end));
interleaved_current = [interleaved_current(2:2:end); D1(ind_end_e, 1); interleaved_current(1:2:end); D1(ind_end_o, 2)];

figure(1); clf;
subplot(2,4,[1,2,5,6]);
imagesc(D1); colormap(inferno); axis off; colorbar; axis image; 

subplot(2,4,[3,4]); 
plot(ascending_current); ylim([0 18]); title(['Sequential Acq., \Sigma Force = ' num2str(sum(ascending_current))]);

subplot(2,4,[7,8]); 
plot(interleaved_current); ylim([0 18]); title(['Interleaved Acq., \Sigma Force = ' num2str(sum(interleaved_current))]);

%% Method 1, sort based on euclidean distance, 
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
