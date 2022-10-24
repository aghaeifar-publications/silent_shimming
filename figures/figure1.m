clear;

% import result of step1
load('../results/vars_zscale=1.mat');

C = coef(:,2:end,:); % exclude Freq.
D = zeros(size(C, 1), size(C, 1), size(C, 3));
for i=1:size(C, 3) % size(C, 3) = number of regularizations
    D(:,:,i)  = pdist2(C(:,:,i), C(:,:,i), 'euclidean');
end

% 
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