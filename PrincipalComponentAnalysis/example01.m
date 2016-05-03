%% Example 01 Using PCA to reduce the demension of face pic
% * used 20 eigvalues of each pic
% * maxinum error rate less than 3%
%
%% *READ DATA*
if (exist('image_set.mat', 'file') ~= 2)
    files = dir('*.happy.gif');
    n = length(files);
    info = imfinfo(files(1).name, 'gif');
    im_map = info.ColorTable;
    im_w = int32(info.Width);
    im_h = int32(info.Height);
    im_set = zeros(n, im_h, im_w, 'uint8');
    for id = 1 : n
        im_set(id, :, :) = imread(files(id).name, 'gif');
    end
    save('image_set.mat', 'im_set', 'im_map');
else
    load('image_set.mat');
end
[n, im_h, im_w] = size(im_set);
im_map_length = size(im_map', 2);
x = zeros(im_h, im_w);
y = zeros(im_h, im_w);

%% *Create Picture with First 20 Eigvalues and Eigvectors*
%  Using SVD
im_svd_set = zeros(n, im_h, im_w, 'uint8');
fprintf('ERROR Summary: Using svd method\n');
name_prefix = 'subject';
name_suffix = '.happy.svd.gif';
for id = 1 : n
    im = reshape(im_set(id, :, :), im_h, im_w);
    x = double(im);
    [D, W, mu] = pca_svd(x, 20);
    y = W * W' * x;
    gtcmap = find(y >= im_map_length);
    ltcmap = find(y < 0);
    if (~isempty(gtcmap) || ~isempty(ltcmap))
        y(gtcmap) = im_map_length - 1;
        y(ltcmap) = 0;
    end
    err = abs(im(:) - y(:)); y = uint8(y);
    err_max = max(err);
    err_avg = mean(err);
    err_ratio = err_avg / double(im_map_length);
    fprintf('sample %02d: maxinum error: %3.0f, average error: %4.2f, error rate: %4.2f%%.\n', ...
        id, err_max, err_avg, err_ratio * 100);
    im_svd_set(id, :, :) = y;
    imwrite(y, im_map, sprintf('%s%02d%s', name_prefix, id, name_suffix));
end
save('image_svd_set.mat', 'im_svd_set');

%% *Create Picture with First 20 Eigvalues and Eigvectors*
%  Using EIGS
fprintf('ERROR Summary: Using eigs method\n');
im_eig_set = zeros(n, im_h, im_w, 'uint8');
name_prefix = 'subject';
name_suffix = '.happy.eig.gif';
for id = 1 : n
    im = reshape(im_set(id, :, :), im_h, im_w);
    x = double(im);
    [D, W, mu] = pca_eig(x, 20);
    y = W * W' * x;
    gtcmap = find(y >= im_map_length);
    ltcmap = find(y < 0);
    if (~isempty(gtcmap) || ~isempty(ltcmap))
        y(gtcmap) = im_map_length - 1;
        y(ltcmap) = 0;
    end
    err = abs(im(:) - y(:));
    err_max = max(err);
    err_avg = mean(err);
    err_ratio = err_avg / double(im_map_length);
    fprintf('sample %02d: maxinum error: %3.0f, average error: %4.2f, error rate: %4.2f%%.\n', ...
        id, err_max, err_avg, err_ratio * 100);
    im_eig_set(id, :, :) = y;
    imwrite(uint8(y), im_map, sprintf('%s%02d%s', name_prefix, id, name_suffix));
end
save('image_eig_set.mat', 'im_eig_set');

%% *Plot All Picture*
for id = 1 : n / 3
    figure(id);
    for ii = 1 : 3
        subplot(3, 3, ii * 3 - 2);
        im = reshape(im_set(id * 3 + ii - 3, :, :), im_h, im_w);
        imshow(im, im_map);
        title('Before PCA');
        subplot(3, 3, ii * 3 - 1);
        im = reshape(im_svd_set(id * 3 + ii - 3, :, :), im_h, im_w);
        imshow(im, im_map);
        title('After PCA, using svd');
        subplot(3, 3, ii * 3);
        im = reshape(im_eig_set(id * 3 + ii - 3, :, :), im_h, im_w);
        imshow(im, im_map);
        title('After PCA, using eigs');
    end
    saveas(gcf, strcat('example01-', num2str(id)), 'png');
end