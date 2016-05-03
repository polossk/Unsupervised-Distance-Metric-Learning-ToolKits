function [D, W, mu] = pca_eig(X, n)
    mu = mean(X(:));
    bsxfun(@minus, X, mu);
    [W, D] = eigs(X * X', n);
end