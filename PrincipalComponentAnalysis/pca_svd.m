function [D, W, mu] = pca_svd(X, n)
    mu = mean(X(:));
    [w, h] = size(X);
    bsxfun(@minus, X, mu);
    [u, s, v] = svd(X, 0);
    [l, l] = size(s);
    sn = zeros(n, n);
    un = zeros(w, n);
    if l < s
        sn(1 : l, 1 : l) = s;
        un(  :  , 1 : l) = u;
    else
        sn = s(1 : n, 1 : n);
        un = u(  :  , 1 : n);
    end
    D = sn;
    W = un;
end