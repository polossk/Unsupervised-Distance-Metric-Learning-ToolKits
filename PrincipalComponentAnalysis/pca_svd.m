%% Principal Component Analysis (PCA), Using Singular Value Decomposition.
% |function[D, W, mu] = pca_svd(X, n)|
%
%% *Input*
% * |X|: N by P data matrix. Rows of X correspond to observations
% and columns to variables.
% * |n|: number of principal component.
%
%% *Output*
% * |D|: the principal component variances
% * |W|: the principal component coefficients
% * |mu|: average of X
%
%% *Source*
function [D, W, mu] = pca_svd(X, n)
    mu = mean(X(:));
    [w, ~] = size(X);
    bsxfun(@minus, X, mu);
    [u, s, ~] = svd(X, 0);
    [~, l] = size(s);
    D = zeros(n, n);
    W = zeros(w, n);
    if l < n
        D(1 : l, 1 : l) = s;
        W(  :  , 1 : l) = u;
    else
        D = s(1 : n, 1 : n);
        W = u(  :  , 1 : n);
    end
end
%% *Note*
% The principal component score didn't show up in this code. Let
% variable |score| be the principal component score. Let variable
% |coeff| be the principal component coefficients. The centered data
% can be reconstructed by |SCORE * COEFF'|.
%
% Like this code below:
%
% |bsxfun(@minus, X, mu); score = X / W';|