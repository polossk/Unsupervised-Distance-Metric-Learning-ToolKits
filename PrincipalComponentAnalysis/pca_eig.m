%% Principal Component Analysis (PCA), Using Eigenvalue Decomposition.
% |function[D, W, mu] = pca_eig(X, n)|
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
function [D, W, mu] = pca_eig(X, n)
    mu = mean(X(:));
    bsxfun(@minus, X, mu);
    [W, D] = eigs(X * X', n);
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