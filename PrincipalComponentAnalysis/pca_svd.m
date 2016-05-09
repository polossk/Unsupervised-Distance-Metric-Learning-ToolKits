%% Principal Component Analysis (PCA), Using Singular Value Decomposition.
% |function[D, W, mu] = pca_svd(X, n)|
%
% * Author:   Shangkun.Shen
% * Method:   Singular Value Decomposition
%
%% *Usage*
%
% * |[D, W, mu] = pca_svd(X, n);|
%
% * |[D, W, mu] = pca_svd(X, n, '-all');|
%
% * |[D, W, mu] = pca_svd(X, n, '-single');|
%
%% *Input Arguments*
% * |X|: N by P data matrix. Rows of X correspond to observations
% and columns to variables.
% * |n|: number of principal component.
% * |flag|: pca option, '-single' by default, setting |mu| to the average
% of each columns, also could be '-all', setting |mu| to the average of the
% whole data.
%
%% *Output Results*
% * |D|: the principal component variances
% * |W|: the principal component coefficients
% * |mu|: average of X, determined by |flag|
%
%% *Source Code*
function [D, W, mu] = pca_svd( varargin )
    if (nargin < 2)
        warning 'Lack of input arguments. Terminates.'; return;
    elseif (nargin > 3)
        warning 'Too many input arguments. Terminates.'; return;
    elseif (isnumeric(varargin{1}) && isnumeric(varargin{2}))
        X = varargin{1}; n = varargin{2};
        if (nargin == 3 && ischar(varargin{3}))
            flag = varargin{3};
        else
            warning 'Wrong input arguments. Terminates.'; return;
        end
    else 
        warning 'Wrong input arguments. Terminates.'; return;
    end

    if strcmp('-all', flag)
        mu = mean(X(:));
    elseif strcmp('-single', flag)
        mu = mean(X, 1);
    else
        warning 'Wrong flag setting. Terminates.'; return;
    end

    [w, ~] = size(X);
    X = bsxfun(@minus, X, mu);
    [u, s, ~] = svd(X, 0);
    [~, l] = size(s);
    D = zeros(n, n); W = zeros(w, n);
    if l < n
        D(1 : l, 1 : l) = s; W(  :  , 1 : l) = u;
    else
        D = s(1 : n, 1 : n); W = u(  :  , 1 : n);
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
% |X = bsxfun(@minus, X, mu); score = X / W';|