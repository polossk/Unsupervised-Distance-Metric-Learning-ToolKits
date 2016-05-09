%% Principal Component Analysis (PCA), Using Eigenvalue Decomposition.
% |function[D, W, mu] = pca_eig(X, eta, flag)|
%
% |function[D, W, mu] = pca_eig(X, num, flag)|
%
% * Author:   Shangkun.Shen
% * Method:   Eigenvalue Decomposition
%
%% *Usage*
%
% * |[D, W, mu] = pca_eig(X, eta);|
% * |[D, W, mu] = pca_eig(X, eta, '-rebuilt');|
% * |[D, W, mu] = pca_eig(X, eta, '-all');|
% * |[D, W, mu] = pca_eig(X, eta, '-single');|
% * |[D, W, mu] = pca_eig(X, num, '-exact');|
% * |[D, W, mu] = pca_eig(X, num, '-exact -all');|
% * |[D, W, mu] = pca_eig(X, num, '-exact -single');|
%
% Example
%
% # |[D, W, mu] = pca_eig(X, 20, '-exact -all');|
% # |[D, W, mu] = pca_eig(X, 0.95, '-single');|
%
%% *Input Arguments*
% * |X|: N by P data matrix. Rows of X correspond to observations
% and columns to variables;
% * |eta|: coefficient of rebuilt data;
% * |num|: number of principal component;
% * |flag|: pca option: '-single' by default, setting |mu| to the average
% of each columns, also could be '-all', setting |mu| to the average of the
% whole data; '-rebuilt' by default, take the second argument as |eta|, or
% as |num| by using '-exact' instead.
%
%% *Output Results*
% * |D|: the principal component variances
% * |W|: the principal component coefficients
% * |mu|: average of X, determined by |flag|
%
%% *Source Code*
function [D, W, mu] = pca_eig( X, num, varargin )
    if (nargin < 2)
        warning 'Lack of input arguments. Terminates.'; return;
    elseif (nargin > 3)
        warning 'Too many input arguments. Terminates.'; return;
    elseif (nargin == 3 && ischar(varargin{1}))
        flag = varargin{1};
    elseif (nargin == 2)
        flag = '-single -exact';
    else
        warning 'Wrong input arguments. Terminates.'; return;
    end
    
    setting = strsplit(flag);

    if (length(setting) > 2)
        warning 'Too many input arguments. Terminates.'; return;
    end

    flags = {'-single', '-all', '-rebuilt', '-exact'};

    if any(strcmp(setting, flags{1})) && any(strcmp(setting, flags{2}))
        warning 'Wrong flag setting. Terminates.'; return;
    elseif ~any(strcmp(setting, flags{1})) && any(strcmp(setting, flags{2}))
        setting_mu = flags{2};
    else
        setting_mu = flags{1};
    end

    if any(strcmp(setting, flags{3})) && any(strcmp(setting, flags{4}))
        warning 'Wrong flag setting. Terminates.'; return;
    elseif ~any(strcmp(setting, flags{3})) && any(strcmp(setting, flags{4}))
        setting_num = flags{4};
    else
        setting_num = flags{3};
    end

    if strcmp(setting_mu, flags{2})
        mu = mean(X(:));
    else
        mu = mean(X, 1);
    end

    X = bsxfun(@minus, X, mu);

    if strcmp(setting_num, flags{4})
        [W, D_] = eigs(X * X', num);
        D = diag(D_);
    else
        [W, D_] = eig(X * X');
        [D, idx] = sort(diag(D_), 'descend');
        n = min(find(cumsum(D ./ sum(D))) > num);
        W = W(:, idx); W = W(:, 1 : n); D = D(:, 1 : n);
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