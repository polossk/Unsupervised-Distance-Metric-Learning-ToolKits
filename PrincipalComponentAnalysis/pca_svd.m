%% Principal Component Analysis (PCA), Using Singular Value Decomposition.
% |function[D, W, mu] = pca_svd(X, eta, flag)|
%
% |function[D, W, mu] = pca_svd(X, num, flag)|
%
% * Author:   Shangkun.Shen
% * Method:   Singular Value Decomposition
%
%% *Usage*
%
% * |[D, W, mu] = pca_svd(X, eta);|
% * |[D, W, mu] = pca_svd(X, eta, '-rebuilt');|
% * |[D, W, mu] = pca_svd(X, eta, '-all');|
% * |[D, W, mu] = pca_svd(X, eta, '-single');|
% * |[D, W, mu] = pca_svd(X, num, '-exact');|
% * |[D, W, mu] = pca_svd(X, num, '-exact -all');|
% * |[D, W, mu] = pca_svd(X, num, '-exact -single');|
%
% Example
%
% # |[D, W, mu] = pca_svd(X, 20, '-exact -all');|
% # |[D, W, mu] = pca_svd(X, 0.95, '-single');|
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
function [D, W, mu] = pca_svd( X, num, varargin )
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

    [w, ~] = size(X);
    X = bsxfun(@minus, X, mu);
    [u, s_, ~] = svd(X, 'econ');
    s = diag(s_);
    l = length(s);

    if strcmp(setting_num, flags{4})
        D = zeros(num, 1); W = zeros(w, num);
        if l < num
            D(1 : l) = s; W(:, 1 : l) = u;
        else
            D = s(1 : num); W = u(:, 1 : num);
        end
    else
        n = min(find(cumsum(s ./ sum(s))) > num);
        W = u(:, 1 : n); D = D(:, 1 : n);
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