%% Multiple Dimensional Scaling (MDS)
% |function[Y, e] = mds(D, p)|
%
% * Author:   Shangkun.Shen
% * Method:   Eigenvalue Decomposition
%
%% *Usage*
% * |[Y, e] = mds(D)|
% * |[Y, e] = mds(D, p)|
%
% Example
%
% # |[Y, e] = mds(pdist([xx, yy, zz]), 2)|
% # |[Y, e] = mds(pdist([xx, yy, zz], 'cityblock'), 2)|
% # |[Y, e] = mds(pdist( [xx, yy, zz], 
% @(Xi, Xj)(sqrt(bsxfun(@minus, Xi, Xj) .^ 2 * [.2 .5 .3]'))), 2)|
%
%% *Input Arguments*
% * |D|: n-by-n distance matrix recommand, or a more general dissimilarity
% matrix;
% * |p|: reduced data demension.
%% *Output Results*
% * |Y|: n-by-p matrix, rows of |Y| are the coordinates of n points in
% p-dimensional space for some p < n, where n equals to the size of |D|;
% * |e|: the eigenvalues of |Y * Y'|.
%
%% *Source Code*
function [ Y, e ] = mds( D, p )
    if (nargin < 1)
        warning 'Lack of input arguments. Terminates.'; return;
    elseif (nargin > 2)
    	warning 'Too many input arguments. Terminates.'; return;
    end
    
    [n, m] = size(D);
    
    if (nargin == 2 && isscalar(p))
    	p = min(p, n);
    elseif (nargin == 1)
    	p = n;
    end
    
    if (n == 1)
    	n = ceil(sqrt(2 * m));
    	if (n * (n - 1) / 2 == m) && all(D >= 0)
    		D = squareform(D);
    	else
    		warning 'Bad input argument D. Terminates.'; return;
    	end
    elseif (n ~= m)
    	warning 'Bad input argument D. Terminates.'; return;
    end

    % P = eye(n) - repmat(1/n,n,n);
    % B = P * (-.5 * D .* D) * P;
    % A more efficient way of doing the same thing.
    D = D .* D; % square elements of D
    B = bsxfun( @plus, ...
            bsxfun( @minus, ...
                bsxfun( @minus, D, sum(D, 1) / n ), ...
                sum(D, 2) / n ), ...
            sum(D(:)) / (n ^ 2) ) * (-0.5);
    [V, E] = eigs( (B + B') / 2, p ); e = diag(E);
    keep = find(e > max(abs(e)) * eps(class(e)) ^ (3 / 4));

    if isempty(keep)
    	Y = zeros(n, 1);
    else
    	Y = bsxfun(@times, V(:, keep), sqrt(e(keep))');
    end
    
    [~, maxind] = max(abs(Y), [], 1);
    d = size(Y, 2);
    colsign = sign(Y(maxind + (0 : n : (d - 1) * n)));
    Y = bsxfun(@times, Y, colsign);
end