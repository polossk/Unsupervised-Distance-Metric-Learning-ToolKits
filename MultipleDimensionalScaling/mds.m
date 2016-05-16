function [ Y, e ] = mds( D )
    [n, ~] = size(D);
    % P = eye(n) - repmat(1/n,n,n);
    % B = P * (-.5 * D .* D) * P;
    % A more efficient way of doing the same thing.
    D = D .* D; % square elements of D
    B = bsxfun( @plus, ...
            bsxfun( @minus, ...
                bsxfun( @minus, D, sum(D, 1) / n ), ...
                sum(D, 2) / n ), ...
            sum(D(:)) / (n ^ 2) ) * (-0.5);
    [V, E] = eig( (B + B') / 2 );
    [e, idx] = sort(diag(E), 'descend');
    keep = find(e > max(abs(e)) * eps(class(e))^(3/4));
    Y = bsxfun(@times, V(:, idx(keep)), sqrt(e(keep))');
    [~, maxind] = max(abs(Y), [], 1);
    d = size(Y, 2);
    colsign = sign(Y(maxind + (0 : n : (d - 1) * n)));
    Y = bsxfun(@times, Y, colsign);
end