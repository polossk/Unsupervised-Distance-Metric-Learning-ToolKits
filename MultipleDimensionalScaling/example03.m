len = 1000;
mark_size = 32;
a = 0.3;
theta = 4 * pi * rand(len, 1);
h = 4 * rand(len, 1);
r = a * theta;
x = r .* cos(theta);
y = r .* sin(theta);
z = h;
colormap('winter');
scatter3(x, y, z, mark_size, r);


alpha1 = pi / 4;
alpha2 = pi / 12;
alpha3 = -pi / 2;
A1 = [
    cos(alpha2), sin(alpha1), 0, 0;
    sin(alpha1), -cos(alpha1), 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1
];
A2 = [
    cos(alpha2), 0, sin(alpha2), 0;
    0, 1, 0, 0;
    sin(alpha2), 0, -cos(alpha2), 0;
    0, 0, 0, 1
];
A3 = [
    1, 0, 0, 0;
    0, cos(alpha3), sin(alpha3), 0;
    0, sin(alpha3), -cos(alpha3), 0;
    0, 0, 0, 1
];

T = A1 * A2 * A3;

PP = [x, y, z, ones(len, 1)] * T;
xx = PP(:, 1);
yy = PP(:, 2);
zz = PP(:, 3);

figure(2);
colormap('winter');
scatter3(xx, yy, zz, mark_size, r);

geoL = @( a, phi ) ( ...
    a .* ( ...
        phi .* 0.5 .* sqrt( 1 + phi .^ 2 ) + ...
        0.5 .* log(phi + sqrt( 1 + phi .^ 2 ) ) ...
        ) ...
    );

pX = geoL(a, theta);
pY = h;
D = squareform(pdist([pX, pY]));
[YY, e] = cmdscale(D);
figure(3);
colormap('winter');
scatter(YY(:, 1), YY(:, 2), mark_size, r);

[YYY, e] = mds(squareform(pdist([xx, yy, zz])));
figure(4);
colormap('winter');
scatter3(YYY(:, 1), YYY(:, 2), YYY(:, 3), mark_size, r);