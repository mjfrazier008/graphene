function [vout, V] = runge_vec(vin, By, y, dx)
n = size(y, 2);
[k, m] = size(vin);
vout = vin;
V = zeros(k, m, n);
for j = 1:n-1
    for k = 1:m
        V(:, k, j) = vout(:, k);
    end
    B0 = By(y(j));
    B1 = By(y(j) + dx/2);
    B2 = By(y(j+1));
    k1 = B0*vout;
    k2 = B1*(vout + (dx/2)*k1);
    k3 = B1*(vout + (dx/2)*k2);
    k4 = B2*(vout + dx*k3);
    vout = vout + (dx/6)*(k1 + 2*k2 + 2*k3 + k4);
end
for k = 1:m
        V(:, k, n) = vout(:, k);
end
%for j = 1:m
%    vout(:, j) = vout(:, j)/norm(vout(:, j));
%end