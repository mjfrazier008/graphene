function [modes, V, Eleft, Eright] = test_modes_vec(By, y, x)
tol = 10^(-8);
[vin, ein] = eig(By(y(1)));
[vout, eout] = eig(By(y(end)));
ein = diag(ein);
eout = diag(eout);
dx = x(2)-x(1);
d = size(vin, 2);
decayleft = zeros([1, d]);
waveleft = zeros([1, d]);
% waveleft_out = zeros([1, d]);
decayright = zeros([1, d]);
waveright = zeros([1, d]);
%waveright_out = zeros([1, d]);
f1 = 0;
f2 = 0;
f3 = 0;
f4 = 0;
for n = 1:d
    if real(ein(n)) > tol
        f1 = f1 + 1;
        decayleft(f1) = n;
    elseif abs(real(ein(n))) < tol && abs(imag(ein(n))) > tol
        f2 = f2 + 1;
        waveleft(f2) = n;
    end
    if real(eout(n)) < -tol
        f3 = f3 + 1;
        decayright(f3) = n;
    elseif abs(real(eout(n))) < tol && abs(imag(eout(n))) > tol
        f4 = f4 + 1;
        waveright(f4) = n;
    end
end

decayleft = decayleft(1:f1);
waveleft = waveleft(1:f2);
decayright = decayright(1:f3);
waveright = waveright(1:f4);
h = round(size(y, 2)/2);
yleft = y(1:h-1);
yright = flip(y(h+1:end));
vin = vin(:, [decayleft, waveleft]);
vout = vout(:, [decayright, waveright]);
[vleft, Vleft] = runge_vec(vin, By, yleft, dx);
[vright, Vright] = runge_vec(vout, By, yright, -dx);
modes = 2*ones([1, 4]);
[Qleft, Rleft] = qr(vleft, "econ");
[Qright, Rright] = qr(vright, "econ");
[p, e] = eigs((Qleft*Qleft')*(Qright*Qright'), 1, 1);
modes(1) = 1-e;
vt_left = Rleft\(Qleft'*p);
vt_right = Rright\(Qright'*p);
vec_left = zeros(size(Vleft, 1), size(Vleft, 3));
vec_right = zeros(size(Vright, 1), size(Vright, 3));
for j = 1:size(Vleft, 3)
    vec_left(:, size(Vleft, 3)+1-j) = Vleft(:, :, size(Vleft, 3)+1-j)*vt_left;
    %vec_left(:, size(Vleft, 3)+1-j) = vec_left(:, size(Vleft, 3)+1-j)/norm(vec_left(:, size(Vleft, 3)+1-j));
    vec_right(:, j) = Vright(:, :, size(Vright, 3) + 1-j)*vt_right;
    %vec_right(:, j) = vec_right(:, j)/norm(vec_right(:, j));
end


Eleft = ein([decayleft, waveleft]);
Eright = eout([decayright, waveright]);
V = [vec_left, p, vec_right];
