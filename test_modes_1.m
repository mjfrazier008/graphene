function modes = test_modes_1(By, y, x)
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
% for n = 1:d
%     if real(ein(n)) > tol
%         f1 = f1 + 1;
%         decayleft(f1) = n;
%     elseif abs(real(ein(n))) < tol && imag(ein(n)) > tol
%         f2 = f2 + 1;
%         waveleft_in(f2) = n;
%     elseif abs(real(ein(n))) < tol && imag(ein(n)) < tol
%         f3 = f3 + 1;
%         waveleft_out(f3) = n;
%     end
%     if real(eout(n)) < -tol
%         f4 = f4 + 1;
%         decayright(f4) = n;
%     elseif abs(real(eout(n))) < tol && imag(eout(n)) < tol
%         f5 = f5 + 1;
%         waveright_in(f5) = n;
%     elseif abs(real(eout(n))) < tol && imag(eout(n)) > tol
%         f6 = f6 +1;
%         waveright_out(f6) = n;
%     end
% end
decayleft = decayleft(1:f1);
waveleft = waveleft(1:f2);
decayright = decayright(1:f3);
waveright = waveright(1:f4);
h = round(size(y, 2)/2);
yleft = y(1:h-1);
yright = flip(y(h+1:end));
vin = vin(:, [decayleft, waveleft]);
vout = vout(:, [decayright, waveright]);
vleft = runge(vin, By, yleft, dx);
vright = runge(vout, By, yright, -dx);
modes = 2*ones([1, 4]);
Qleft = orth(vleft);
Qright = orth(vright);
modes(1) = 1-eigs((Qleft*Qleft')*(Qright*Qright'), 1, 1);

% Qleft_in = orth(vleft(:, [waveleft_in, waveleft_out]));
% Qleft_out = orth(vleft(:, [decayleft, waveleft_out]));
% Qright_in = orth(vright(:, [waveright_in, waveright_out]));
% Qright_out = orth(vright(:, [waveright_out, decayright]));
% Qleft_decay = orth(vleft(:, decayleft));
% Qright_decay = orth(vright(:, decayright));
% modes = 2*ones([1, 4]);
% modes(1) = 1 - eigs((Qright_decay*Qright_decay')*(Qleft_decay*Qleft_decay'), 1, 1);
% modes(2) = 1 - eigs((Qright_in*Qright_in')*(Qleft_out*Qleft_out'), 1, 1);
% modes(3) = 1 - eigs((Qleft_in*Qleft_in')*(Qright_out*Qright_out'), 1, 1);
