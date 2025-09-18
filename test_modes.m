function modes = test_modes(By, y, x)
tol = 10^(-8);
[vin, ein] = eig(By(y(1)));
[vout, eout] = eig(By(y(end)));
ein = diag(ein);
eout = diag(eout);
dx = x(2)-x(1);
d = size(vin, 2);
decayleft = zeros([1, d]);
waveleft = zeros([1, d]);
decayright = zeros([1, d]);
waveright = zeros([1, d]);
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
vleft = runge(vin, By, yleft, dx);
vright = runge(vout, By, yright, -dx);
Qdecayr = orth(vright(:, decayright));
Qwaver = orth(vright(:, waveright));
Qdecayl = orth(vleft(:, decayleft));
Qwavel = orth(vleft(:, waveleft));
modes = 2*ones([1, 4]);
modes(1) = abs(eigs((Qdecayr*Qdecayr')*(Qdecayl*Qdecayl'), 1, 1) -1);
modes(2) = abs(eigs((Qwaver*Qwaver')*(Qdecayl*Qdecayl'), 1, 1) -1);
modes(3) = abs(eigs((Qdecayr*Qdecayr')*(Qwavel*Qwavel'), 1, 1) -1);
modes(4) = abs(eigs((Qwaver*Qwaver')*(Qwavel*Qwavel'), 1, 1) -1);
% v1 = Qdecay*Qdecay'*vleft;
% v2 = Qwave*Qwave'*vleft;
% if f1 > 0
%     for i = 1:f1
%         if f3 > 0
%             modes(1) = norm(vleft(:, decayleft(i)) - v1(:, decayleft(i)));
%         end
%         if f4 > 0
%             modes(2) = norm(vleft(:, decayleft(i)) - v2(:, decayleft(i)));
%         end
%     end
% end
% if f2 > 0
%     for i = 1:f2
%         if f3 > 0
%             modes(3) = norm(vleft(:, waveleft(i)) - v1(:, waveleft(i)));
%         end
%         if f4 > 0
%             modes(4) = norm(vleft(:, waveleft(i)) - v2(:, waveleft(i)));         
%         end
%     end
% end
