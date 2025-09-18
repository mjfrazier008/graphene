function [Ae, Al, Ar, Alr] = grid_adapt(kbound, Ebound, Nk_base, NE_base, tolbase, L,...
    B, f, x)
% Assuming odd number of grid points for simplicity:
if mod(Nk_base, 2) == 1
    Nk_base = Nk_base - 1;
end
if mod(NE_base, 2) == 1 
    NE_base = NE_base - 1;
end
% Calculate base case:
[Ae, Al, Ar, Alr] = residue_map(kbound, Ebound, Nk_base+1, NE_base+1, B, f, x);
% Overall acceptable error:
tol = tolbase*2^(-L);
% Total number of grid points for tracking:
Ntot = (Nk_base*2^L+1)*(NE_base*2^L+1);
for j = 1:L
    Nk = Nk_base*2^j+1;
    NE = NE_base*2^j+1;
    kx = linspace(kbound(1), kbound(2), Nk);
    e = linspace(Ebound(1), Ebound(2), NE);
    Aenew = zeros(NE, Nk);
    Alnew = zeros(NE, Nk);
    Arnew = zeros(NE, Nk);
    Alrnew = zeros(NE, Nk);
    % Initialize new grid with twice as many points:
    Aenew(1:2:NE, 1:2:Nk) = Ae;
    Alnew(1:2:NE, 1:2:Nk) = Al;
    Arnew(1:2:NE, 1:2:Nk) = Ar;
    Alrnew(1:2:NE, 1:2:Nk) = Alr;
    Ae = Aenew;
    Al = Alnew;
    Ar = Arnew;
    Alr = Alrnew;
    % Upper tolerance for level:
    tol_l = tolbase*2^(1-j);
    for i = 1:2:NE_base*2^(j)
        for k = 1:2:Nk_base*2^(j)
            % 2 is maximum difference in unit vectors so a sum of 2*16
            % possile residues means no valid states
            S = Alr(i, k) + Alr(i+2, k) + Alr(i, k+2) + Alr(i+2, k+2) +...
                Al(i, k) + Al(i+2, k) + Al(i, k+2) + Al(i+2, k+2)+ ...
                Ar(i, k) + Ar(i+2, k) + Ar(i, k+2) + Ar(i+2, k+2) +...
                Ae(i, k) + Ae(i+2, k) + Ae(i, k+2) + Ae(i+2, k+2);
            % For each old grid point, check if (old) neighbors are within
            % overall tolerance or above level tolerance and interpolate if
            % so.
            if S == 32 ||...
                    (Alr(i, k) > tol_l && Alr(i+2, k) > tol_l &&...
                    Alr(i, k+2) > tol_l && Alr(i+2, k+2) > tol_l &&...
                    Al(i, k) > tol_l && Al(i+2, k) > tol_l && Al(i, k+2) > tol_l...
                    && Al(i+2, k+2) > tol_l && Ar(i, k) > tol_l &&...
                    Ar(i+2, k) > tol_l && Ar(i, k+2) > tol_l && Ar(i+2, k+2) > tol_l...
                    && Ae(i, k) > tol_l && Ae(i+2, k) > tol_l && Ae(i, k+2) > tol_l...
                    && Ae(i+2, k+2) > tol_l) ||...
                    (Alr(i,k) < tol && Alr(i+2, k) < tol && Alr(i, k+2) < tol && Alr(i+2, k+2) < tol) ||...
                    (Al(i,k) < tol && Al(i+2, k) < tol && Al(i, k+2) < tol && Al(i+2, k+2) < tol) ||...
                    (Ar(i,k) < tol && Ar(i+2, k) < tol && Ar(i, k+2) < tol && Ar(i+2, k+2) < tol) ||...
                    (Ae(i,k) < tol && Ae(i+2, k) < tol && Ae(i, k+2) < tol && Ae(i+2, k+2) < tol)
                Alr(i+1, k) = mean([Alr(i,k), Alr(i+2, k)]);
                Alr(i, k+1) = mean([Alr(i, k), Alr(i, k+2)]);
                Alr(i+1, k+1) = mean([Alr(i,k), Alr(i+2, k+2), Alr(i+2, k),...
                    Alr(i, k+2)]);
                Al(i+1, k) = mean([Al(i,k), Al(i+2, k)]);
                Al(i, k+1) = mean([Al(i, k), Al(i, k+2)]);
                Al(i+1, k+1) = mean([Al(i,k), Al(i+2, k+2), Al(i+2, k),...
                    Al(i, k+2)]);
                Ar(i+1, k) = mean([Ar(i,k), Ar(i+2, k)]);
                Ar(i, k+1) = mean([Ar(i, k), Ar(i, k+2)]);
                Ar(i+1, k+1) = mean([Ar(i,k), Ar(i+2, k+2), Ar(i+2, k),...
                    Ar(i, k+2)]);
                Ae(i+1, k) = mean([Ae(i,k), Ae(i+2, k)]);
                Ae(i, k+1) = mean([Ae(i, k), Ae(i, k+2)]);
                Ae(i+1, k+1) = mean([Ae(i,k), Ae(i+2, k+2), Ae(i+2, k),...
                    Ae(i, k+2)]);
                if i == NE_base*2^j
                    Alr(i+2, k+1) = mean([Alr(i+2, k), Alr(i+2, k+2)]);
                    Al(i+2, k+1) = mean([Al(i+2, k), Al(i+2, k+2)]);
                    Ar(i+2, k+1) = mean([Ar(i+2, k), Ar(i+2, k+2)]);
                    Ae(i+2, k+1) = mean([Ae(i+2, k), Ae(i+2, k+2)]);
                end
                if k == Nk_base*2^j
                    Alr(i+2, k+2) = mean([Alr(i, k+2), Alr(i+2, k+2)]);
                    Al(i+2, k+2) = mean([Al(i, k+2), Al(i+2, k+2)]);
                    Ar(i+2, k+2) = mean([Ar(i, k+2), Ar(i+2, k+2)]);
                    Ae(i+2, k+2) = mean([Ae(i, k+2), Ae(i+2, k+2)]);
                end
            end
        end
    end
    Level = j-1
    percent_done = nnz(Ae)/Ntot
    % Calculate residues that weren't interpolated:
    % CHANGE TO FOR LOOP HERE IF NO PARALLEL CAPABILITY:
    parfor i = 1:NE
        for k = 1:Nk
            if Alr(i, k) == 0
                En = e(i);
                kn = kx(k);
                H = B(En, kn);
                modes = test_modes_1(H, f, x);
                Ae(i, k) = modes(1);
                Ar(i, k) = modes(2);
                Al(i, k) = modes(3);
                Alr(i,k) = modes(4);
            end
        end
    end
end

    