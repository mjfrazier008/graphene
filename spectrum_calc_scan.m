
% Choose energy and wavenumber bounds to calculate spectrum over:
kbound = [-3, 3];
Ebound = [-0.1, 0.1];

% Discretization parameters. Nk, Ne; initial size of (k, E) grid to
% calculate over. L; number of adaptive grid steps. Number of grid points
% in each direction is doubled each step so that the final grid size is
% Nk*2^(L-1) by Ne*2^(L-1). N is number of discretization points for the
% ODE solver to use at each point.
Nk = 750;
Ne = 300;
N = 1000;
L = 3;

% Define a switch funtion from -1 to 1:
x = linspace(-1, 1, N+1);
m = -1+2./(1+exp(-10*x));
%m = spline([x(1), x(N/2+1), x(N+1)], [0, x(1), x(N/2+1), x(N+1), 0], x);

% Choose to vary epsilon or Omega by choosing @graphene_eps or @graphene_om
% respectively:
B = @graphene_eps;
tic;

% Choose initial error e and calculate spectrum. Final error guaranteed is e*2^(1-L).
e = 0.01;
[Ae, AL, AR, ALR] = grid_adapt(kbound, Ebound, Nk, Ne, .01, L, B, m, x);
% Functions "residue_map.m" and "uniform_grid.m" may also be used but do
% not perform as well.
t = toc/60

% Plot spectrum:
figure, imshow(Ae, [0, .1]);
k = linspace(kbound(1), kbound(2), Nk*2^L+1);
E = linspace(Ebound(1), Ebound(2), Ne*2^L+1);
plot_grid(Ae, AL, AR, ALR, k, E, 10^(-4));