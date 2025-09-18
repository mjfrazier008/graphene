function B = graphene_om(E, kx)

% Change parameters as desired:
nlayers = 5;
Om = 1.2;
eps = 1;

% Do not change:
s1 = [[0, 1];[1, 0]];
s2 = [[0, -1j];[1j, 0]];
t = zeros([1, nlayers]);
t(2) = 1;
T1 = toeplitz(zeros([1, nlayers]), t);
T2 = toeplitz(t, zeros([1, nlayers]));
A = [[0, 1];[0, 0]];
om = Om/2*diag(linspace(1, -1, nlayers));
% if mod(nlayers, 2) == 1
%         om = Om*diag(linspace((nlayers-1)/2, -(nlayers-1)/2, nlayers));
%     else
%         om = Om*diag(linspace((nalyers-1)/2))
% end
a = kron(eye(nlayers), s2);
B = @(m) -1j*a*(kron(T1, eps*A') +...
    kron(T2, eps*A) + kron(eye(nlayers), kx*s1) +...
    kron(0.5*((1-m)*om - (1+m)*om), eye(2)) - E*eye(2*nlayers));
end