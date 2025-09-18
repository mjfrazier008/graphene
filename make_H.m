function [H, v, e] = make_H(kx, ky, om, eps, nlayers)
    U = [[0, 1];[0, 0]];
    s0 = eye(2);
    s1 = [[0, 1];[1, 0]];
    s2 = [[0, -1j];[1j, 0]];
    t = zeros([1, nlayers]);
    t(2) = 1;
    T1 = toeplitz(zeros([1, nlayers]), t);
    T2 = toeplitz(t, zeros([1, nlayers]));
    om = om/2*diag(linspace(1, -1, nlayers));
    H = kron(T1, eps*U') + kron(T2, eps*U) + kron(om, s0) +...
        kron(diag(1:nlayers), kx*s1+ky*s2);
    [v, e] = eig(H);