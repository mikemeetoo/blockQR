function [R, V, Tau] = simpleQR(A)
    [m, n] = size(A);
    
    V = zeros(m,n);
    Tau = zeros(1,n);
    
    for k = 1:n,
        x = A(k:end,k);
        v = x + sign(x(1))*norm(x)*eye(m-k+1,1);
        tau = v'*v/2;
        V(k:end,n-k+1) = v;
        Tau(1,n-k+1) = tau;
        A(k:end,k:end) = A(k:end,k:end) - (1/tau)*v*(v'*A(k:end,k:end));
    end
    R = A;
end