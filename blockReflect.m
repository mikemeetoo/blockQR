function [Q R] = blockReflect(A)
    m = size(A,1);
    half = m/2; 
 
    if m == 1,
        R = A;
        Q = 1;
    else 
        E1 = A(1:half,1:half);
        E2 = A(half+1:end, 1:half);

        X = (E1'\E2')';

        [V D] = eig(eye(half)+X'*X);
        F = V*sqrt(D)*V'*E1;

        D = [F-E1;-E2];
        H = eye(half*2) - 2*D*pinv(D);
        A(:,1:half) = [F; zeros(half,half)];
        A(:,half+1:end) = A(:,half+1:end) - 2*D*pinv(D)*A(:,half+1:end);

        [Q1 A(1:half,:)] = blockReflect(A(1:half,:));
        [Q2 A(half+1:end,half+1:end)] = blockReflect(A(half+1:end,half+1:end));
        
        R = A;
        Q = H'*[Q1 zeros(half,half);zeros(half,half) Q2];
    end
    
end
 