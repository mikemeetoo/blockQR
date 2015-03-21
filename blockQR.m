function [Q, R] = blockQR(A, block_size)
    [m,n] = size(A);
    Q = eye(m,n);
    
    for ptr = 1:block_size:n,
        if ptr+block_size-1 > n,

            [A(ptr:end,ptr:end), V, Tau] = simpleQR(A(ptr:end,ptr:end));
            S = CWYTransform(V, Tau);
            Q(ptr:end,:) = Q(ptr:end,:) - V*(S*(V'*Q(ptr:end,:))); 
        else 
            A2 = A(ptr:end, ptr+block_size:end);
            [A(ptr:end,ptr:ptr+block_size-1), V, Tau] = simpleQR(A(ptr:end,ptr:ptr+block_size-1));
            S = CWYTransform(V, Tau);
            A(ptr:end, ptr+block_size:end) = A2 - V*(S*(V'*A2)) ;
            Q(ptr:end,:) = Q(ptr:end,:) - V*(S*(V'*Q(ptr:end,:))); 
        end 
    end
    Q = Q';
    R=A;
end