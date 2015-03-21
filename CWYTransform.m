function S = CWYTransform(V, Tau)
    col = size(V,2);
    S = zeros(col, col);
    S(1,1) = 1/Tau(1);
    for k = 2:col,
        S00 = S(1:k-1,1:k-1);
        s01 = V(:,1:k-1)'*V(:,k);
        sig11 = 1/Tau(k);
        s01 = -S00*s01*sig11;

        S(1:k-1,k) = s01;
        S(k,k) = sig11;
    end
end