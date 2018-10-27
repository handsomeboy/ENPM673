function [Fnorm] = FindFundamentalMatrix(X,X_nxt)

for i = 1: size(X,1)
    
    xl = X(i, 1);
    yl = X(i, 2);
    xr = X_nxt(i, 1);
    yr = X_nxt(i, 2);
    A(i, :) = [xl * xr yl * xr xr xl * yr yl * yr yr xl yl 1]; 
end

[U, S, V] = svd(A);
temp = reshape(V(:,end), [3 3])';
[U, S, V] = svd(temp);
S(3,3) = 0;
S2 = S;

Fnorm = U * S2 * V';

end

