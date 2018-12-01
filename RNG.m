function x = RNG(A)
%
% The forward elimination is attained by the following loops
%
[M N] = size(A);
%
for n=1:N-2
    for m=n+1:M
        A(m,:) = A(m,:) - A(m,n)/A(n,n)*A(n,:);
    end
end
% and the backward substitution by the following loop
x(M) = A(M,N)/A(M,M);
for m=M-1:-1:1
    x(m) = ( A(m,N) - sum( A(m,m+1:M).*x(m+1:M) ) )/A(m,m);
end
