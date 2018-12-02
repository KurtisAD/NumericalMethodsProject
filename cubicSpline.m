function S = cubicSpline(xn,yn)
    % construct the matrix for zn
    % z1 = zn+1 = 0
    hn = xn(2:end) - xn(1:end-1);
    wn = (yn(2:end) - yn(1:end-1))./hn;

    A = diag(hn(2:end-1),1) + diag(hn(2:end-1),-1) + 2.*(diag(hn(1:end-1) + hn(2:end)));
    A = [A,6.*(wn(2:end)-wn(1:end-1))'];
    zn = [0,RNG(A),0];

    an = (zn(2:end) - zn(1:end-1))./(6.*hn);
    bn = zn./2;
    cn = wn-hn./6.*(zn(2:end) + 2.*zn(1:end-1));
    S = [an',bn(1:end-1)',cn',yn(1:end-1)'];
end