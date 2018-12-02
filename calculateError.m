function E = calculateError(SU,u,dx,dt,xmax,tmax)
    % t is bounded between 0:tmax, but at t=0 the error is 0
    % x is bounded between 1:xmax-1

    % initializing
    E = zeros(1,tmax);

    for i = 1:tmax
        diff = u(dx*(1:xmax-1),(i)*dt) - SU(i,:);
        E(i) = norm(diff,inf);
    end
end