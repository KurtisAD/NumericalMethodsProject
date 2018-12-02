function S2R = coarseRefined(S2,S2_0,dx,dt,xmax,tmax)
    % Appending the base condition to spline between
    S2 = [S2_0;S2];
    [M N] = size(S2);    
    
    % We need an intial condition between t=0 and dt
    for tau = 1:M
        S2R_0(tau,:) = splineMidpoints(dx*(0:xmax),[0,S2(tau,:),0]);
    end
    [M N] = size(S2R_0);
    
    for chi = 1:N
        S2R(:,chi) = splineMidpoints(dt*(0:tmax),S2R_0(:,chi)')';
    end
    S2R = S2R(2:end,:);
end