function PP2
    % The exact solution to the heat equation
    u = @(x,t) exp(-t).* sin(x);
    % u(x,0)
    f = @(x) sin(x);
    
    % Scheme parameters
    N = 10; L = pi; dx = L/N; dt = 0.1; tmax = 1/dt;
    
    % The initial value for the schemes
    U_0 = f(dx * (1:N-1));
    U2_0 = f(dx/2 * (1:2*N-1));
    
    % First derivation
    Phih = generateS2(U_0,dx,dt,N,tmax);
    Phih2 = generateS2(U2_0,dx/2,dt/2,2*N,2*tmax);
    PhihR = coarseRefined(Phih,U_0,dx,dt,N,tmax);
    SR = (4*Phih2-PhihR)/3;
    
    % Errors
    E1 = calculateError(Phih,u,dx,dt,N,tmax);
    E2 = calculateError(Phih2,u,dx/2,dt/2,2*N,2*tmax-1);
    E3 = calculateError(PhihR,u,dx/2,dt/2,2*N,2*tmax-1);
    E4 = calculateError(SR,u,dx/2,dt/2,2*N,2*tmax-1);
    
    
    hold on;
    plot(1:tmax,E1);
    plot((1:(2*tmax-1))/2,E2);
    plot((1:(2*tmax-1))/2,E3);
    plot((1:(2*tmax-1))/2,E4);
    
    legend('phi(h)','phi(h/2)','refined phi(h)','SR');
    

end