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
    
    figure(1);
    hold on;
    title('Scheme plot at t=0, 0.1, 0.5, and 1.0');
    ylabel('Heat');
    xlabel('x position');
    chi = dx/2 * (0:2*N);
    % Initial value
    L1 = plot(chi,[0,U2_0(1,:),0],'k');
    % PhiR(h)
    L2 = plot(chi,[0,PhihR(0.2/dt,:),0],'b');
    plot(chi,[0,PhihR(1/dt,:),0],'b');
    plot(chi,[0,PhihR(2/dt,:),0],'b');
    % Phi(h/2)
    L3 = plot(chi,[0,Phih2(0.2/dt,:),0],'r');
    plot(chi,[0,Phih2(1/dt,:),0],'r');
    plot(chi,[0,Phih2(2/dt,:),0],'r');
    % SR
    L4 = plot(chi,[0,SR(0.2/dt,:),0],'g');
    plot(chi,[0,SR(1/dt,:),0],'g');
    plot(chi,[0,SR(2/dt,:),0],'g');
    legend([L1,L2,L3,L4],'f(x)','Phi(h) Refined','Phi(h/2)','SR');
    hold off;
    
    % Error plot
    figure(2);
    hold on;
    plot(1:tmax,E1);
    plot((1:(2*tmax-1))/2,E2);
    plot((1:(2*tmax-1))/2,E3);
    plot((1:(2*tmax-1))/2,E4);
    
    legend('phi(h)','phi(h/2)','refined phi(h)','SR');
    set(gca,'yscale','log')

end