function PP1
    % The exact solution to the heat equation
    u = @(x,t) exp(-t).* sin(x);
    % u(x,0)
    f = @(x) sin(x);
    
    % Scheme parameters
    N = 10; L = pi; dx = L/N; dt = 0.1; tmax = 1/dt+1;
    
    % The initial value for the schemes
    U_0 = f(dx * (1:N-1));
    
    S1_U = generateS1(U_0,dx,dt,N,tmax);
    S2_U = generateS2(U_0,dx,dt,N,tmax);

    E1 = calculateError(S1_U,u,dx,dt,N,tmax);
    E2 = calculateError(S2_U,u,dx,dt,N,tmax);
    
    % Plotting
    % We need to plot for timeslots 0.1, 0.5 and 1
    chi = dx * (0:N);
    
    L1 = plot(chi,[0,U_0(1,:),0],'k');

    
    % Scheme plot
    figure(1);
    hold on;
    title('Scheme plot at t=0, 0.1, 0.5, and 1.0');
    ylabel('Heat');
    xlabel('x position');
    % S1
    L2 = plot(chi,[0,S1_U(0.1/dt+1,:),0],'b');
    plot(chi,[0,S1_U(0.5/dt+1,:),0],'b');
    plot(chi,[0,S1_U(1/dt+1,:),0],'b');
    
    % S2
    L3 = plot(chi,[0,S2_U(0.1/dt+1,:),0],'r');
    plot(chi,[0,S2_U(0.5/dt+1,:),0],'r');
    plot(chi,[0,S2_U(1/dt+1,:),0],'r');
    hold off;
    
    legend([L1,L2,L3],'f(x)','S1','S2');
    % Error plot
    figure(2);
    hold on;
    title('Error plot');
    ylabel('Average error');
    xlabel('Time');
    semilogy(dt * (0:tmax-1),E1, '--b');
    semilogy(dt * (0:tmax-1),E2, '--r');
    legend('S1','S2');
    hold off;
    
    %
    % Function Definitions
    %
    function S1_U = generateS1(S1_0,dx,dt,xmax,tmax)
        % The iterative matrix for generating S1
        A = eye(xmax-1) + dt/(dx^2)*(-2*eye(xmax-1)+diag(ones(xmax-2,1),1) + diag(ones(xmax-2,1),-1));
        
        % Initilize S1_U
        S1_U = zeros(tmax,xmax-1);
        S1_U(1,:) = S1_0;
        
        for i = 2:tmax
            S1_U(i,:) = S1_U(i-1,:) * A;
        end
    end
    function S2_U = generateS2(S2_0,dx,dt,xmax,tmax)
        % The left hand side iterative matrix
        A_L = eye(xmax-1)-dt/(2*dx^2)*(-2*eye(xmax-1)+diag(ones(xmax-2,1),1) + diag(ones(xmax-2,1),-1));
        % The right hand side
        A_R = eye(xmax-1)+dt/(2*dx^2)*(-2*eye(xmax-1)+diag(ones(xmax-2,1),1) + diag(ones(xmax-2,1),-1));
        
        S2_U = zeros(tmax,xmax-1);
        S2_U(1,:) = S2_0;
        
        for i = 2:tmax
            matrix = [A_L,(S2_U(i-1,:)*A_R)'];
            S2_U(i,:) = RNG(matrix);
        end
    end
    function E = calculateError(SU,u,dx,dt,xmax,tmax)
        % t is bounded between 0:tmax, but at t=0 the error is 0
        % x is bounded between 1:xmax-1
        
        % initializing
        E = zeros(1,tmax);
        
        for i = 1:tmax
            diff = u(dx*(1:xmax-1),(i-1)*dt) - SU(i,:);
            E(i) = norm(diff,inf);
        end
    end
end