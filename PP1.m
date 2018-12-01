function PP1
    % The exact solution to the heat equation
    u = @(x,t) exp(-t).* sin(x);
    % u(x,0)
    f = @(x) sin(x);
    
    % Scheme parameters
    N = 10; L = pi; dx = L/N; dt = 0.01; tmax = 1/dt;
    
    % The initial value for the schemes
    U_0 = f(dx * (1:N-1));
    
    S1_U = generateS1(U_0,dx,dt,N,tmax);
    S2_U = generateS2(U_0,dx,dt,N,tmax);

    E1 = calculateError(S1_U,u,dx,dt,N,tmax);
    E2 = calculateError(S2_U,u,dx,dt,N,tmax);
    

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
            diff = u(dx*(1:xmax-1),i*dt) - S1_U(i,:);
            E(i) = norm(diff,inf);
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
        E = zeros(1,100);
        
        for i = 1:tmax
            diff = u(dx*(1:xmax-1),i*dt) - SU(i,:);
            E(i) = norm(diff,inf);
        end
    end
end