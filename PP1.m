function PP1
    % The exact solution to the heat equation
    u = @(x,t) exp(-t).* sin(x);
    % u(x,0)
    f = @(x) sin(x);
    
    % Scheme parameters
    N = 10; L = pi; dx = L/N; dt = 0.01; tmax = 100; x = dx * (1:N-1);
    

    function S1_U = generateS1()
        % The iterative matrix for generating S1
        A = eye(N-1) + dt/(dx^2)*(-2*eye(N-1)+diag(ones(N-2,1),1) + diag(ones(N-2,1),-1));
        
        S1_U(1,:) = f(x);
        hold on;
        for i = 2:tmax
            S1_U(i,:) = S1_U(i-1,:) * A;
        end
    end
end