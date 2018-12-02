function PP2
    % The exact solution to the heat equation
    u = @(x,t) exp(-t).* sin(x);
    % u(x,0)
    f = @(x) sin(x);
    
    % Scheme parameters
    N = 10; L = pi; dx = L/N; dt = 0.1; tmax = 1/dt+1;
    
    % The initial value for the schemes
    U_0 = f(dx * (1:N-1));

end