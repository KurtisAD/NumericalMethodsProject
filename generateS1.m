function S1_U = generateS1(S1_0,dx,dt,xmax,tmax)
    % The iterative matrix for generating S1
    A = eye(xmax-1) + dt/(dx^2)*(-2*eye(xmax-1)+diag(ones(xmax-2,1),1) + diag(ones(xmax-2,1),-1));

    % Initilize S1_U
    S1_U = zeros(tmax,xmax-1);
    S1_U(1,:) = S1_0*A;

    for i = 2:tmax
        S1_U(i,:) = S1_U(i-1,:) * A;
    end
end
