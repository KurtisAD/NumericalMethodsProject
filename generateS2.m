function S2_U = generateS2(S2_0,dx,dt,xmax,tmax)
    % The left hand side iterative matrix
    A_L = eye(xmax-1)-dt/(2*dx^2)*(-2*eye(xmax-1)+diag(ones(xmax-2,1),1) + diag(ones(xmax-2,1),-1));
    % The right hand side
    A_R = eye(xmax-1)+dt/(2*dx^2)*(-2*eye(xmax-1)+diag(ones(xmax-2,1),1) + diag(ones(xmax-2,1),-1));

    S2_U = zeros(tmax,xmax-1);
    matrix = [A_L,(S2_0*A_R)'];
    S2_U(1,:) = RNG(matrix);

    for i = 2:tmax
        matrix = [A_L,(S2_U(i-1,:)*A_R)'];
        S2_U(i,:) = RNG(matrix);
    end
end
