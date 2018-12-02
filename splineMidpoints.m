function set = splineMidpoints(xn,yn)
    S = cubicSpline(xn,yn);
    set = [];

    % yn(1) == 0 if we're dealing with dx, we want to add the first
    % position if we're dealing with dt
    if (yn(1) ~= 0)
        set = [set,yn(1)];
    end
    for alpha = 1:length(xn)-1
        % the average value, we want approximations of the midpoint from
        % the phi(h) estimation
        chi = (xn(alpha)+xn(alpha+1))/2;
        % x-xi
        delta = chi - xn(alpha);

        gamma = S(alpha,:) * [delta^3;delta^2;delta;1];

        % Inserting the non-edge values


        set = [set,gamma];
        if (alpha ~= length(xn)-1)
            set = [set,yn(alpha+1)];
        end
    end
    % yn(end) == 0 if we're dealing with dx, we want to add the first
    % position if we're dealing with dt
    if (yn(end) ~= 0)
        set = [set,yn(end)];
    end
end