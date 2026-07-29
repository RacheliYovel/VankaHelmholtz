function symbol = symbol5by5(s0,s1,s2,theta)
    % the stencil is symmetric and is a Kroncker product of 1D stencils
    % [s2 s1 s0 s1 s2]
    % theta is a vector (theta1;theta2)

    symbol1 = (s0 + 2*s1*cos(theta(1)) + 2*s2*cos(2*theta(1)));
    symbol2 = (s0 + 2*s1*cos(theta(2)) + 2*s2*cos(2*theta(2)));
    symbol = symbol1 * symbol2;

end