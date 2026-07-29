function SmootherSymbol = AdditiveVanka(V,W,invHi,h,sigma,theta,w,smoother,stencil)

    theta1 = theta(1);
    theta2 = theta(2);

    if smoother == "ElementVanka"
        Phi = diag([1, exp(1i * theta1), exp(1i * theta2), exp(1i * (theta1+theta2))]);
    elseif smoother == "PlusVanka"
        Phi = diag([exp(-1i * theta2), exp(-1i * theta1), 1, exp(1i * theta1), exp(1i * theta2)]);
    elseif smoother == "RBVanka"
            Phi = diag([exp(-1i * (theta1 + theta2)), exp(-1i * (theta2 - theta1)), 1, exp(1i * (theta2 - theta1)), exp(1i * (theta1 + theta2))]);
    elseif smoother == "FullVanka"
        Phi = exp(-1i * (theta1 + theta2));
        Phi = [Phi exp(-1i * theta1)];
        Phi = [Phi exp(-1i * (theta1 - theta2))];
        Phi = [Phi exp(-1i * theta2)];
        Phi = [Phi 1];
        Phi = [Phi exp(1i * theta2)];
        Phi = [Phi exp(1i * (theta1 - theta2))];
        Phi = [Phi exp(1i * theta1)];
        Phi = [Phi exp(1i * (theta1 + theta2))];
        Phi = diag(Phi);
    elseif smoother == "Jacobi"
        Phi = 1;
    end

    if stencil == 5
        Hsymbol = (1/h^2)*(4-(h^2)*sigma-2*cos(theta1)-2*cos(theta2));
    elseif stencil == 9
        Asymbol = (1/(3*h^2))*(10-4*cos(theta1)-4*cos(theta2)-2*cos(theta1)*cos(theta2));
        Msymbol = (sigma/6)*(4 + cos(theta1) + cos(theta2));
        Hsymbol = Asymbol - Msymbol;
    end

    M = ((V'*W)*Phi')*(invHi*(Phi*V));
    f = M*Hsymbol;

    SmootherSymbol = 1 - w * f;

end