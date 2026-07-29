include("MGcycle.jl")

function MGsolver(A,b,guess,w,nu1,nu2,levels,recursive_calls,R_arr,P_arr,Ac_arr,LUAcoarsest,maxit,tol; relaxType = "Jacobi", M_arr = [], printres = true)

    x = copy(guess);
    r = b - A * x;
    initialResNorm = norm(r); relres = 1;
    iter = 0;
    r_arr = [norm(r)/initialResNorm];
    x_arr = [x];

    while relres > tol && iter < maxit
        if printres == true
            print(relres,", ")
        end
        x = MGcycle(A,b,x,w,nu1,nu2,levels,recursive_calls,R_arr,P_arr,Ac_arr,LUAcoarsest; relaxType, M_arr);
        iter = iter + 1;
        relres = norm(b - A*x)/initialResNorm;
        r = b - A * x;
        r_arr = [r_arr norm(r)/initialResNorm];
        temp1 = copy(x_arr);
        temp2 = copy(x);
        x_arr = [temp1 [temp2]];
        ######################################
        # x = reshape(x,size(m));

        # local z = zeros(ComplexF64,size(x));
        # z[2:end-1,2:end-1] = x[2:end-1,2:end-1];
        # x = z;

        # x[1,2:end-1] .= 0.0; # left
        # x[end,2:end-1] .= 0.0; # right
        # x[2:end-1,1] .= 0.0; # top
        # x[2:end-1,end] .= 0.0; # bottom

        # x[1,1] = 0.0; # left top corner
        # x[1,end] = 0.0; # left bottom corner
        # x[end,1] = 0.0; # right top corner
        # x[end,end] = 0.0; # right bottom corner

        x = vec(x);

        # zero right left top bottom separately
        # take only one size and fix it untill it works
        # then zero each corner and treat each one separately
        ######################################
    end

    r_arr = vec(r_arr);

    r_arr = r_arr[2:end];
    x_arr = x_arr[2:end];
    
    return x, iter, r_arr, x_arr

end