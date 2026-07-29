include("myVanka.jl")

function MGcycle(A,b,x,w,nu1,nu2,levels,recursive_calls,R_arr,P_arr,Ac_arr,LUAcoarsest; coarseSolve = "exact", relaxType = "Jacobi", M_arr = [])

	# for Vcycle, recurcive_calls = 1
    # for Wcycle, recursive_calls = 2
    # for two-grid, levels = 2

	if levels == 1
        if coarseSolve == "exact"
		    return LUAcoarsest \ b
        elseif coarseSolve == "gmres"
            # println("coarse solve iterations")
            inner = 10;
            e = fgmres(A, b, inner, maxIter = 15, out = 1, tol = 0.01 , flexible = true)[1]
            # println("finished coarse solve iterations")
            return e
        elseif coarseSolve == "jac"
            tol = 1e-6;
            e,_ = dampedJac(A,b,0.2,(0.0 + 0.0*1im)*zeros(size(b)),tol,10);
            return e
        elseif coarseSolve == "GS"
            tol = 1e-6;
            e,_ = GaussSeidel(A,b,(0.0 + 0.0*1im)*zeros(size(b)),tol,5);
            return e
        end
	end
	
	# pre-smoothing
    tol = 1e-5;
    x = Vanka(A,b,w[1],x,M_arr[1],tol,nu1);



    # compute and restrict the residual
    r = b - A * x;

    R = R_arr[1];
    P = P_arr[1];
    Ac = Ac_arr[1]

    R_arr = R_arr[2:end];
    P_arr = P_arr[2:end];
    Ac_arr = Ac_arr[2:end];

    rc = R * r;

    # solve the error-residual equation directly or recursively
    ec = 0.0 .* rc; # initial guess
	if levels == 2
		recursive_calls = 1;
	end
    for j=1:recursive_calls
        ec = MGcycle(Ac,rc,ec,w[2:end],nu1,nu2,levels-1,recursive_calls,R_arr,P_arr,Ac_arr,LUAcoarsest ; coarseSolve, relaxType, M_arr = M_arr[2:end]);
    end

    e = P * ec;
    x = x + e;

    x = Vanka(A,b,w[1],x,M_arr[1],tol,nu2);

    return x

end
