using KrylovMethods
using Helmholtz
using jInv.Mesh
using PyPlot

close("all")

include("MGsetup.jl")
include("MGcycle.jl")
include("MGsolver.jl")
include("getModels.jl")


function getAcousticHelmholtzMGVankaSetup(n,model,intergrid,relaxType,patch,levels,shift;order=4,source="mid",omega_factor=1)

    # println("================= grid size is ",n," =================")

    nodes = n .+ 1;
    dim = length(n)

    if model == "const"
        m = 1.0 * ones(nodes...);
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            Omega = [0.0;1.0;0.0;1.0;0.0;1.0];
        end
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "linear"
        m = getLinearModel(1,2,nodes); # println("linear model")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            Omega = [0.0;1.0;0.0;1.0;0.0;1.0];
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "linear025"
        m = getLinearModel(0.25,1,nodes); # println("linear model [0.25,1]")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            Omega = [0.0;1.0;0.0;1.0;0.0;1.0];
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "linear005"
        m = getLinearModel(0.05,1,nodes); # println("linear model [0.05,1]")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            Omega = [0.0;1.0;0.0;1.0;0.0;1.0];
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "wedge025"
        m = getWedge(0.25,1,nodes); # println("wedge model [0.25,1]")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            println("wedge is 2D")
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "wedge005"
        m = getWedge(0.05,1,nodes); # println("wedge model [0.05,1]")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            println("wedge is 2D")
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "Overthrust"
        _,_,_,M,Vp = getModel("OverthrustAcoustic",nodes); 
        m = 1 ./ (Vp .^ 2); 
        # println("Overthrust model")
        n = M.n;
        h = M.h;
        Omega = M.domain;
        nodes = n .+ 1;
        # println("grid size (cells) after extension: ",n)
    end

    omega = omega_factor * getMaximalFrequency(m,M);
    # println("omega is ",omega/pi," times pi")

    pad = 20;
    if dim == 2
        pad_vec = [pad;pad]
    elseif dim == 3
        pad_vec = [pad;pad;pad]
    end
    aten = 0.0;
    # println("shift is ",shift)
    alpha = aten + shift*omega;
    neumanOnTop = false;
    gamma = getABL(nodes, neumanOnTop, pad_vec, omega) .+ aten;
    gamma_s = getABL(nodes, neumanOnTop, pad_vec, omega) .+ alpha;
    param = HelmholtzParam(M,gamma,m,omega,neumanOnTop,false);
    param_s = HelmholtzParam(M,gamma_s,m,omega,neumanOnTop,false);
    if order == 4
        if dim == 2
            H = GetHelmholtzOperatorHO(param,[2/3;2/3]); 
            H_s = GetHelmholtzOperatorHO(param_s,[2/3;2/3]); 
        elseif dim == 3
            H = GetHelmholtzOperatorHO(param,[1/3;0.5]); 
            H_s = GetHelmholtzOperatorHO(param_s,[1/3;0.5]); 
        end
    elseif order == 2
        H = GetHelmholtzOperator(param); 
        H_s = GetHelmholtzOperator(param_s); 
    end

    # right hand side 
    if source == "mid"
        q,src = getAcousticPointSource(M,Float64,getMidPointSrc(M));
    elseif source == "top"
        q,src = getAcousticPointSource(M,Float64);
    end
    b = vec(q);

    #### multigrid preconditioner

    nodal = [true;true;true];

    if intergrid == "trilin" || intergrid == "bilin"
        intergridTypeArr = ["BI";"BI";"BI";"BI";"BI";"BI"];
    elseif intergrid == "tricub" || intergrid == "bicub"
        intergridTypeArr = ["high";"high";"high";"high";"high";"high"];
    elseif intergrid == "mixed"
        intergridTypeArr = ["mixed_high";"mixed_high";"mixed_high";"mixed_high";"mixed_high";"mixed_high"];
    elseif intergrid == "levdep"
        intergridTypeArr = ["high";"mixed_high";"mixed_high";"mixed_high";"mixed_high";"mixed_high";"mixed_high";"mixed_high"];
    end
    R_arr,P_arr,Ac_arr,LUcoarsest = myMGsetup(H_s,n,levels,nodal; intergridTypeArr);

    # println(intergrid," intergrid")
    # if relaxType == "Jacobi"
    #     println("Jacobi")
    # elseif relaxType == "Vanka"
    #     println(patch," patch Vanka")
    # end

    if relaxType == "Vanka"
        M_fine = VankaSetup(H_s,n,patch[1]) 
        M_arr = [M_fine];
        for i=1:levels-2
            M_arr_temp = VankaSetup(Ac_arr[i],div.(n,2^i),patch[i])
            M_arr = [M_arr; [M_arr_temp]];
        end
    elseif relaxType == "Jacobi"
        M_fine = 1 ./ diag(H_s);
        M_arr = [M_fine];
        for i=1:levels-2
            M_arr_temp = 1 ./ diag(Ac_arr[i])
            M_arr = [M_arr; [M_arr_temp]];
        end
    end


    if dim == 2
        relaxParamJac = [0.89;0.9;0.3;0.71;0.79;0.78]; 
        relaxParamElement = [0.97;0.66;0.48;0.88;0.88];
        relaxParamPlus = [0.87;0.57;0.55;0.74;1.09;1.1]; 
        relaxParamRB = [0.83;0.5;0.4;0.65;0.7;0.8];
    elseif dim == 3
        relaxParamJac = [0.6;0.4;0.3;0.5;0.5];
        relaxParamElement = [1.1;0.7;0.45;0.6;0.7;0.7];
        relaxParamPlus = [0.92;0.55;0.45;0.55];
        relaxParamRB = [1.0;1.0;1.0;1.0;1.0;1.0];
    end

    relaxParam = zeros(levels-1)
    if relaxType == "Jacobi"
        relaxParam = relaxParamJac;
    else
        for i=1:levels-1
            if patch[i] == "Element"
                relaxParam[i] = relaxParamElement[i];
            elseif patch[i] == "Plus"
                relaxParam[i] = relaxParamPlus[i];
            elseif patch[i] == "Jac"
                relaxParam[i] = relaxParamJac[i];
            elseif patch[i] == "RB"
                relaxParam[i] = relaxParamRB[i];
            end
        end
    end
    # println("relaxParam = ",relaxParam)


    return H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam

end



function solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels,nu,relaxParam,relaxType; inner = 5, tol = 1e-6)

    nu1 = nu[1];
    nu2 = nu[2];
    # println("nu is ",nu)
    # println(levels," levels")
    # if recursive_calls == 1
    #     println("V-cycle")
    # elseif recursive_calls == 2
    #     println("W-cycle")
    # end

    function PrecFuncSL(r) 

        e = MGcycle(H_s,r,0.0*r,relaxParam,nu1,nu2,levels,recursive_calls,R_arr,P_arr,Ac_arr,LUcoarsest; relaxType, M_arr);

        return e
    end


    elapsed_time = @elapsed begin
      x, flag, err, iter, resvec = fgmres(H, (1.0 + 0.0*im)*b, inner, maxIter = 200, M = PrecFuncSL, out = -2, tol = tol , flexible = true);
    end

    return length(resvec), elapsed_time

end