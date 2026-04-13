using KrylovMethods
using Helmholtz
using jInv.Mesh
using PyPlot

close("all")

include("MGsetup.jl")
include("MGcycle.jl")
include("getModels.jl")


function getAcousticHelmholtzMGVankaSetup(n,model,intergrid,relaxType,patch,levels,shift;order=4,source="mid",omega_factor=1)

    println("================= grid size is ",n," =================")

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
        m = getLinearModel(1,2,nodes); println("linear model")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            Omega = [0.0;1.0;0.0;1.0;0.0;1.0];
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "linear025"
        m = getLinearModel(0.25,1,nodes); println("linear model [0.25,1]")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            Omega = [0.0;1.0;0.0;1.0;0.0;1.0];
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "linear005"
        m = getLinearModel(0.05,1,nodes); println("linear model [0.05,1]")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            Omega = [0.0;1.0;0.0;1.0;0.0;1.0];
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "wedge025"
        m = getWedge(0.25,1,nodes); println("wedge model [0.25,1]")
        if dim == 2
            Omega = [0.0;1.0;0.0;1.0];
        elseif dim == 3
            println("wedge is 2D")
        end        
        M = getRegularMesh(Omega,n);
        h = M.h;
    elseif model == "wedge005"
        m = getWedge(0.05,1,nodes); println("wedge model [0.05,1]")
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
        println("Overthrust model")
        n = M.n;
        h = M.h;
        Omega = M.domain;
        nodes = n .+ 1;
        println("grid size (cells) after extension: ",n)
    end

    omega = omega_factor * getMaximalFrequency(m,M);
    println("omega is ",omega/pi," times pi")

    pad = 20;
    if dim == 2
        pad_vec = [pad;pad]
    elseif dim == 3
        pad_vec = [pad;pad;pad]
    end
    aten = 0.0;
    println("shift is ",shift)
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

    println(intergrid," intergrid")
    if relaxType == "Jacobi"
        println("Jacobi")
    elseif relaxType == "Vanka"
        println(patch," patch Vanka")
    end

    @time begin
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
    println("relaxParam = ",relaxParam)


    return H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam

end



function solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels,nu,relaxParam,relaxType; inner = 5, tol = 1e-6)

    nu1 = nu[1];
    nu2 = nu[2];
    println("nu is ",nu)
    println(levels," levels")
    if recursive_calls == 1
        println("V-cycle")
    elseif recursive_calls == 2
        println("W-cycle")
    end

    function PrecFuncSL(r) 

        e = MGcycle(H_s,r,0.0*r,relaxParam,nu1,nu2,levels,recursive_calls,R_arr,P_arr,Ac_arr,LUcoarsest; relaxType, M_arr);

        return e
    end


    @time begin
      e = fgmres(H, (1.0 + 0.0*im)*b, inner, maxIter = 200, M = PrecFuncSL, out = 2, tol = tol , flexible = true)[1];
    end

end






println("=================================================================")
println("                           2D experiments                        ") 
println("=================================================================")

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 3                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

n_arr = [[[128;128]];[[256;256]];[[512;512]]];
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]
plus = ["Plus";"Plus";"Plus";"Plus";"Plus";"Plus"]
rb = ["RB";"RB";"RB";"RB";"RB";"RB"]
patches_arr = [[jac];[element];[plus];[rb]];
shifts_arr = [0.3;0.25;0.25;0.18]
intergrid_arr = ["bicub";"mixed";"levdep"]
inners = [5;100]
levels = 4
recursive_calls = 2

for i=1:length(n_arr)
    for j=1:length(patches_arr)
        for k=1:length(intergrid_arr)
            for l=1:length(inners)
                println("GMRES(",inners[l],")")
                H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n_arr[i],"const",intergrid_arr[k],"Vanka",patches_arr[j],levels,shifts_arr[j])
                solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels,[1;1],relaxParam,"Vanka"; inner = inners[l])
            end
        end
    end
end

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                          Fig. 7                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

n = [256;256]
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
rb = ["RB";"RB";"RB";"RB";"RB";"RB";"RB"]
patches_arr = [[jac];[rb];[rb]];
nu_arr = [[[1;1]];[[1;1]];[[1;0]]]
levels_arr = [2;3;4;5;6;7]
Jacobi11_shift_arr = [0.0;0.1;0.7;2.5;2.5;2.5]
Vanka11_shift_arr = [0.0;0.1;0.18;0.18;0.18;0.18];
Vanka10_shift_arr = [0.0;0.08;0.15;0.15;0.15;0.15];
shifts_arr = [Jacobi11_shift_arr';Vanka11_shift_arr';Vanka10_shift_arr']
recursive_calls = 1

for i=1:length(levels_arr)
    for j=1:length(patches_arr)
        H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n,"const","levdep","Vanka",patches_arr[j],levels_arr[i],shifts_arr[j,i])
        solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels_arr[i],nu_arr[j],relaxParam,"Vanka")
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 4                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n_arr = [[[128;128]];[[256;256]];[[512;512]];[[1024;1024]]];
rb = ["RB";"RB";"RB";"RB";"RB";"RB"]
levels_arr = [2;3;4]
shifts_arr = [0.0 0.1 0.25;
              0.0 0.1 0.25;
              0.0 0.1 0.25;
              0.0 0.1 0.3]
models_arr = ["linear025";"wedge025";"linear005";"wedge005"]
recursive_calls = 1

for i=1:length(n_arr)
    for j=1:length(models_arr)
        for k=1:length(levels_arr)
            H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n_arr[i],models_arr[j],"levdep","Vanka",rb,levels_arr[k],shifts_arr[j,k])
            solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels_arr[k],[1;1],relaxParam,"Vanka")
        end
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                          Fig. 9                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n = [256;256];
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
rb = ["RB";"RB";"RB";"RB";"RB";"RB"]
patches_arr = [[jac];[rb]];
levels = 4
shifts_arr = [0.17;0.19;0.21;0.23;0.25;0.27]
models_arr = ["linear025";"wedge025";"linear005";"wedge005"]
recursive_calls = 1

for i=1:length(patches_arr)
    for j=1:length(models_arr)
        for k=1:length(shifts_arr)
            H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n,models_arr[j],"levdep","Vanka",patches_arr[i],levels,shifts_arr[k])
            solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels,[1;1],relaxParam,"Vanka")
        end
    end
end


println("=================================================================")
println("                           3D experiments                        ") 
println("=================================================================")


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 5                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

n_arr = [[[48;48;48]];[[64;64;64]];[[96;96;96]];[[128;128;128]]];
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]
plus = ["Plus";"Plus";"Plus";"Plus";"Plus";"Plus"]
patches_arr = [[jac];[element];[plus]];
shifts_arr = [0.5;0.4;0.65]
levels = 5
recursive_calls = 1

for i=1:length(n_arr)
    for j=1:length(patches_arr)
        H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n_arr[i],"const","levdep","Vanka",patches_arr[j],levels,shifts_arr[j])
        solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels,[1;1],relaxParam,"Vanka")
    end
end



println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                        Fig. 10a                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n_arr = [[[48;48;48]];[[64;64;64]];[[96;96;96]];[[128;128;128]]];
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]
patches_arr = [[jac];[element]];
intergrid_arr = ["trilin";"levdep"]
shift = 0.5
levels = 5
recursive_calls = 1

for i=1:length(n_arr)
    for j=1:length(patches_arr)
        H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n_arr[i],"linear",intergrid_arr[j],"Vanka",patches_arr[j],levels,shift)
        solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels,[1;1],relaxParam,"Vanka")
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                        Fig. 10b                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


n = [128;128;128];
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]
patches_arr = [[jac];[element]];
intergrid_arr = ["levdep";"levdep"]
levels_arr = [3;4;5;6;7]
shifts_arr = [0.2;0.5;0.5;0.5;0.5]
recursive_calls = 1

for i=1:length(levels_arr)
    for j=1:length(patches_arr)
        H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n,"linear",intergrid_arr[j],"Vanka",patches_arr[j],levels_arr[i],shifts_arr[i])
        solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels_arr[i],[1;1],relaxParam,"Vanka")
    end
end



println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                        Fig. 12                          ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

### need to hold a folder "BenchmarkModels" with the Overthrust file

n = [256;256;80]

levels_arr = [3;4;5;6];
shifts_arr = [0.2;0.4;0.5;0.5]
patch = ["Jac";"Element";"Element";"Element";"Element";"Element"]
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
element = ["Element";"Element";"Element";"Element";"Element";"Element"]
patches_arr = [[patch];[jac];[element]]

recursive_calls = 1

for i=1:length(patches_arr)
    for j=1:length(levels_arr)
        H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n,"Overthrust","trilin","Vanka",patches_arr[i],levels_arr[j],shifts_arr[j])
        solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels_arr[j],[1;1],relaxParam,"Vanka")
    end
end



println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 6                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

## need to hold a folder "BenchmarkModels" with the Overthrust file

n_arr = [[[128;128;40]];[[192;192;60]];[[256;256;80]];[[320;320;100]]];
levels_arr = [3;4];
shifts_arr = [0.2;0.4]
patch = ["Jac";"Element";"Element";"Element";"Element";"Element"]
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
patches_arr = [[patch];[jac]]

recursive_calls = 1

for i=1:length(n_arr)
    for j=1:length(levels_arr)
        for k=1:length(patches_arr)
            H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n_arr[i],"Overthrust","trilin","Vanka",patches_arr[k],levels_arr[j],shifts_arr[j])
            solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels_arr[j],[1;1],relaxParam,"Vanka")
        end
    end
end

n_arr = [[[256;256;80]];[[320;320;100]]];
levels_arr = [5];
shifts_arr = [0.5]
patch = ["Jac";"Element";"Element";"Element";"Element";"Element"]
jac = ["Jac";"Jac";"Jac";"Jac";"Jac";"Jac";"Jac"]
patches_arr = [[patch];[jac]]

recursive_calls = 1

for i=1:length(n_arr)
    for j=1:length(levels_arr)
        for k=1:length(patches_arr)
            H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,relaxParam = getAcousticHelmholtzMGVankaSetup(n_arr[i],"Overthrust","trilin","Vanka",patches_arr[k],levels_arr[j],shifts_arr[j])
            solveMGVanka(H,H_s,b,R_arr,P_arr,Ac_arr,LUcoarsest,M_arr,recursive_calls,levels_arr[j],[1;1],relaxParam,"Vanka")
        end
    end
end






println("all done!")