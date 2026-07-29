include("intergrid.jl")

function myMGsetup(A,n,levels,isnodal;intergridTypeArr,coarseSolve = "exact")

    n_arr = (0.5 .^ (0:levels-1))' * 1.0 .* n;
    n_arr = Int.(n_arr);


    R1,P1 = intergrid(n,isnodal ; intergridType = intergridTypeArr[1]);
    Ac1 = R1 * A * P1;

    R_arr = [R1];
    P_arr = [P1];
    Ac_arr = [Ac1];
    for i=2:levels-1
        R_temp,P_temp = intergrid(n_arr[:,i],isnodal ; intergridType = intergridTypeArr[i]);
        R_arr = [R_arr ; [R_temp]];
        P_arr = [P_arr ; [P_temp]];
        
        Ac_temp = R_temp * Ac_arr[i-1] * P_temp;
        Ac_arr = [Ac_arr ; [Ac_temp]];
    end
    
    if coarseSolve == "exact"
        LUAcoarsest = lu(Ac_arr[end]);
    else 
        LUAcoarsest = 0.0;
    end

    return R_arr,P_arr,Ac_arr,LUAcoarsest

end