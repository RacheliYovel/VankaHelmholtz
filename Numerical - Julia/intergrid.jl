function intergrid1D(n,isnodal; intergridType = "BI") # Neumann, put the number of cells as n, isnodal is a scalar
    
    if isnodal == true

        P_BI = spdiagm(-1 => 0.5 .* ones(n), 0 => ones(n+1), 1 => 0.5 .* ones(n));
        P_HO = spdiagm(-2 => 0.125 .* ones(n-1), -1 => 0.5 .* ones(n), 0 => 0.75 .* ones(n+1), 1 => 0.5 .* ones(n), 2 => 0.125 .* ones(n-1));
        

        if intergridType == "mixed_high"
            P = P_HO[:,1:2:end];
            R = 0.5 .* P_BI[:,1:2:end]';
        elseif intergridType == "BI"
            P = P_BI[:,1:2:end];
            R = 0.5 .* P';
        elseif intergridType == "high"
            P = P_HO[:,1:2:end];
            R = 0.5 .* P';
        end

        R[1,:] = R[1,:] ./ sum(R[1,:]);
        R[end,:] = R[end,:] ./ sum(R[end,:]);
        P[1,:] = P[1,:] ./ sum(P[1,:]);
        P[end,:] = P[end,:] ./ sum(P[end,:]);


    else

        P = spdiagm(-1 => 0.25 .* ones(n-1), 0 => 0.75 .* ones(n), 1 => 0.75 .* ones(n-1), 2 => 0.25 .* ones(n-2));
        P = P[:,2:2:end];

        P[1,1] = 1;
        P[end,end] = 1;

        if intergridType == "mixed"
            R = spdiagm(0 => 0.5 .* ones(n), 1 => 0.5 .* ones(n-1));
            R = R[:,2:2:end];
            R = 1.0 .* R';
        else
            R = 0.5 .* P';
        end

    end
    
    return R,P
end


function intergrid(n,isnodal ; intergridType = "BI") # Neumann, put an array of number of cells as n, isnodal is an array

    dim = length(n);

    R,P = intergrid1D(n[1],isnodal[1] ; intergridType);
    R = [R];
    P = [P];
    for i=2:dim
        temp,_ = intergrid1D(n[i],isnodal[i] ; intergridType);
        R = [[temp] ; R];
        _,temp = intergrid1D(n[i],isnodal[i] ; intergridType);
        P = [[temp] ; P];
    end
    if dim > 1
        R = kron(R...);
        P = kron(P...);
    else
        R = R[1];
        P = P[1];
    end
    
    return R,P
end
