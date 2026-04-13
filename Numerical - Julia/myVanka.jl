# using StaticArrays

eye(n) = Matrix{Float64}(I, n, n)


function getDenseBlockFromAT(AT::SparseMatrixCSC{VAL,IND},Idxs::Array,Acc::Array{VAL,2}) where {VAL,IND}
    Acc[:] .= 0.0;
    for t = 1:length(Idxs)
        ii = AT.colptr[Idxs[t]];
        jj = 1;		
        while ii < AT.colptr[Idxs[t]+1] && jj <= length(Idxs)  
            if AT.rowval[ii] ==  Idxs[jj]
                Acc[t,jj] = conj(AT.nzval[ii]);
                ii+=1;
                jj+=1;
            elseif AT.rowval[ii] >  Idxs[jj]
                jj+=1;
            else
                ii+=1;
            end
        end
    end
    return Acc;
end


function VankaSetup(A,n,patch)

    nodes = n .+ 1;
    t = 1:prod(nodes);
    t = reshape(t,nodes...);
    
    if length(n) == 2
        if patch == "Jac"

            indices = vec(t);
            block_id = [0];
            W_grid = ones(size(t))
            
        elseif patch == "Element" # 4 points

            indices = vec(t[1:end-1,1:end-1]); # indices for upper left node of each 4-patch (cell)
            block_id = [0;1;nodes[1];nodes[1]+1]; # indices to add for each upper left index to get the entire cell

            W_grid = 2 * ones(size(t)); # number of patches that include the points
            W_grid[2:end-1,2:end-1] .= 4; 
            W_grid[1,1] = 1; W_grid[1,end] = 1; W_grid[end,1] = 1; W_grid[end,end] = 1;

        elseif patch == "Plus" # 5 points

            indices = vec(t[2:end-1,2:end-1]);; # indices for central node of each inner 5-patch
            block_id = [-nodes[1];-1;0;1;nodes[1]];

            W_grid = 4*ones(size(t)); # number of patches that include the node
            W_grid[2:end-1,2:end-1] .= 5; 
            W_grid[1,1] = 3; W_grid[1,end] = 3; W_grid[end,1] = 3; W_grid[end,end] = 3; 


        elseif patch == "RB" # 5 points rotated

            indices = vec(t[2:end-1,2:end-1]);; # indices for central node of each inner 5-patch
            block_id = [-nodes[1]-1;-nodes[1]+1;0;nodes[1]-1;nodes[1]+1];

            W_grid = 3 * ones(size(t)); # number of patches that include the points
            W_grid[2:end-1,2:end-1] .= 5; 
            W_grid[1,1] = 2; W_grid[1,end] = 2; W_grid[end,1] = 2; W_grid[end,end] = 2;

        elseif patch == "Full" # 9 points

            indices = vec(t[2:end-1,2:end-1]); # indices for central node of each 9-patch 
            block_id = [-nodes[1]-1;-nodes[1];-nodes[1]+1;-1;0;1;nodes[1]-1;nodes[1];nodes[1]+1];

            W_grid = 6 * ones(size(t)); # number of patches that include the points
            W_grid[2:end-1,2:end-1] .= 9; 
            W_grid[1,1] = 4; W_grid[1,end] = 4; W_grid[end,1] = 4; W_grid[end,end] = 4;

        elseif patch == "RB_wide" # 9 points streched rotated

            indices = vec(t[3:end-2,3:end-2]); # indices for central node of each 9-patch 
            block_id = [-2*nodes[1];-nodes[1]-1;-nodes[1]+1;-2;0;2;nodes[1]-1;nodes[1]+1;2*nodes]; 

            W_grid = ones(size(t)); # number of patches that include the points
            W_grid[1,2] = 2.0; W_grid[2,1] = 2.0; W_grid[2,end] = 2.0; W_grid[end,2] = 2.0; 
            W_grid[end,end-1] = 2.0; W_grid[end-1,end] = 2.0; W_grid[end-1,1] = 2.0; W_grid[1,end-1] = 2.0;
            W_grid[1,3:end-2] .= 3.0; W_grid[end,3:end-2] .= 3.0; W_grid[3:end-2,1] .= 3.0; W_grid[3:end-2,end] .= 3.0;
            W_grid[2,2] = 4.0; W_grid[2,end-1] = 4.0; W_grid[end-1,2] = 4.0; W_grid[end-1,end-1] = 4.0;
            W_grid[2,3] = 5.0; W_grid[3,2] = 5.0; W_grid[2,end-2] = 5.0; W_grid[end-2,2] = 5.0;
            W_grid[3,end-1] = 5.0; W_grid[end-1,3] = 5.0; W_grid[end-1,end-2] = 5.0; W_grid[end-2,end-1] = 5.0;
            W_grid[2,4:end-3] .= 6.0; W_grid[end-1,4:end-3] .= 6.0; W_grid[4:end-3,2] .= 6.0; W_grid[4:end-3,end-1] .= 6.0;
            W_grid[3,3] = 7.0; W_grid[3,end-2] = 7.0; W_grid[end-2,3] = 7.0; W_grid[end-2,end-2] = 7.0;
            W_grid[3,4:end-3] .= 8.0; W_grid[4:end-3,3] .= 8.0; W_grid[end-2,4:end-3] .= 8.0; W_grid[4:end-3,end-2] .= 8.0;
            W_grid[4:end-3,4:end-3] .= 9.0;

        end

    elseif length(n) == 3

        if patch == "Jac"

            indices = vec(t);
            block_id = [0];
            W_grid = ones(size(t))

        elseif patch == "Plus" # 7 points standard vertex wise
            indices = vec(t[2:end-1,2:end-1,2:end-1]); # indices for central node of each 7-patch 
            block_id = [-nodes[1]*nodes[2];-nodes[1];-1;0;1;nodes[1];nodes[1]*nodes[2]];

            W_grid = 5*ones(size(t)); # 12 edges (the rest will be calculated next)
            W_grid[2:end-1,2:end-1,2:end-1] .= 7; # 1 inner
            W_grid[1,2:end-1,2:end-1] .= 6; W_grid[end,2:end-1,2:end-1] .= 6; W_grid[2:end-1,1,2:end-1] .= 6;
            W_grid[2:end-1,end,2:end-1] .= 6; W_grid[2:end-1,2:end-1,1] .= 6; W_grid[2:end-1,2:end-1,end] .= 6; # 6 faces
            W_grid[1,1,1] = 4; W_grid[1,1,end] = 4; W_grid[end,1,1] = 4; W_grid[1,end,1] = 4;
            W_grid[end,end,end] = 4; W_grid[end,end,1] = 4; W_grid[1,end,end] = 4; W_grid[end,1,end] = 4; # 8 corners


        elseif patch == "Element" # 8 points (cell)

            indices = vec(t[1:end-1,1:end-1,1:end-1]); # indices for upper left node of each 8-patch (cell)
            block_id = [0;1;nodes[1];nodes[1]+1;nodes[1]*nodes[2];nodes[1]*nodes[2]+1;nodes[1]*nodes[2]+nodes[1];nodes[1]*nodes[2]+nodes[1]+1]; # indices to add for each upper left index to get the entire cell

            W_grid = 2*ones(size(t)); # 12 edges
            W_grid[2:end-1,2:end-1,2:end-1] .= 8; # 1 inner
            W_grid[1,2:end-1,2:end-1] .= 4; W_grid[end,2:end-1,2:end-1] .= 4; W_grid[2:end-1,1,2:end-1] .= 4;
            W_grid[2:end-1,end,2:end-1] .= 4; W_grid[2:end-1,2:end-1,1] .= 4; W_grid[2:end-1,2:end-1,end] .= 4; # 6 faces
            W_grid[1,1,1] = 1; W_grid[1,1,end] = 1; W_grid[end,1,1] = 1; W_grid[1,end,1] = 1;
            W_grid[end,end,end] = 1; W_grid[end,end,1] = 1; W_grid[1,end,end] = 1; W_grid[end,1,end] = 1; # 8 corners


        elseif patch == "Corners" # 9 points: corners and middle

            indices = vec(t[2:end-1,2:end-1,2:end-1]); # indices for central node of each 9-patch 
            block_id = []


        elseif patch == "RB" # 13 points

            indices = vec(t[2:end-1,2:end-1,2:end-1]); # indices for central node of each 13-patch 
            block_id = [-nodes[1]-nodes[2]*nodes[1]; -1-nodes[2]*nodes[1]; 1-nodes[2]*nodes[1]; nodes[1]-nodes[2]*nodes[1]; -nodes[1]-1; -nodes[1]+1; 0; nodes[1]-1; nodes[1]+1; -nodes[1]+nodes[2]*nodes[1]; -1+nodes[2]*nodes[1]; 1+nodes[2]*nodes[1]; nodes[1]+nodes[2]*nodes[1]];
            
            W_grid = 6*ones(size(t)); # 12 edges (the rest will be calculated next)
            W_grid[2:end-1,2:end-1,2:end-1] .= 13; # 1 inner
            W_grid[1,2:end-1,2:end-1] .= 9; W_grid[end,2:end-1,2:end-1] .= 9; W_grid[2:end-1,1,2:end-1] .= 9;
            W_grid[2:end-1,end,2:end-1] .= 9; W_grid[2:end-1,2:end-1,1] .= 9; W_grid[2:end-1,2:end-1,end] .= 9; # 6 faces
            W_grid[1,1,1] = 4; W_grid[1,1,end] = 4; W_grid[end,1,1] = 4; W_grid[1,end,1] = 4;
            W_grid[end,end,end] = 4; W_grid[end,end,1] = 4; W_grid[1,end,end] = 4; W_grid[end,1,end] = 4; # 8 corners

        end
    end

    W_grid = 1 ./ W_grid; # a grid matrix with the weights
    W_grid = vec(W_grid);


    if patch == "Element"
        nnz_max = length(indices)*length(block_id)^2
    else
        nnz_max = length(t)*length(block_id)^2
    end
    rows = Vector{Int}(undef, nnz_max)
    cols = Vector{Int}(undef, nnz_max)
    vals = Vector{ComplexF64}(undef, nnz_max)
    p = 1

    if patch == "Element"
        if length(n) == 2
            ind = [1;2;3;4];
        elseif length(n) == 3
            ind = [1;2;3;4;5;6;7;8];
        end
        running_index = indices
    elseif patch == "Jac"
        ind = [1]
        running_index = indices
    else
        running_index = vec(t)
    end

    for i in running_index
        # indices of the block
        if length(n) == 2 # 2D
            if patch == "Plus"
                if i in indices
                    ind = [1;2;3;4;5];
                elseif i in t[1,2:end-1]
                    ind = [1;3;4;5];
                elseif i in t[end,2:end-1]
                    ind = [1;2;3;5];
                elseif i in t[2:end-1,1]
                    ind = [2;3;4;5];
                elseif i in t[2:end-1,end]
                    ind = [1;2;3;4]
                elseif i == t[1,1]
                    ind = [3;4;5];
                elseif i == t[1,end]
                    ind = [1;3;4];
                elseif i == t[end,1]
                    ind = [2;3;5];
                elseif i == t[end,end]
                    ind = [1;2;3];
                end
            elseif patch == "RB"
                if i in indices
                    ind = [1;2;3;4;5];
                elseif i in t[1,2:end-1]
                    ind = [2;3;5];
                elseif i in t[end,2:end-1]
                    ind = [1;3;4];
                elseif i in t[2:end-1,1]
                    ind = [3;4;5];
                elseif i in t[2:end-1,end]
                    ind = [1;2;3];
                elseif i == t[1,1]
                    ind = [3;5];
                elseif i == t[1,end]
                    ind = [2;3];
                elseif i == t[end,1]
                    ind = [3;4];
                elseif i == t[end,end]
                    ind = [1;3];
                end
            elseif patch == "Full"
                if i in indices
                    ind = [1;2;3;4;5;6;7;8;9];
                elseif i in t[1,2:end-1]
                    ind = [2;3;5;6;8;9];
                elseif i in t[end,2:end-1]
                    ind = [1;2;4;5;7;8];
                elseif i in t[2:end-1,1]
                    ind = [4;5;6;7;8;9];
                elseif i in t[2:end-1,end]
                    ind = [1;2;3;4;5;6];
                elseif i == t[1,1]
                    ind = [5;6;8;9];
                elseif i == t[1,end]
                    ind = [2;3;5;6];
                elseif i == t[end,1]
                    ind = [4;5;7;8];
                elseif i == t[end,end]
                    ind = [1;2;4;5];
                end
            elseif patch == "RB_wide"
                if i in indices
                    ind = [1;2;3;4;5;6;7;8;9];
                elseif i in t[3:end-2,2]
                    ind = [2;3;4;5;6;7;8;9];
                elseif i in t[3:end-2,end-1]
                    ind = [1;2;3;4;5;6;7;8];
                elseif i in t[2,3:end-2]
                    ind = [1;2;3;5;6;7;8;9];
                elseif i in t[end-1,3:end-2]
                    ind = [1;2;3;4;5;7;8;9]
                elseif i == t[2,2]
                    ind = [2;3;5;6;7;8;9];
                elseif i == t[2,end-1]
                    ind = [1;2;3;5;6;7;8];
                elseif i == t[end-1,2]
                    ind = [2;3;4;5;7;8;9];
                elseif i == t[end-1,end-1]
                    ind = [1;2;3;4;5;7;8];
                else
                    ind = [];
                end
            end

        elseif length(n) == 3 # 3D
            if patch == "Plus" 
                if i in indices
                    ind = [1;2;3;4;5;6;7];
                ### faces ###
                elseif i in t[2:end-1,1,2:end-1]; # top face
                    ind = [1;3;4;5;6;7];
                elseif i in t[2:end-1,end,2:end-1]; # bottom face
                    ind = [1;2;3;4;5;7];
                elseif i in t[2:end-1,2:end-1,end]; # front face
                    ind = [1;2;3;4;5;6];
                elseif i in t[2:end-1,2:end-1,1]; # rear face
                    ind = [2;3;4;5;6;7];
                elseif i in t[1,2:end-1,2:end-1]; # left face
                    ind = [1;2;4;5;6;7];
                elseif i in t[end,2:end-1,2:end-1]; # right face
                    ind = [1;2;3;4;6;7];
                ### edges ###
                elseif i in t[2:end-1,1,end]; # top-front edge
                    ind = [1;3;4;5;6];
                elseif i in t[2:end-1,1,1]; # top-rear edge
                    ind = [3;4;5;6;7];
                elseif i in t[2:end-1,end,end]; # bottom-front edge
                    ind = [1;2;3;4;5];
                elseif i in t[2:end-1,end,1]; # bottom-rear edge
                    ind = [2;3;4;5;7];
                elseif i in t[1,1,2:end-1]; # top-left edge
                    ind = [1;4;5;6;7];
                elseif i in t[end,1,2:end-1]; # top-right edge
                    ind = [1;3;4;6;7];
                elseif i in t[1,end,2:end-1]; # bottom-left edge
                    ind = [1;2;4;5;7];
                elseif i in t[end,end,2:end-1]; # bottom-right edge
                    ind = [1;2;3;4;7];
                elseif i in t[1,2:end-1,end]; # left-front edge
                    ind = [1;2;4;5;6];
                elseif i in t[1,2:end-1,1]; # left-rear edge
                    ind = [2;4;5;6;7];
                elseif i in t[end,2:end-1,end]; # right-front edge
                    ind = [1;2;3;4;6];
                elseif i in t[end,2:end-1,1]; # right-rear edge
                    ind = [2;3;4;6;7];  
                ### corners ###
                elseif i == t[1,1,1]; # left-top-rear corner
                    ind = [4;5;6;7];
                elseif i == t[end,1,1]; # right-top-rear corner
                    ind = [3;4;6;7];
                elseif i == t[1,end,1]; # left-bottom-rear corner
                    ind = [2;4;5;7];
                elseif i == t[1,1,end]; # left-top-front corner
                    ind = [1;4;5;6];
                elseif i == t[end,end,1]; # right-bottom-rear corner
                    ind = [2;3;4;7];
                elseif i == t[end,1,end]; # right-top-front corner
                    ind = [1;3;4;6];
                elseif i == t[1,end,end]; # left-bottom-front corner
                    ind = [1;2;4;5];
                elseif i == t[end,end,end]; # right-bottom-front corner
                    ind = [1;2;3;4];
                end

            elseif patch == "RB"
                if i in indices
                ind = [1;2;3;4;5;6;7;8;9;10;11;12;13];
                ### faces ###
                elseif i in t[2:end-1,1,2:end-1]; # top face
                    ind = [2;3;4;7;8;9;11;12;13];
                elseif i in t[2:end-1,end,2:end-1]; # bottom face
                    ind = [1;2;3;5;6;7;10;11;12];
                elseif i in t[2:end-1,2:end-1,end]; # front face
                    ind = [1;2;3;4;5;6;7;8;9];
                elseif i in t[2:end-1,2:end-1,1]; # rear face
                    ind = [5;6;7;8;9;10;11;12;13];
                elseif i in t[1,2:end-1,2:end-1]; # left face
                    ind = [1;3;4;6;7;9;10;12;13];
                elseif i in t[end,2:end-1,2:end-1]; # right face
                    ind = [1;2;4;5;7;8;10;11;13];
                ### edges ###
                elseif i in t[2:end-1,1,end]; # top-front edge
                    ind = [2;3;4;7;8;9];
                elseif i in t[2:end-1,1,1]; # top-rear edge
                    ind = [7;8;9;11;12;13];
                elseif i in t[2:end-1,end,end]; # bottom-front edge
                    ind = [1;2;3;5;6;7];
                elseif i in t[2:end-1,end,1]; # bottom-rear edge
                    ind = [5;6;7;10;11;12];
                elseif i in t[1,1,2:end-1]; # top-left edge
                    ind = [3;4;7;9;12;13];
                elseif i in t[end,1,2:end-1]; # top-right edge
                    ind = [2;4;7;8;11;13];
                elseif i in t[1,end,2:end-1]; # bottom-left edge
                    ind = [1;3;6;7;10;12];
                elseif i in t[end,end,2:end-1]; # bottom-right edge
                    ind = [1;2;5;7;10;11];
                elseif i in t[1,2:end-1,end]; # left-front edge
                    ind = [1;3;4;6;7;9];
                elseif i in t[1,2:end-1,1]; # left-rear edge
                    ind = [6;7;9;10;12;13];
                elseif i in t[end,2:end-1,end]; # right-front edge
                    ind = [1;2;4;5;7;8];
                elseif i in t[end,2:end-1,1]; # right-rear edge
                    ind = [5;7;8;10;11;13];   
                ### corners ###
                elseif i == t[1,1,1]; # left-top-rear corner
                    ind = [7;9;12;13];
                elseif i == t[end,1,1]; # right-top-rear corner
                    ind = [7;8;11;13];
                elseif i == t[1,end,1]; # left-bottom-rear corner
                    ind = [6;7;10;12];
                elseif i == t[1,1,end]; # left-top-front corner
                    ind = [3;4;7;9];
                elseif i == t[end,end,1]; # right-bottom-rear corner
                    ind = [5;7;10;11];
                elseif i == t[end,1,end]; # right-top-front corner
                    ind = [2;4;7;8];
                elseif i == t[1,end,end]; # left-bottom-front corner
                    ind = [1;3;6;7];
                elseif i == t[end,end,end]; # right-bottom-front corner
                    ind = [1;2;5;7];
                end
            end
        end
        Idxs = block_id[ind] .+ i;

        # creating the inverses of local Blocks
        k = length(Idxs);
    
        Mi = @inbounds [A[Idxs[ii], Idxs[jj]] for ii=1:k, jj=1:k]


        ##########################################################

        W = (1/k) * ones(size(Mi,1)); # equal weights
        for j = 1:k
            W[j] = W_grid[Idxs[j]]; # weights from preprepared array
        end

        Mi = W .* inv(Mi)

        # Store contributions in a large matrix as sparse CSC
        for ii = 1:k
            for jj = 1:k
                rows[p] = Idxs[ii]
                cols[p] = Idxs[jj]
                vals[p] = Mi[ii,jj]
                p += 1
            end
        end
    end

    rows[p:end] .= 1
    cols[p:end] .= 1
    vals[p:end] .= 0
    M_all = sparse(rows, cols, vals, size(A,1), size(A,2))

    return M_all

end


function Vanka(A,b,w,x,M_all,tol,maxit)

    wMall = w .* M_all;
    r = b - A * x;
    iter = 0;

    for i=1:maxit
        if (norm(r) > tol)
            x += wMall * r;
            r = b - A * x;
            iter = iter + 1;
        end
    end

    return x
end


############## test Vanka ################

# n = [4;4;4];
# nodes = n .+ 1;

# A = spdiagm(ones(prod(nodes)));
# b = ones(prod(nodes));
# x = 0.0*b;
# patch = "RB";

# M_all = VankaSetup(A,n,patch);
# tol = 1e-5;
# maxit = 1;
# w = 1;

# x = Vanka(A,b,w,x,M_all,tol,maxit)

# x = real(reshape(x,nodes...))