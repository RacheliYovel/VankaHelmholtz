function [rho_loc,mu_loc,w] = optimizeDamping(h,gppw,gamma,smoother,samples,stencil,nu,CGA,intergrid)

    wre = 0.8:0.01:1.1;
    wim = -0.01:0.01:0.01;
    rho_loc_arr = zeros(length(wre),length(wim));
    mu_loc_arr = zeros(length(wre),length(wim));
    for i = 1:length(wre)
        for j=1:length(wim)
            w = wre(i) + wim(j) * 1i;
            [rho_loc_arr(i,j),mu_loc_arr(i,j)] = AcousticHelmholtzLFA(h,gppw,gamma,w,smoother,samples,stencil,nu,CGA,intergrid);
        end
    end

    % disp(rho_loc_arr)
    
    % figure()
    % contour(wre,wim,rho_loc_arr')
    % colorbar
    % title(intergrid," intergrid, ",smoother," smoother")
    % xlabel("Real part")
    % ylabel("Imaginary part")
    
    % [rho_loc, n] = min(rho_loc_arr(:));
    % [x, y] = ind2sub(size(rho_loc_arr),n);
    % mu_loc = mu_loc_arr(x,y);
    % w = wre(x) + (wim(y))*1i;


    [mu_loc, n] = min(mu_loc_arr(:));
    [x, y] = ind2sub(size(mu_loc_arr),n);
    rho_loc = rho_loc_arr(x,y);
    w = wre(x) + (wim(y))*1i;

end