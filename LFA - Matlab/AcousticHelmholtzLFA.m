function [rho_loc,mu_loc] = AcousticHelmholtzLFA(h,gppw,gamma,w,smoother,samples,stencil,nu,CGA,intergrid)

[V,W,invHi,sigma] = AdditiveVankaSetup(h,gppw,gamma,smoother,stencil);
h2inv = 1/h^2;
omega = (2 * pi)/(h * gppw);
Sigma = omega^2 * (1 - 2* gamma * 1i);

theta1 = linspace(-pi/2, 3*pi/2, samples);
theta2 = linspace(-pi/2, 3*pi/2, samples);

spec_rad_arr = zeros(length(theta1),length(theta2));
smoother_eig_arr = zeros(length(theta1),length(theta2));

for i = 1:length(theta1)
    for j = 1:length(theta2)

        % space of harmonics
        theta = [theta1(i) theta2(j)];
        theta_arr = [theta; (theta+[pi pi]); (theta+[pi 0]); (theta+[0 pi])];

        % fine stencil
        if (stencil == 5)
            s1 = 0; s2 = -1*h2inv; s3 = 4*h2inv-sigma;
        elseif (stencil == 9)
            s1 = -1*h2inv; s2 = -1*h2inv - (1/6)*sigma; s3 = 8*h2inv - (1/3)*sigma;  
        end

        % symbol matrices
        S = zeros(4);
        R_bilin = zeros(1,4);
        R_bicub = zeros(1,4);
        R_hanshaw = zeros(1,4);  
        Af = zeros(4);
        for l=1:4
            t = theta_arr(l,:);

            S(l,l) = AdditiveVanka(V,W,invHi,h,sigma,t,w,smoother,stencil);

            R_bilin(l) = (1/4)*(1+cos(t(1)))*(1+cos(t(2))); 
            R_bicub(l) = symbol5by5(6/256,4/256,1/256,t);
            R_hanshaw(l) = (1/256)*(-8-9*cos(t(1))+cos(3*t(1)))*(-8-9*cos(t(2))+cos(3*t(2)));

            % symbol matrix of the fine grid operator
            Af(l,l) = s3 + 2*s2*cos(t(1)) + 2*s2*cos(t(2)) + 4*s1*cos(t(1))*cos(t(2));
        end

        smoother_eig_arr(i,j) = S(1,1);

        if intergrid == "FW"
            R = R_bilin;
            P = R';
        elseif intergrid == "cubic"
            R = R_bicub;
            P = R';
        elseif intergrid == "HO"
            R = R_hanshaw;
            P = R';
        elseif intergrid == "mixed"
            R = R_bilin;
            P = R_bicub';
        end
       

        % Coarse grid operator
        if (CGA == "GCA")
            Ac = R*Af*P;
        elseif (CGA == "DCA")
            Theta = 2 .* theta;
            H = 2 .* h;
            if (stencil == 5)
                S1 = 0; S2 = -1/(H^2); S3 = 4/(H^2)-Sigma;
            elseif (stencil == 9)
                S1 = -(1/6)*(1/H^2); S2 = -(2/3)*(1/H^2) - (1/12)*Sigma; S3 = (10/3)*(1/H^2) - (2/3)*Sigma;
            end
            Ac = S3 + 2*S2*cos(Theta(1)) + 2*S2*cos(Theta(2)) + 4*S1*cos(Theta(1))*cos(Theta(2));
        elseif (CGA == "GCAlikeDCA")
            Theta = 2 .* theta;
            H = 2 .* h;
            s0 = (1/H^2)*(1/256)*(-98) + (Sigma/64^2)*70;
            s1 = (1/H^2)*(1/256)*(-44) + (Sigma/64^2)*28;
            s2 = (1/H^2)*(1/256)*(-3) + (Sigma/64^2)*1;
            s0 = s0/2; s1 = s1/2; s2 = s2/2;
            Ac = symbol5by5(s0,s1,s2,Theta);
        end
        Ac_inv = 1 / Ac;

        % symbol matrix of the two-grid
        K = eye(4) - P * Ac_inv * R * Af;
        M =  S^nu * K ; 
        spec = eig(M);
        spec_rad_arr(i,j) = (abs(max(spec)));     

    end
end

ind_pihalf = floor(samples/2);
rho_loc = max(max(spec_rad_arr(1:ind_pihalf,1:ind_pihalf)));

mu1 = max(max(abs(smoother_eig_arr(1:ind_pihalf,ind_pihalf+1:end))));
mu2 = max(max(abs(smoother_eig_arr(ind_pihalf+1:end,1:ind_pihalf))));
mu3 = max(max(abs(smoother_eig_arr(ind_pihalf+1:end,ind_pihalf+1:end))));

mu_loc = max([mu1 mu2 mu3]);

end