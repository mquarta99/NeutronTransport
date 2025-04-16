function [Phi,it_vec,err_vec] = solve_Sn_oneG(H,Nz,N,sigma_t,sigma_sc,S)
    
    dz = H / Nz;
    z = linspace(dz/2, H-dz/2, Nz); 
    [mu, weig] = lgwt(N,-1, 1);
    Phi = ones(Nz, 1); % Constant flux
    Phi_old = Phi;
    toll = 1e-5;
    err = 1;
    it = 0;
    it_vec = [];
    err_vec = [];
    while err > toll
        it = it + 1;
        phi = zeros(Nz, N); 
        
        for n = 1:N
            if mu(n) > 0 % Transport Sweep from left to right (mu > 0)
                phi_inlet = 0; % Zero incoming flux at the left boundary
                for i = 1:Nz
                    % Diamond Difference Scheme
                    Q_tilde = S(i)/2 + sigma_sc/2*Phi(i);
                    phi_in = (2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t * dz + 2 * mu(n));
                    phi(i, n) = phi_in;
                    phi_inlet = 2 * phi_in - phi_inlet; % Update for the next step
                end
            end
    
            if mu(n) < 0
                phi_inlet = 0; % Zero incoming flux at the right boundary
                for i = Nz:-1:1
                    % Diamond Difference Scheme
                    Q_tilde = S(i)/2 + sigma_sc/2*Phi(i);
                    phi_in = (- 2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t * dz - 2 * mu(n));
                    phi(i, n) = phi_in;
                    phi_inlet = 2 * phi_in - phi_inlet; % Update for the next step
                end
            end
        end
    
        % Update scalar flux (sum over all angular directions)
        Phi = zeros(Nz, 1);
        for i = 1:Nz
            for n = 1:N
                Phi(i) = Phi(i) + weig(n) * phi(i, n);
            end
        end 
        
        diff_Phi = Phi - Phi_old;
        % Compute relative error
        err = sqrt(sum(diff_Phi.^2)) / sqrt(sum(Phi.^2));
        Phi_old = Phi;
        it_vec(end+1) = it;
        err_vec(end+1) = err;
        %fprintf('Iteration %d - Error: %.5e\n', it, err);
    end
end
