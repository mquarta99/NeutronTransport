clear all
close all
clc

%develop solver for 1D slab SN with isotropic scattering
%uniform mesh
%zero incoming flux BC
%use cross section without fission:
%%%%%CASE 1: USE 1-GROUP CROSS SECTIONS --> SOURCE ITERATION ONLY
%%%%%CASE 2: USE G-GROUP CROSS SECTIONS --> ADD THERMAL ITERATIONS

%% 4G case

% Problem parameters
H = 100;            
Nz = 200;           
dz = H / Nz;        
z = linspace(dz/2, H-dz/2, Nz); 

N = 4;              % Order of S_N method
[mu, weig] = lgwt(N, -1, 1); 

% Nuclear properties (1-group)
ENRG = [2.00000000e+01, 8.21000000e-01, 9.11800000e-03, 1.02000000e-06];
N_groups = 4;
sigma_t = [0.272215751922190, 0.596492773574977, 0.938697442352708, 1.30784498748829];   

sigma_sc = [0.188207704213408,     0,                     0,                   0;
            0.0767803256491297,    0.525825002305149,     0,                   0;
            0.000515991069970616,  0.0669072168907535,    0.851079520721167,   0.000749233230385878;
            4.98756447025729e-08,  6.86364083003606e-06,  0.0438703502594663,  1.14783206728197];

% nuSigma_f = [0.0150412957883228, 0.00201452404189487, 0.0257648566047470, 0.289323130396279];  
% chi = [0.751496026100000, 0.248281023300000, 0.000222950628000000, 0];                    

% External source
SS = zeros(Nz, 1);
SS(Nz/4+1:3*Nz/4) = 1; % Localized source in the center

% hyp: All emitted in first group, with energy > 0.821 MeV 

% Initial guess
Phi = ones(Nz,N_groups);

%% Thermal iterations
for jj = 1:N_groups 
    
    up_scattering = false;
    for kk = jj+1:N_groups
        if sigma_sc(jj,kk) ~= 0
            up_scattering = true;
            break
        else
        end
    end
    fprintf("Group %d - Upscattering %s\n", jj, string(up_scattering));
    S = zeros(Nz,1); % Initialize the source
    if up_scattering == false
        % no upscattering -> no thermal iterations
        % directly source iterations
        if jj == 1
            S = SS; %*ENRG(1)/(ENRG(1)-ENRG(2)); % First group: no downscattering from higher energies, but external source
        else
            for tt = 1:jj-1 % Other groups: downscattering from higher energies
                S = S + sigma_sc(jj,tt)*Phi(:,tt);
                fprintf("Downscattering from group %d\n",tt);
            end
        end
        toll_src = 1e-5;
        err_src = 1;
        it_src = 0;
        Phi_j = Phi(:,jj); % Initial guess
        Phi_j_old = Phi_j;
        while err_src > toll_src
            it_src = it_src + 1;
            phi = zeros(Nz, N); 
            % Here I discretize using Sn and FV + Diamond difference
            for n = 1:N
                if mu(n) > 0 % Transport Sweep from left to right (mu > 0)
                    phi_inlet = 0; % Zero incoming flux at the left boundary
                    for i = 1:Nz
                        % Diamond Difference Scheme
                        Q_tilde = S(i)/2 + sigma_sc(jj,jj)/2*Phi_j(i);
                        phi_in = (2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t(jj) * dz + 2 * mu(n));
                        phi(i, n) = phi_in;
                        phi_inlet = 2 * phi_in - phi_inlet; % Update for the next step
                    end
                end
        
                if mu(n) < 0
                    phi_inlet = 0; % Zero incoming flux at the right boundary
                    for i = Nz:-1:1
                        % Diamond Difference Scheme
                        Q_tilde = S(i)/2 + sigma_sc(jj,jj)/2*Phi_j(i);
                        phi_in = (- 2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t(jj) * dz - 2 * mu(n));
                        phi(i, n) = phi_in;
                        phi_inlet = 2 * phi_in - phi_inlet; % Update for the next step
                    end
                end
            end
        
            % Update scalar flux (weighted average)
            Phi_j = zeros(Nz, 1);
            for i = 1:Nz
                for n = 1:N
                    Phi_j(i) = Phi_j(i) + weig(n) * phi(i, n);
                end
            end 
        
            % Compute relative error
            err_src = max(abs(Phi_j - Phi_j_old)) / max(abs(Phi_j));
            Phi_j_old = Phi_j;
        
            fprintf('Source Iteration %d - Error: %.5e\n', it_src, err_src);
        end
    fprintf("\nSource Iterations for group %d - Convergence reached \n",jj);
    Phi(:,jj) = Phi_j; % Source iterations finished, I stock the result
    else 
        % upscatterig true -> thermal iterations needed
        toll_thermal = 1e-5;
        err_thermal = 1;
        it_thermal = 0;
        thermal_groups = jj:N_groups; % Groups which are included in the thermal
        Phi_thermal = Phi; % Converged groups + Initial guesses for the thermal groups
        
        while err_thermal > toll_thermal
            it_thermal = it_thermal + 1;
            
            for gg = thermal_groups(1):thermal_groups(end)
                S = zeros(Nz,1);
                for tt = 1:gg-1 % Downscattering from higher energies
                    S = S + sigma_sc(gg,tt)*Phi_thermal(:,tt);
                    fprintf("Downscattering from group %d\n",tt);
                end
                for tt = gg+1:N_groups % Upscattering from lower groups
                    S = S + sigma_sc(gg,tt)*Phi_thermal(:,tt); % This is a guess
                    fprintf("Upscattering from group %d\n",tt);
                end
                % I can use now the source iterations
                toll_src = 1e-5;
                err_src = 1;
                it_src = 0;
                Phi_g = Phi_thermal(:,gg); % Initial guess
                Phi_g_old = Phi_g;
                while err_src > toll_src
                    it_src = it_src + 1;
                    phi = zeros(Nz, N); 
                    % Here I discretize using Sn and FV + Diamond difference
                    for n = 1:N
                        if mu(n) > 0 % Transport Sweep from left to right (mu > 0)
                            phi_inlet = 0; % Zero incoming flux at the left boundary
                            for i = 1:Nz
                                % Diamond Difference Scheme
                                Q_tilde = S(i)/2 + sigma_sc(gg,gg)/2*Phi_g(i);
                                phi_in = (2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t(gg) * dz + 2 * mu(n));
                                phi(i, n) = phi_in;
                                phi_inlet = 2 * phi_in - phi_inlet; % Update for the next step
                            end
                        end
                
                        if mu(n) < 0
                            phi_inlet = 0; % Zero incoming flux at the right boundary
                            for i = Nz:-1:1
                                % Diamond Difference Scheme
                                Q_tilde = S(i)/2 + sigma_sc(gg,gg)/2*Phi_g(i);
                                phi_in = (- 2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t(gg) * dz - 2 * mu(n));
                                phi(i, n) = phi_in;
                                phi_inlet = 2 * phi_in - phi_inlet; % Update for the next step
                            end
                        end
                    end
                
                    % Update scalar flux (weighted average)
                    Phi_g = zeros(Nz, 1);
                    for i = 1:Nz
                        for n = 1:N
                            Phi_g(i) = Phi_g(i) + weig(n) * phi(i, n);
                        end
                    end 
                
                    % Compute relative error
                    err_src = max(abs(Phi_g - Phi_g_old)) / max(abs(Phi_g));
                    Phi_g_old = Phi_g;
                
                    fprintf('Source Iterations %d - Error: %.5e\n', it_src, err_src);
                end
                fprintf("\nSource Iterations for group %d - Convergence reached \n",gg);
                Phi_thermal(:,gg) = Phi_g; % Source iterations finished, I stock the result
            end
            err_thermal = max(max(abs(Phi_thermal - Phi)) ./ max(abs(Phi)));
            Phi = Phi_thermal;
            fprintf('Thermal Iterations %d - Error: %.5e\n', it_thermal, err_thermal);
        end
        fprintf("\nThermal Iterations - Convergence reached \n");
        break
    end
end

Phi_tot = sum(Phi,2);

% Plot scalar flux distribution and external source
figure;
hold on;
for jj = 1:N_groups
    plot(z, Phi(:, jj), 'LineWidth', 2, 'DisplayName', sprintf('\\phi_{%d}(z)', jj));
end
plot(z, sum(Phi, 2), 'k--', 'LineWidth', 2, 'DisplayName', '\phi_{tot}(z)');
ylabel('Neutron Flux \phi(z)');
xlabel('Position z');
title('Scalar Flux Distribution in 1D Slab');
grid on;
yyaxis right
plot(z, SS, 'r--', 'LineWidth', 2, 'DisplayName', 'S(z) - External Source');
ylabel('External Source S(z)');
ylim([0, 1.2]);
legend();
hold off;
