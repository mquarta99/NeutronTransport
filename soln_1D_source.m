clear all
close all
clc

%develop solver for 1D slab SN with isotropic scattering
%uniform mesh
%zero incoming flux BC
%use cross section without fission:
%%%%%CASE 1: USE 1-GROUP CROSS SECTIONS --> SOURCE ITERATION ONLY
%%%%%CASE 2: USE G-GROUP CROSS SECTIONS --> ADD THERMAL ITERATIONS


% Problem parameters
H = 100;            
Nz = 200;           
dz = H / Nz;        
z = linspace(dz/2, H-dz/2, Nz); 

N = 4;              % Order of S_N method
[mu, weig] = lgwt(N, -1, 1); 

% Nuclear properties (1-group)
sigma_t = 0.669618339332488;     
sigma_sc = 0.641172644495600;     
% nu = 2.504747808339256;         
% sigma_f = 0.014233685639664;   
% chi = 1.0;                    

% External source
S = zeros(Nz, 1);
S(Nz/4+1:3*Nz/4) = 1; % Localized source in the center

% Initial guess of the scalar flux
Phi = ones(Nz, 1); % Constant flux
%Phi = zeros(Nz,1); Phi(Nz/4:3*Nz/4)=35; % Better guess
Phi_old = Phi;




% Source iteration parameters
toll = 1e-5;
err = 1;
it = 0;

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

    % Compute relative error
    err = max(abs(Phi - Phi_old)) / max(abs(Phi));
    Phi_old = Phi;

    fprintf('Iteration %d - Error: %.5e\n', it, err);
end

% Plot scalar flux distribution and external source
figure;
yyaxis left
plot(z, Phi, 'b', 'LineWidth', 2);
ylabel('Neutron Flux \phi(z)');
xlabel('Position z');
title('Scalar Flux and Source Distribution in 1D Slab');
grid on;

yyaxis right
plot(z, S, 'r--', 'LineWidth', 2);
ylabel('External Source S(z)');
ylim([0, 1.2])
legend('\phi(z) - Neutron Flux', 'S(z) - External Source');
