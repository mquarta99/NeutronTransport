clear all
close all

%develop solver for 1D slab SN with isotropic scattering
%uniform mesh
%zero incoming flux BC
%use cross section without fission:
%%%%%CASE 1: USE 1-GROUP CROSS SECTIONS --> SOURCE ITERATION ONLY
%%%%%CASE 2: USE G-GROUP CROSS SECTIONS --> ADD THERMAL ITERATIONS

%source localized somewhere
H=100;

Nz=200;
dz=H/Nz;
z=[dz/2:dz:H-dz/2];
N=4; %sn order

%generation of S_N directions and weights
[mu,weig]=lgwt(N,-1.0,1.0);

%definition of source in the central part of the system
% SS=zeros(Nz,1);
% SS(Nz/4+1:3*Nz/4)=1;

% Nuclear properties (1-group)
sigma_t = 0.669618339332488;     
sigma_sc = 0.641172644495600;     
nuSigma_f = 0.035033062356715;
%chi = 1.0;          


%% POWER METHOD 

% Initial guess
k = 1.0;
Phi_pow = ones(Nz,1); 
%Phi_pow = sin(pi * z / H); % Better guess

F = nuSigma_f; % Fission operator

toll = 1e-5;
err_pow = 1;
it_pow = 0;

while err_pow > toll
    
    it_pow = it_pow + 1;
    
    S = 1/k * F * Phi_pow; % New source term
    Phi = Phi_pow; % I use this as initial guess for the source iteration
    Phi_old = Phi;
    err_src = 1;
    it_src = 0;
    while err_src > toll
        it_src = it_src + 1;
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
        Phi = sum(weig' .* phi, 2); 
        
        % Compute relative error
        err_src = max(abs(Phi - Phi_old)) / max(abs(Phi));
        Phi_old = Phi;
    end

    k_new = k * sum(Phi) / sum(Phi_pow);
    err_pow = abs(k_new - k) / k_new;
    k = k_new;
    Phi_pow = Phi;
    fprintf('Iteration %d - Error: %.5e\n', it_pow, err_pow);
end

disp(['k_effective: ', num2str(k)]);

figure;
plot(z, Phi_pow, 'b', 'LineWidth', 2);
xlabel('Position z');
ylabel('Neutron Flux \Phi(z)');
title('Scalar Flux Distribution in 1D Slab');
grid on;