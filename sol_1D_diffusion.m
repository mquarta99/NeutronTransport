clear all
close all
clc

% Problem parameters
H = 100;            
Nz = 200;           
dz = H / Nz;        
z = linspace(dz/2, H-dz/2, Nz); 

sigma_t = 0.66962;     
sigma_sc = 0.64117;             

% External source
S = zeros(Nz, 1);
S(Nz/4+1:3*Nz/4) = 1; % Localized source in the center

Phi = solve_diffusion_1D(Nz,dz,1/3/sigma_t,sigma_t-sigma_sc,S);

figure;
yyaxis left
plot(z, Phi, 'b', 'LineWidth', 2);
ylabel('Neutron Flux \phi(z)');
xlabel('Position z');
grid on;
title('Scalar Flux and Source Distribution in 1D Slab');
yyaxis right
plot(z, S, 'r--', 'LineWidth', 2);
ylabel('External Source S(z)');
ylim([0, 1.2])
legend('Neutron Flux','External Source');