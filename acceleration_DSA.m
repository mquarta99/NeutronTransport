clear all
close all
clc

%% Source iteration method

% Problem parameters
H = 100;            
Nz = 100;           
N = 4;    % Order of S_N method
N_low = 2; % Order lower order Sn solver
dz = H/Nz;
z = dz/2 : dz : H-dz/2;

% Nuclear properties (1-group)
sigma_t =  0.66962;     
sigma_sc = 0.64117;     

% External source
S = zeros(Nz, 1);
S(Nz/4+1:3*Nz/4) = 1; % Localized source in the center

tic;
[Phi, it1,err1] = solve_Sn_oneG(H,Nz,N,sigma_t,sigma_sc,S);
toc;

tic;
[Phi_acc, it2, err2] = solve_Sn_oneG_accelerated(H,Nz,N,N_low,sigma_t,sigma_sc,S,"transport");
toc;

% Plot scalar flux distribution and external source
figure;
yyaxis left
plot(z, Phi, 'b', 'LineWidth', 2);
hold on;
plot(z,Phi_acc,'g','LineWidth',2);
ylabel('Neutron Flux \phi(z)');
xlabel('Position z');
title('Scalar Flux and Source Distribution in 1D Slab');
grid on;
yyaxis right
plot(z, S, 'r--', 'LineWidth', 2);
ylabel('External Source S(z)');
ylim([0, 1.2])
legend('Neutron Flux','Neutron Flux - acceleration','External Source');

figure;
semilogy(it1, err1, 'LineWidth', 2); hold on;
semilogy(it2, err2, 'LineWidth', 2); 
xlabel('Iteration');
ylabel('Relative error on Phi');
title('Convergence of Phi');
legend('without acceleration','with acceleration')
grid on;
