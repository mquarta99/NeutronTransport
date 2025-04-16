clear all
close all

%% Acceleration method for the power method

H = 100;
Nz = 200;
dz = H / Nz;
z = dz/2 : dz : H - dz/2;
N = 4; % S_N order

% Generazione delle direzioni e dei pesi S_N
[mu, weig] = lgwt(N, -1.0, 1.0);

% Proprietà nucleari (1 gruppo)
sigma_t = 0.66962;
sigma_sc = 0.64117;
nuSigma_f = 0.03503;

%% Without acceleration

% Initial guess
k = 1.0;
Phi_pow = ones(Nz, 1); 

F = nuSigma_f; % Operatore di fissione

toll = 1e-5; % 1 pcm
err_k = 1;

% Inizializzazione errori
err_k_vec = [];
err_phi_L2 = [];
err_phi_Linf = [];

it_pow = 0;

while err_k > toll
    it_pow = it_pow + 1;

    S = 1/k * F * Phi_pow; % Sorgente nuova
    Phi = Phi_pow; 
    Phi_old = Phi;
    err_src = 1;
    it_src = 0;

    while err_src > toll
        it_src = it_src + 1;
        phi = zeros(Nz, N); 

        for n = 1:N
            if mu(n) > 0
                phi_inlet = 0; 
                for i = 1:Nz
                    Q_tilde = S(i)/2 + sigma_sc/2 * Phi(i);
                    phi_in = (2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t * dz + 2 * mu(n));
                    phi(i, n) = phi_in;
                    phi_inlet = 2 * phi_in - phi_inlet;
                end
            end

            if mu(n) < 0
                phi_inlet = 0;
                for i = Nz:-1:1
                    Q_tilde = S(i)/2 + sigma_sc/2 * Phi(i);
                    phi_in = (-2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t * dz - 2 * mu(n));
                    phi(i, n) = phi_in;
                    phi_inlet = 2 * phi_in - phi_inlet;
                end
            end
        end

        % Aggiorna flusso scalare
        Phi = sum(weig' .* phi, 2);

        % Errore sorgente
        err_src = max(abs(Phi - Phi_old)) / max(abs(Phi));
        Phi_old = Phi;
    end

    k_new = k * sum(F * Phi) / sum(F * Phi_pow);

    err_k = abs(k_new - k) / k_new;
    diff_phi = Phi - Phi_pow;
    
    err_phi_L2(end+1) = sqrt(sum (diff_phi.^2))/ sqrt(sum(Phi.^2));
    err_phi_Linf(end+1) = max(abs(diff_phi)) / max(abs(Phi));
    err_k_vec(end+1) = err_k;

    k = k_new;
    Phi_pow = Phi;

    %fprintf('Iterazione %d - Errore k: %.5e\n', it_pow, err_k);
end

disp(['Without acceleration - k_eff =  ', num2str(k) ]);


% figure;
% plot(z, Phi_pow, 'b', 'LineWidth', 2);
% xlabel('z');
% ylabel('Neutron Flux \Phi(z)');
% title('Neutron Flux Distribution in a 1D Slab');
% grid on;


% figure;
% semilogy(1:it_pow, err_k_vec, 'LineWidth', 2); hold on;
% semilogy(1:it_pow, err_phi_L2, 'LineWidth', 2); 
% semilogy(1:it_pow, err_phi_Linf, 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Relative error');
% title('Convergence of k  and Phi');
% legend('k_{eff}','Phi L_{2}','Phi L_{\infty}')
% grid on;


%% With acceleration

% Initial guess
k = 1.0;
Phi_pow = ones(Nz, 1); 

F = nuSigma_f; % Operatore di fissione

toll = 1e-5; % 1 pcm
err_k = 1;

% Inizializzazione errori
err_k_vec_acc = [];
err_phi_L2_acc = [];
err_phi_Linf_acc = [];

it_pow_acc = 0;

while err_k > toll
    it_pow_acc = it_pow_acc + 1;

    S = 1/k * F * Phi_pow; % Sorgente nuova
    Phi = Phi_pow; 
    Phi_old = Phi;
    err_src = 1;
    it_src = 0;

    while err_src > toll
        it_src = it_src + 1;
        phi = zeros(Nz, N); 

        for n = 1:N
            if mu(n) > 0
                phi_inlet = 0; 
                for i = 1:Nz
                    Q_tilde = S(i)/2 + sigma_sc/2 * Phi(i);
                    phi_in = (2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t * dz + 2 * mu(n));
                    phi(i, n) = phi_in;
                    phi_inlet = 2 * phi_in - phi_inlet;
                end
            end

            if mu(n) < 0
                phi_inlet = 0;
                for i = Nz:-1:1
                    Q_tilde = S(i)/2 + sigma_sc/2 * Phi(i);
                    phi_in = (-2 * mu(n) * phi_inlet + Q_tilde * dz) / (sigma_t * dz - 2 * mu(n));
                    phi(i, n) = phi_in;
                    phi_inlet = 2 * phi_in - phi_inlet;
                end
            end
        end

        % Aggiorna flusso scalare
        Phi = sum(weig' .* phi, 2);

        % Errore sorgente
        err_src = max(abs(Phi - Phi_old)) / max(abs(Phi));
        Phi_old = Phi;
    end

    k_new = k * sum(Phi .* (F * Phi)) / sum(Phi .* (F * Phi_pow));

    err_k = abs(k_new - k) / k_new;
    diff_phi = Phi - Phi_pow;
    
    err_phi_L2_acc(end+1) = sqrt(sum(diff_phi.^2)) / sqrt(sum(Phi.^2));
    err_phi_Linf_acc(end+1) = max(abs(diff_phi)) / max(abs(Phi));
    err_k_vec_acc(end+1) = err_k;

    k = k_new;
    Phi_pow = Phi;

    %fprintf('Iterazione %d - Errore k: %.5e\n', it_pow_acc, err_k);
end

disp(['Acdelerated - k_eff =  ', num2str(k) ]);

figure;
plot(z, Phi_pow, 'b', 'LineWidth', 2);
xlabel('z');
ylabel('Neutron Flux \Phi(z)');
title('Neutron Flux Distribution in a 1D Slab');
grid on;


% figure;
% semilogy(1:it_pow_acc, err_k_vec_acc, 'LineWidth', 2); hold on;
% semilogy(1:it_pow_acc, err_phi_L2_acc, 'LineWidth', 2); 
% semilogy(1:it_pow_acc, err_phi_Linf_acc, 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Relative error');
% title('Convergence of k  and Phi');
% legend('k_{eff}','Phi L_{2}','Phi L_{\infty}')
% grid on;



figure;
semilogy(1:it_pow_acc, err_k_vec_acc, 'LineWidth', 2); hold on;
semilogy(1:it_pow_acc, err_phi_L2_acc, 'LineWidth', 2); 
%semilogy(1:it_pow_acc, err_phi_Linf_acc, 'LineWidth', 2);
semilogy(1:it_pow, err_k_vec, 'LineWidth', 2); hold on;
semilogy(1:it_pow, err_phi_L2, 'LineWidth', 2); 
%semilogy(1:it_pow, err_phi_Linf, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Relative error');
title('Convergence of k  and Phi');
% legend('acc - k_{eff}','acc - Phi L_{2}','acc - Phi L_{\infty}','k_{eff}','Phi L_{2}','Phi L_{\infty}')
legend('acc - k_{eff}','acc - Phi L_{2}','k_{eff}','Phi L_{2}')
grid on;
