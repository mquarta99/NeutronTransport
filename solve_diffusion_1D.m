function Phi = solve_diffusion_1D(Nz, dz, D, sigma_a, source)
    A = zeros(Nz, Nz);
    b = source * dz^2;

    for i = 1:Nz
        if i == 1 % null flux b.c.
            A(i, i) = 1;
            b(i) = 0;
        elseif i == Nz % null flux b.c.
            A(i, i) = 1;
            b(i) = 0;
        % if i == 1  % Robin b.c.
        %     A(i, i)   = 1 + 2*D/dz;
        %     A(i, i+1) = -2*D/dz;
        %     b(i) = 0;
        % elseif i == Nz  % Robin b.c.
        %     A(i, i)   = 1 + 2*D/dz;
        %     A(i, i-1) = -2*D/dz;
        %     b(i) = 0;
        else
            A(i, i-1) = -D;
            A(i, i)   = 2*D + sigma_a*dz^2;
            A(i, i+1) = -D;
        end
    end
    Phi = A \ b;
end