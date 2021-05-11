function [T_new] = UnsteadyAdvDiff1D_consC(N,x,dx,dt,f,u,M,S,T_old,T_new,int_T,...
                            BC_low_a_old,BC_low_b_old,BC_low_c_old,BC_up_a_old,BC_up_b_old,BC_up_c_old,...
                            BC_low_a_new,BC_low_b_new,BC_low_c_new,BC_up_a_new,BC_up_b_new,BC_up_c_new)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                   %
    %   "UnsteadyAdvDiff1D_consC" is modified from "UnsteadyAdvDiff1D" %
    %                              due to the dual Robin condition      %
    %                              applied at both ends                   %
    %    - last condition is specified by the conservation of           %
    %      concentration, using Simpson's integration                   %
    %    - solved using sparse matrix solver instead of band solver     %
    %                                                                   %
    %                                              Yeh (3-18-2010)       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Matrix_B = zeros(1,N);
    Matrix_C = zeros(1,N);
    Matrix_D = zeros(1,N);
    Matrix_E = zeros(1,N);
    dV = zeros(1,N);
    S = S .* ones(1,N); % copy S into S elelemnts?

    M_mid(1:N-1) = (M(1:N-1) + M(2:N)) / 2d0; % 2d0*M(1:N-1)*M(2:N)/(M(1:N-1)+M(2:N)) %
    u_mid(1:N-1) = (u(1:N-1) + u(2:N)) / 2d0;
    Pe_mid(1:N-1) = u_mid(1:N-1) * dx / M_mid(1:N-1);
    
    for i = 2:(N-1) % i indicate the actual grid point %
       a_W(i) = M_mid(i-1) / dx * max([0d0, (1d0 - 0.1d0 * abs(Pe_mid(i-1))).^5]) + max([u_mid(i-1), 0d0]);
       a_E(i) = M_mid(i) / dx * max([0d0, (1d0 - 0.1d0 * abs(Pe_mid(i))).^5]) + max([-u_mid(i), 0d0]);
       a_P(i) = a_W(i) + a_E(i) + (u_mid(i) - u_mid(i-1));
    end
    
    a_W(N) = M_mid(N-1) / dx * max([0d0, (1d0 - 0.1d0 * abs(Pe_mid(N-1))).^5]) + max([u_mid(N-1), 0d0]);
    a_W(1) = 0d0; % Grid 1 doesn't have a_W
    a_E(N) = 0d0; % Grid N doesn't have a_E
    a_E(1) = M_mid(1) / dx * max([0d0, (1d0 - 0.1d0 * abs(Pe_mid(1))).^5]) + max([-u_mid(1), 0d0]);

    Matrix_B(2:N-1) = (1d0-f) * (a_W(2:N-1) .* T_old(1:N-2) - a_P(2:N-1) .* T_old(2:N-1) + a_E(2:N-1) .* T_old(3:N)) + dx / dt * T_old(2:N-1) + S(2:N-1) * dx;
    Matrix_C(2:N-1) = -f * a_W(2:N-1);
    Matrix_D(2:N-1) = f * a_P(2:N-1) + dx / dt;
    Matrix_E(2:N-1) = -f * a_E(2:N-1);

    
    
    if (abs(BC_low_a_new) < 10d-12) % Dirichlet %
        Matrix_B(1) = BC_low_c_new / BC_low_b_new * Matrix_D(2);
        Matrix_C(1) = 0d0;
        Matrix_D(1) = 1d0 * Matrix_D(2);
        Matrix_E(1) = 0d0;
    else % Neumann or Robin %
        Matrix_B(1) = (1d0-f) * ((-(a_E(1) + u_mid(1)) - BC_low_b_old / BC_low_a_old + u(1)) * T_old(1) + a_E(1) * T_old(2)) + dx / dt / 2d0 * T_old(1) + ...
                      f * BC_low_c_new / BC_low_a_new + (1d0-f) * BC_low_c_old / BC_low_a_old + ...
                      S(1) * dx / 2d0;
        Matrix_C(1) = 0d0;
        Matrix_D(1) = f * (a_E(1) + u_mid(1) + BC_low_b_new / BC_low_a_new - u(1)) + dx / dt / 2d0;
        Matrix_E(1) = -f * a_E(1);
    end

    if (abs(BC_up_a_new) < 10d-12) % Dirichlet %
        Matrix_B(N) = BC_up_c_new / BC_up_b_new * Matrix_D(N-1);
        Matrix_C(N) = 0d0;
        Matrix_D(N) = 1d0 * Matrix_D(N-1);
        Matrix_E(N) = 0d0;
    else % Neumann or Robin %
        Matrix_B(N) = (1d0-f) * ((-(a_W(N) - u_mid(N-1)) + BC_up_b_old / BC_up_a_old - u(N)) * T_old(N) + a_W(N) * T_old(N-1)) + dx / dt / 2d0 * T_old(N) - ...
                      f * BC_up_c_new / BC_up_a_new - (1d0-f) * BC_up_c_old / BC_up_a_old + ...
                      S(N) * dx / 2d0;
        Matrix_C(N) = -f * a_W(N);
        Matrix_D(N) = f * (a_W(N) - u_mid(N-1) - BC_up_b_new / BC_up_a_new + u(N)) + dx / dt / 2d0;
        Matrix_E(N) = 0d0;
    end
    
    % preallocate the nonzero things
    nonzero_row = zeros(1, 4*N-4);
    nonzero_col = zeros(1, 4*N-4);
    nonzero_A = zeros(1, 4*N-4);
    
    nonzero_row(1) = 1;
    nonzero_row(2) = 1;
    nonzero_col(1) = 1;
    nonzero_col(2) = 2;
    nonzero_A(1) = Matrix_D(1);
    nonzero_A(2) = Matrix_E(1);
    
    % doesn't the following line just overwrite the boundary condition set above???
    Matrix_B(N) = int_T / dx * 3d0;
    
    for i = 2:(N-1)
        for j = (i-1):(i+1)
           nonzero_row(3*i-3+(j-(i-1))) = i;
           nonzero_col(3*i-3+(j-(i-1))) = j;
        end
        nonzero_A(3*i-3) = Matrix_C(i);
        nonzero_A(3*i-2) = Matrix_D(i);
        nonzero_A(3*i-1) = Matrix_E(i);
    end


    for j = 1:N
        nonzero_row(3*N-4+j) = N;
        nonzero_col(3*N-4+j) = j;
        if ( (j == 1) || (j == N) )
            nonzero_A(3*N-4+j) = 1d0;
        elseif ( mod(j,2) == 1 )
            nonzero_A(3*N-4+j) = 2d0;
        elseif ( mod(j,2) == 0)
            nonzero_A(3*N-4+j) = 4d0;
        end
    end

    Matrix_S = sparse(nonzero_row, nonzero_col, nonzero_A, N, N);
    new_Matrix_B = Matrix_S \ Matrix_B';
    T_new = new_Matrix_B';

    
end
