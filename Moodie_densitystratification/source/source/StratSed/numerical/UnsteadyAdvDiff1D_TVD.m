function [T_new] = UnsteadyAdvDiff1D_TVD(scheme,beta_TVD,N,x,dx,dt,f,u,M,S,T_old,T_new,...
							BC_low_a_old,BC_low_b_old,BC_low_c_old,BC_up_a_old,BC_up_b_old,BC_up_c_old, ...
							BC_low_a_new,BC_low_b_new,BC_low_c_new,BC_up_a_new,BC_up_b_new,BC_up_c_new)
    %%%%%%%%%%%%%%%%%%%%%% No1Seed Toolbox %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %																				%
    %   "UnsteadyAdvDiff1D" integrates the adv-diffusion equation					%
    %				        dT/dt+d(uT)/dx = d(MdT/dx)/dx + S						%
    %						over one time step dt									%
    %	- scheme 0=UD 1=CD 2=LUD 3=QUICK 4=van Leer 5=van Albada					%
    %			 6= Min-Mod 7=SUPERBEE 8=Sweby 9=QUICK(Mod) 10=UMIST				%
    %   - u = u(x) is advection velocity											%
    %   - M = M(x) is conductivity													%
    %   - S = S(x) is source term													%
    %	- T_new = T_new(x) is main variable (input/output)							%
    %	- u,M,S are piecewise constant over dt										%				
    %   - f = implicit ratio  1=fully implicit										%
    %						  0=fully explicit										% 			
    %                         0.5=Crank-Nicolson									%
    %   - u_mid = u at midpoint, using arithmetic mean (linear)						%
    %   - M_mid = M at midpoint, using harmonic mean								%
    %   - Boundary condition BC_a(t)*q+BC_b(t)*T=BC_c(t)							%
    %   - q = -M(dT/dx) is diffusive flux term										%
    %	- Boundary condition type must be consistent although						%
    %						 may change in value over time							%													
    %   - Non-uniform grid size dx													%
    %   - Remarks: 1. This solver may NOT work when Robin and/or					%
    %				  Neumann conditions are specified at both   					%
    %			      boundaries													%
    %			   2. Can work for both u=0 (pure diffusion)						%
    %									M=0 (pure advection)						%
    %			   3. When scheme=8 (Sweby), beta_TVD needs to be specified.		%
    %				  1<=beta_TVD<=2 guarantees second order %						%
    %			   4. Schemes 1,2,3 are not TVD schemes								%
    %			   5. T_new as an input gives the initial guess to compute the		%
    %				  deferred correction terms. Problem will occur if it is set to %
    %				  be uniform. To avoid this problem, can perform a UD			%
    %				  computation prior to other schemes.							%
    %																				%
    %											(Yeh 1-21-2011)						%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % my preallocations
    tol = 1d-10;
    psi_TVD = zeros(1,N-1);
    Matrix_B = zeros(1,N);
    Matrix_C = zeros(1,N);
    Matrix_D = zeros(1,N);
    Matrix_E = zeros(1,N);
    dV = zeros(1,N);

    u_mid(1:N-1) = (u(1:N-1) + u(2:N)) ./ 2d0;
    M_mid(1:N-1) = (M(1:N-1) + M(2:N)) ./ 2d0;

    dV(2:N-1) = (dx(1:N-2) + dx(2:N-1)) ./ 2d0;
    dV(1) = dx(1) ./ 2d0;
    dV(N) = dx(N-1) ./ 2d0;

    % compute a_W,a_P,a_E %
    for i = 2:(N-1) % i indicate the actual grid point %
        a_W(i) = M_mid(i-1) / dx(i-1) + max([u_mid(i-1), 0d0]);
        a_E(i) = M_mid(i) / dx(i) + max([-u_mid(i), 0d0]);
        a_P(i) = a_W(i) + a_E(i) + (u_mid(i) - u_mid(i-1));
    end
    a_W(1) = 0d0;
    a_W(N) = M_mid(N-1) / dx(N-1) + max([u_mid(N-1), 0d0]); % Grid 1 doesn't have a_W
    a_E(1) = M_mid(1) / dx(1) + max([-u_mid(1), 0d0]); % Grid N doesn't have a_E
    a_E(N) = 0d0;

    % compute g_old %
	T(2:N+1) = T_old;
	T(1) = 2d0*T(2)-T(3);
	T(N+2) = 2d0*T(N+1)-T(N);

	for i = 1:(N-1)
		if (u_mid(i) > 0d0)
			r(i) = (T(i+1)-T(i)) ./ (T(i+2)-T(i+1));
		else
			r(i) = (T(i+2)-T(i+3)) ./ (T(i+1)-T(i+2));
        end
    end

	switch scheme
        case 0	% UD
			psi_TVD(:) = 0d0;
        case 1		% CD
			psi_TVD(:) = 1d0;
        case 2		% LUD
			psi_TVD(:) = r;
        case 3		% QUICK
			psi_TVD(:) = (3d0+r) / 4d0;
        case 4		% van Leer
			psi_TVD(:) = (r+dabs(r)) / (1d0+r);
        case 5	% van Albada
			psi_TVD(:) = (r+r.^2d0)/(1d0+r.^2d0);
        case 6		% Min-Mod
			for i = 1:(N-1)
				if (r(i) > 0d0)
					psi_TVD(i) = min([r(i), 1d0]);
				else
					psi_TVD(i) = 0d0;
                end
            end
        case 7		% SUPERBEE
			for i = 1:(N-1)
				psi_TVD(i) = max( [0d0, min([2d0*r(i), 1d0]), min([r(i), 2d0])] );
            end
        case 8		% Sweby
			for i = 1:(N-1)
				psi_TVD(i) = max( [0d0, min([beta_TVD*r(i), 1d0]), min([r(i), beta_TVD])] );
            end
        case 9		% QUICK(Mod)
			for i = 1:(N-1)
				psi_TVD(i) = max( [0d0, min([2d0*r(i), (3d0+r(i)) / 4d0, 2d0])] );
            end
        case 10	%UMIST
			for i = 1:(N-1)
				psi_TVD(i) = max( [0d0, min([2d0*r(i), (1d0+3d0*r(i)) / 4d0, (3d0+r(i)) / 4d0, 2d0])] );
            end
    end

    % compute secondary terms (?) %
	for i = 1:(N-1)
		if (u_mid(i) > 0d0)
			g(i) = 0.5d0 * psi_TVD(i) * (T(i+2)-T(i+1));
		else
			g(i) = 0.5d0 * psi_TVD(i) * (T(i+1)-T(i+2));
        end
		
		if isnan(g(i))
			g(i) = 0d0; % turn into UD scheme
        end
    end

	Matrix_B_inv(2:N-1) = (1d0-f) * (a_W(2:N-1) .* T(2:N-1) - a_P(2:N-1) .* T(3:N) + a_E(2:N-1) .* T(4:N+1)) + ...
                          (1d0-f) * (-u_mid(2:N-1) .* g(2:N-1) + u_mid(1:N-2) .* g(1:N-2)) + ...
						  dV(2:N-1) ./ dt .* T(3:N) + S(2:N-1) .* dV(2:N-1); % the invariant part of B (exclude deferred correction)
	Matrix_C(2:N-1) = -f * a_W(2:N-1);
	Matrix_D(2:N-1) = f * a_P(2:N-1) + dV(2:N-1) / dt;
	Matrix_E(2:N-1) = -f * a_E(2:N-1);


	if (abs(BC_low_a_new) < 10d-12)	% Dirichlet %
		Matrix_B_inv(1) = BC_low_c_new / BC_low_b_new * Matrix_D(2);
		Matrix_C(1) = 0d0;
		Matrix_D(1) = 1d0*Matrix_D(2);
		Matrix_E(1) = 0d0;
    else % Neumann or Robin %
        Matrix_B_inv(1) = (1d0-f) * ((-(a_E(1) + u_mid(1)) - BC_low_b_old / BC_low_a_old+u(1)) * T(2) + a_E(1) * T(3)) + ...
						 (1d0-f) * (-u_mid(1) * g(1)) + ...
						 dV(1) / dt * T(2) + f * BC_low_c_new / BC_low_a_new + (1d0 - f) * BC_low_c_old / BC_low_a_old + ...
					     S(1) * dV(1);
		Matrix_C(1) = 0d0;
		Matrix_D(1) = f*(a_E(1)+u_mid(1)+BC_low_b_new/BC_low_a_new-u(1))+dV(1)/dt;
		Matrix_E(1) = -f*a_E(1);
    end

	if (abs(BC_up_a_new) < 10d-12) % Dirichlet %
		Matrix_B_inv(N) = BC_up_c_new / BC_up_b_new * Matrix_D(N-1);
		Matrix_C(N) = 0d0;
		Matrix_D(N) = 1d0 * Matrix_D(N-1);
		Matrix_E(N) = 0d0;
    else % Neumann or Robin %
		Matrix_B_inv(N) = (1d0-f) * ((-(a_W(N) - u_mid(N-1)) + BC_up_b_old / BC_up_a_old - u(N)) * T(N+1) + a_W(N) * T(N)) + ...
						  (1d0-f) * (u_mid(N-1) * g(N-1)) + ...
						  dV(N) / dt * T(N+1) - f * BC_up_c_new / BC_up_a_new - (1d0 - f) * BC_up_c_old / BC_up_a_old + ...
						  S(N) * dV(N);
		Matrix_C(N) = -f * a_W(N);
		Matrix_D(N) = f * (a_W(N) - u_mid(N-1) - BC_up_b_new / BC_up_a_new + u(N)) + dV(N) / dt;
		Matrix_E(N) = 0d0;
    end

    % iteration %
	err = 10d0;
	while (err > tol)
        T(2:N+1) = T_new;
		T(1) = 2d0*T(2)-T(3);
		T(N+2) = 2d0*T(N+1)-T(N);

        for i = 1:(N-1)
			if (u_mid(i) > 0d0)
				r(i) = (T(i+1)-T(i)) ./ (T(i+2)-T(i+1));
			else
				r(i) = (T(i+2)-T(i+3)) ./ (T(i+1)-T(i+2));
            end
        end


		switch (scheme)
            case (0)		% UD
                psi_TVD(:) = 0d0;
            case (1)		% CD
                psi_TVD(:) = 1d0;
            case (2)		% LUD
                psi_TVD(:) = r;
            case (3)		% QUICK
                psi_TVD(:) = (3d0+r) / 4d0;
            case (4)		% van Leer
                psi_TVD(:) = (r+dabs(r)) / (1d0+r);
            case (5)		% van Albada
                psi_TVD(:) = (r+r.^2d0) / (1d0+r.^2d0);
            case (6)		% Min-Mod
                for i = 1:(N-1)
                    if (r(i) > 0d0)
                        psi_TVD(i) = min([r(i), 1d0]);
                    else
                        psi_TVD(i) = 0d0;
                    end
                end
            case (7)		% SUPERBEE
                for i = 1:(N-1)
                    psi_TVD(i) = max( [0d0, min([2d0*r(i), 1d0]), min([r(i), 2d0])] );
                end
            case (8)		% Sweby
                for i = 1:(N-1)
                    psi_TVD(i) = max( [0d0, min([beta_TVD*r(i), 1d0]), min([r(i), beta_TVD])] );
                end
            case (9)		% QUICK(Mod)
                for i = 1:(N-1)
                    psi_TVD(i) = max( [0d0, min([2d0*r(i), (3d0+r(i)) / 4d0, 2d0])] );
                end
            case (10)	%UMIST
                for i = 1:(N-1)
                    psi_TVD(i) = max( [0d0, min([2d0*r(i), (1d0+3d0*r(i)) / 4d0, (3d0+r(i)) / 4d0, 2d0])] );
                end
        end

		for i = 1:(N-1)
			if (u_mid(i) > 0d0)
                g(i) = 0.5d0 * psi_TVD(i) * (T(i+2)-T(i+1));                
			else
                g(i) = 0.5d0 * psi_TVD(i) * (T(i+1)-T(i+2));
            end
			
			if isnan(g(i))
				g(i) = 0d0; % turn into UD scheme
            end
        end

		Matrix_B(2:N-1) = Matrix_B_inv(2:N-1) + f * (-u_mid(2:N-1) .* g(2:N-1) + u_mid(1:N-2) .* g(1:N-2));

		if (abs(BC_low_a_new) < 10d-12) % Dirichlet %
			Matrix_B(1) = Matrix_B_inv(1);
        else % Neumann or Robin %
			Matrix_B(1) = Matrix_B_inv(1) - f * u_mid(1) * g(1);
        end

		if (abs(BC_up_a_new) < 10d-12) % Dirichlet %
			Matrix_B(N) = Matrix_B_inv(N);
        else % Neumann or Robin %
			Matrix_B(N) = Matrix_B_inv(N) + f * u_mid(N-1) * g(N-1);
        end

            
            
		if (abs(f) < 10d-12) % Fully explicit %
			T_new = Matrix_B / Matrix_D; % should be elementwise division?
        else % Implicit %
			Matrix_C_temp = Matrix_C;
			Matrix_D_temp = Matrix_D;
			Matrix_E_temp = Matrix_E;
            
            %% thomas soln
            new_Matrix_B = tridiag( Matrix_D, Matrix_C, Matrix_E, Matrix_B );
    
            T_new = new_Matrix_B';
            
        end

		err = max(abs((T(2:N+1)-T_new)./T_new));
    end

    

end
