function [T_new] = SteadyAdvDiff1D_TVD(scheme,beta_TVD,N,x,dx,u,M,S,T_new,...
                                       BC_low_a,BC_low_b,BC_low_c,BC_up_a,BC_up_b,BC_up_c)
    %%%%%%%%%%%%%%%%%%%%%% No1Seed Toolbox %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %																				%
    %   "SteadyAdvDiff1D_TVD" solves the steady state advection-diffusion equation 	%
    %						  d(uT)/dx = d(MdT/dx)/dx + S							%
    %                         using TVD schemes										%
    %	- scheme 0=UD 1=CD 2=LUD 3=QUICK 4=van Leer 5=van Albada					%
    %			 6= Min-Mod 7=SUPERBEE 8=Sweby 9=QUICK(Mod) 10=UMIST				%
    %   - u = u(x) is advection velocity											%
    %   - M = M(x) is conductivity													%
    %   - S = S(x) is source term													%
    %	- T_new = T_new(x) is main variable (input/output)							%
    %   - u_mid = u at midpoint, using arithmetic mean (linear)						%
    %   - M_mid = M at midpoint, using arithmetic mean (linear)						%
    %   - Boundary condition BC_a*q+BC_b*T=BC_c										%
    %   - q = -M(dT/dx) is diffusive flux term										%
    %   - Non-uniform grid size 													%
	%   - Remarks: 1. This solver may NOT work when Robin and/or Neumann conditions	%
    %				  are specified at both boundaries								%
    %			   2. Can work for both u=0	(pure diffusion)						%
    %								and M=0 (pure advection)						%
    %			   3. When scheme=8 (Sweby), beta_TVD needs to be specified.		%
    %				  1<=beta_TVD<=2 guarantees second order %						%
    %			   4. Schemes 1,2,3 are not TVD schemes								%
    %			   5. T_new as an input gives the initial guess to compute the		%
    %				  deferred correction terms. Problem will occur if it is set to %
    %				  be uniform. To avoid this problem, can perform a UD			%
    %				  computation prior to other schemes.							%
    %																				%
    %															(Yeh 1-14-2011)		%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % my preallocations
    tol = 1d-10;
    psi_TVD = zeros(1,N-1);
    Matrix_B = zeros(1,N);
    Matrix_C = zeros(1,N);
    Matrix_D = zeros(1,N);
    Matrix_E = zeros(1,N);
    T = zeros(1,N+2);

    u_mid(1:N-1) = (u(1:N-1)+u(2:N)) ./ 2d0;
    M_mid(1:N-1) = (M(1:N-1)+M(2:N)) ./ 2d0;
    
    %%% Iteration %%%

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
                psi_TVD(:) = (3d0+r)/4d0;
            case (4)		% van Leer
                psi_TVD(:) = (r+dabs(r))/(1d0+r);
            case (5)		% van Albada
                psi_TVD(:) = (r+r.^2d0)/(1d0+r.^2d0);
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
                    psi_TVD(i) = max( [0d0, min([2d0*r(i), (3d0+r(i))/4d0, 2d0])] );
                end
            case (10)	%UMIST
                for i = 1:(N-1)
                    psi_TVD(i) = max( [0d0, min([2d0*r(i), (1d0+3d0*r(i))/4d0, (3d0+r(i))/4d0, 2d0])] );
                end
        end

    	% compute secondary terms %
        for i = 1:(N-1)
    		if (u_mid(i) >  0d0)
    			g(i) = 0.5d0 * psi_TVD(i) * ( T(i+2)-T(i+1) );
    		else
    			g(i) = 0.5d0 * psi_TVD(i) * ( T(i+1)-T(i+2) );
            end
    
            if isnan(g(i))
                g(i) = 0d0; % turn into UD scheme
            end
        end
    
    	% formulate matrix %
    	for i = 2:(N-1)
    		Matrix_B(i) = S(i) * (dx(i-1) + dx(i)) / 2d0 - u_mid(i) * g(i) + u_mid(i-1) * g(i-1);
    		Matrix_C(i) = -( M_mid(i-1)/dx(i-1) + max([u_mid(i-1), 0d0]) ); % this = -a_W
    		Matrix_E(i) = -( M_mid(i)/dx(i) + max([-u_mid(i), 0d0]) ); % this = -a_E
    		Matrix_D(i) = -Matrix_C(i) - Matrix_E(i) + (u_mid(i) - u_mid(i-1)); % this = a_P
        end
    
    	if (abs(BC_low_a) < 10d-12) % Dirichlet %
            Matrix_C(1) = 0d0;
            Matrix_D(1) = 1d0*Matrix_D(2);
            Matrix_E(1) = 0d0;
            Matrix_B(1) = BC_low_c / BC_low_b * Matrix_D(2);
        else % Neumann or Robin %
            Matrix_C(1) = 0d0;				
            Matrix_E(1) = -( M_mid(1) / dx(1) + max([-u_mid(1), 0d0]) ); % this = -a_E
            Matrix_D(1) = -Matrix_E(1) + (u_mid(1) - u(1)) + BC_low_b/BC_low_a; % this = a_P
            Matrix_B(1) = BC_low_c / BC_low_a + S(1) * dx(1) / 2d0 - u_mid(1) * g(1);
        end
    
    	if (abs(BC_up_a) < 10d-12) % Dirichlet %
            Matrix_C(N) = 0d0;
            Matrix_D(N) = 1d0 * Matrix_D(N-1);
            Matrix_E(N) = 0d0;
            Matrix_B(N) = BC_up_c/BC_up_b * Matrix_D(N-1);
        else % Neumann or Robin %
            Matrix_C(N) = -(M_mid(N-1) / dx(N-1) + max([u_mid(N-1), 0d0])); % this = -a_W
            Matrix_E(N) = 0d0;
            Matrix_D(N) = -Matrix_C(N) + u(N) - u_mid(N-1) - BC_up_b / BC_up_a;	% this = a_P
            Matrix_B(N) = -BC_up_c / BC_up_a + S(N) * dx(N-1) / 2d0 + u_mid(N-1) * g(N-1);
        end
        
        %% THOMAS ALGO
        new_Matrix_B = tridiag( Matrix_D, Matrix_C, Matrix_E, Matrix_B );
    
    	T_new = new_Matrix_B';
    	err = max(abs((T(2:N+1)-T_new)./T_new));
    end % end while loop



end









