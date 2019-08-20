function [eta] = calculate_new_bed_exner(eta, au, qs, qs_u, Bc, Be, nx, dx, dt, Int, con)
	% applies exner equation for mass conservation, evolving the best for dt
	
    eta0 = eta;
    Qs = qs .* Bc; % changed this to only account for the width of the channel for Q calc, but use Be for exner calc
    Qsu = qs_u * Bc(1);
    
    dQsdx(nx+1) = (Qs(nx+1) - Qs(nx)) / dx;
    dQsdx(1) = au*(Qs(1)-Qsu)/dx + (1-au)*(Qs(2)-Qs(1))/dx;
    dQsdx(2:nx) = au*(Qs(2:nx)-Qs(1:nx-1))/dx + (1-au)*(Qs(3:nx+1)-Qs(2:nx))/dx; % use winding coefficient in central portion
    
    eta = eta0 - ((1/(1-con.phi)) .* dQsdx .* dt .* (1 ./ Be) .* Int);

end